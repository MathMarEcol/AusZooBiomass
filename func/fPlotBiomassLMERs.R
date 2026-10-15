# fPlotBiomassLMERs.R
# Covariate effect plots for lm / lmer / glmer models
# Uses marginaleffects + ggplot2 + patchwork
#
# Author: Jason Everett
# Updated: October 2026
#
# Dependencies:
#   marginaleffects, ggplot2, patchwork, insight
#   (all available on CRAN)
#
# Usage:
#   p <- fPlotBiomassLMERs(model, covariates = list(SST = "SST (°C)", ...))
#   ggsave("Figures/MyModel.pdf", p, width = 420, height = 297, units = "mm", device = cairo_pdf)

library(marginaleffects)
library(ggplot2)
library(patchwork)
library(rlang) # `:=` walrus operator for tidy evaluation in datagrid()


# ── Helper: auto-detect interaction pairs from model formula ──────────────────
#
# Parses the fixed-effect formula and returns a list of 2-element character
# vectors, one per interaction term (terms joined by * or :).
#
# Strategy:
#   1. Extract all terms from the formula via terms().
#   2. Keep only those whose label contains ":" (interaction indicator used
#      by R regardless of whether the user wrote * or :).
#   3. Split each label on ":" to get the two (or more) component term strings.
#   4. For each component, extract the *first* bare symbol that looks like a
#      variable name — i.e. the first word-token that is not a function name
#      (we heuristically treat the outermost call's first argument as the
#      variable).  This handles wrappers like harmonic(HarmTOD, k = 1),
#      ns(Bathy, df = 2), exp(-Depth / 1000), log10(Chl), etc.
#   5. Return only pairs (length-2 interactions); higher-order interactions
#      are skipped with a warning.
#
# Returns NULL if no interactions are found.

.extract_interactions <- function(model) {
  fml <- stats::formula(model)

  # Drop the response and random-effect parts; keep only fixed-effect RHS.
  # nobars() strips random-effect bars; it moved from lme4 to reformulas.
  .nobars <- if (requireNamespace("reformulas", quietly = TRUE)) {
    reformulas::nobars
  } else {
    lme4::nobars
  }
  rhs <- if (is(model, "merMod")) {
    .nobars(fml)[-2] # remove response (LHS)
  } else {
    stats::delete.response(stats::terms(fml))
  }

  trms <- stats::terms(rhs)
  labels <- attr(trms, "term.labels") # e.g. "harmonic(HarmTOD, k = 1):exp(-Depth/1000)"
  int_labels <- labels[grepl(":", labels, fixed = TRUE)]

  if (length(int_labels) == 0L) {
    return(NULL)
  }

  # Extract the raw variable name from a term string such as:
  #   "harmonic(HarmTOD, k = 1)"  → "HarmTOD"
  #   "exp(-Depth / 1000)"        → "Depth"
  #   "ns(Bathy, df = 2)"         → "Bathy"
  #   "log10(Chl)"                → "Chl"
  #   "SST"                       → "SST"   (bare variable)
  .raw_var <- function(term_str) {
    term_str <- trimws(term_str)
    # If it contains a '(' it is a function call — parse and walk the AST
    if (grepl("(", term_str, fixed = TRUE)) {
      expr <- tryCatch(str2lang(term_str), error = function(e) NULL)
      if (is.null(expr)) {
        return(term_str)
      }
      # Recursively find the first SYMBOL that is not a function name
      .first_sym <- function(node) {
        if (is.symbol(node)) {
          return(as.character(node))
        }
        if (is.call(node)) {
          # node[[1]] is the function name; node[[2]], node[[3]], ... are args
          for (i in seq(2L, length(node))) {
            hit <- .first_sym(node[[i]])
            if (!is.null(hit)) return(hit)
          }
        }
        NULL
      }
      hit <- .first_sym(expr)
      return(if (!is.null(hit)) hit else term_str)
    }
    term_str # bare variable name
  }

  pairs <- lapply(int_labels, function(lbl) {
    parts <- strsplit(lbl, ":", fixed = TRUE)[[1]]
    if (length(parts) != 2L) {
      warning(
        "Skipping higher-order interaction (>2 terms): ",
        lbl,
        call. = FALSE
      )
      return(NULL)
    }
    unname(vapply(parts, .raw_var, character(1L)))
  })

  Filter(Negate(is.null), pairs)
}


# ── Helper: detect model family and link ──────────────────────────────────────
#
# Returns a list with:
#   $is_glm  — TRUE if the model is a glm/glmer (non-Gaussian or non-identity link)
#   $link    — character link name, e.g. "log", "logit", "identity"
#   $family  — character family name, e.g. "Gamma", "gaussian"
#
# For lm / lmer (Gaussian identity), is_glm = FALSE and link = "identity".
# family() throws an error for lm/lmer, which we catch and treat as Gaussian.

.model_family_info <- function(model) {
  fam <- tryCatch(family(model), error = function(e) NULL)
  if (is.null(fam)) {
    return(list(is_glm = FALSE, link = "identity", family = "gaussian"))
  }
  list(
    is_glm = !(fam$family == "gaussian" && fam$link == "identity"),
    link = fam$link,
    family = fam$family
  )
}


# ── Helper: compute partial residuals for a focal variable ────────────────────
#
# Partial residual for variable `var`:
#   e_partial = residuals(model) + X_var %*% beta_var
#
# For lm / lmer the residuals are on the link (= response) scale, so the
# partial residuals are also on the link scale.  Back-transformation to the
# response scale (10^x) is handled in .plot_one_covariate() when
# response_scale = TRUE.
#
# For glmer with a non-identity link, partial residuals are suppressed entirely
# because the working-residual decomposition is on the link scale and is
# misleading when plotted on the response scale.  The caller checks is_glm
# before calling this function.

.partial_residuals <- function(model, var) {
  mm <- model.matrix(model) # full fixed-effect design matrix
  cf <- lme4::fixef(model) # named fixed-effect coefficients
  # fall back for plain lm objects
  if (inherits(model, "lm") && !inherits(model, "lmerMod")) {
    cf <- stats::coef(model)
  }

  # Identify columns of the model matrix that belong to this variable.
  # Column names for spline / harmonic terms look like:
  #   "ns(Bathy, df = 2)1", "harmonic(HarmDOY, k = 1)sin1", etc.
  # We match any column whose name contains the raw variable name.
  col_idx <- grep(var, colnames(mm), fixed = TRUE)

  if (length(col_idx) == 0L) {
    warning(
      "No model-matrix columns found for variable '",
      var,
      "'. Partial residuals will not be shown."
    )
    return(NULL)
  }

  partial_fit <- as.numeric(mm[, col_idx, drop = FALSE] %*% cf[col_idx])
  resids <- stats::residuals(model)

  partial_fit + resids
}


# ── Helper: axis labels for harmonic variables ─────────────────────────────────

.harmonic_breaks_tod <- function() {
  list(
    breaks = seq(0, 2 * pi, by = pi / 2),
    labels = c("00:00", "06:00", "12:00", "18:00", "24:00")
  )
}

.harmonic_breaks_doy <- function() {
  # Month midpoints as fraction of year × 2π
  mid_days <- c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349)
  list(
    breaks = (mid_days / 365) * 2 * pi,
    labels = month.abb
  )
}


# ── Helper: build a single covariate panel ────────────────────────────────────

.plot_one_covariate <- function(
  model,
  var,
  xlab,
  re.form,
  response_scale,
  show_rug,
  show_partial_residuals,
  fam_info
) {
  mf <- insight::get_data(model)
  is_factor <- is.factor(mf[[var]]) || is.character(mf[[var]])

  if (is_factor) {
    nd <- rlang::exec(
      datagrid,
      model = model,
      !!!setNames(list(unique(mf[[var]])), var)
    )
  } else {
    nd <- rlang::exec(
      datagrid,
      model = model,
      !!!setNames(
        list(seq(
          min(mf[[var]], na.rm = TRUE),
          max(mf[[var]], na.rm = TRUE),
          length.out = 200
        )),
        var
      )
    )
  }

  # type = "response" behaviour differs by model class:
  #   lmer/lm  — returns predictions on the model's response scale
  #              (e.g. log10(Biomass) if the response was log10-transformed).
  #              Apply 10^ manually if response_scale = TRUE.
  #   glmer    — marginaleffects back-transforms through the link automatically,
  #              returning predictions in mg/m³.  Never apply 10^ manually.
  preds <- predictions(
    model,
    newdata = nd,
    re.form = re.form,
    type = "response"
  )

  # Issue 2 fix: only apply 10^ back-transform for lmer/lm, not glmer
  if (!fam_info$is_glm && response_scale) {
    preds$estimate <- 10^preds$estimate
    preds$conf.low <- 10^preds$conf.low
    preds$conf.high <- 10^preds$conf.high
  }

  # Issue 4 fix: y-axis label depends on model type and response_scale
  #   glmer (any response_scale) → always on response scale → "Biomass (mg/m³)"
  #   lmer + response_scale TRUE → back-transformed          → "Biomass (mg/m³)"
  #   lmer + response_scale FALSE → still on log10 scale     → "log10(Biomass)"
  ylab <- if (fam_info$is_glm || response_scale) {
    expression("Biomass (mg/m"^3 * ")")
  } else {
    expression(log[10] * "(Biomass)")
  }

  if (is_factor) {
    p <- ggplot(
      preds,
      aes(
        x = .data[[var]],
        y = estimate,
        ymin = conf.low,
        ymax = conf.high
      )
    ) +
      geom_pointrange(colour = "steelblue", linewidth = 0.7, size = 0.5) +
      labs(x = xlab, y = ylab) +
      theme_bw(base_size = 10) +
      theme(
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.text.y = element_text(),
        axis.ticks = element_line()
      )
  } else {
    p <- ggplot(preds, aes(x = .data[[var]])) +
      geom_ribbon(
        aes(ymin = conf.low, ymax = conf.high),
        fill = "steelblue",
        alpha = 0.25
      ) +
      geom_line(aes(y = estimate), colour = "steelblue", linewidth = 0.9) +
      labs(x = xlab, y = ylab) +
      theme_bw(base_size = 10) +
      theme(
        axis.text.x = element_text(),
        axis.text.y = element_text(),
        axis.ticks = element_line()
      )

    # Issue 5 fix: suppress partial residuals entirely for glmer models.
    # For glmer, the working-residual decomposition is on the link scale and
    # is misleading when plotted on the response scale.
    if (show_partial_residuals && !fam_info$is_glm) {
      pr <- .partial_residuals(model, var)
      if (!is.null(pr)) {
        pr_df <- data.frame(
          x = mf[[var]],
          pr = if (response_scale) 10^pr else pr
        )
        p <- p +
          geom_point(
            data = pr_df,
            aes(x = x, y = pr),
            inherit.aes = FALSE,
            colour = "grey30",
            alpha = 0.25,
            size = 0.8,
            shape = 16
          )
      }
    }

    if (show_rug) {
      p <- p +
        geom_rug(
          data = mf,
          aes(x = .data[[var]]),
          inherit.aes = FALSE,
          colour = "grey40",
          alpha = 0.4,
          length = unit(0.015, "npc"),
          sides = "b"
        )
    }

    if (var == "HarmTOD") {
      brks <- .harmonic_breaks_tod()
      p <- p + scale_x_continuous(breaks = brks$breaks, labels = brks$labels)
    }
    if (var == "HarmDOY") {
      brks <- .harmonic_breaks_doy()
      p <- p +
        scale_x_continuous(
          breaks = brks$breaks,
          labels = brks$labels,
          guide = guide_axis(angle = 45)
        )
    }
    if (var == "Bathy") {
      p <- p + scale_x_reverse()
    }
  }

  p
}


# ── Helper: build an interaction panel ────────────────────────────────────────

.plot_interaction <- function(model, vars, re.form, response_scale, fam_info) {
  stopifnot(length(vars) == 2)
  var1 <- vars[1] # x-axis variable
  var2 <- vars[2] # grouping variable

  mf <- insight::get_data(model)
  is_factor2 <- is.factor(mf[[var2]]) || is.character(mf[[var2]])

  v2_vals <- if (is_factor2) {
    unique(mf[[var2]])
  } else {
    quantile(mf[[var2]], probs = c(0.1, 0.33, 0.67, 0.9), na.rm = TRUE)
  }

  nd <- rlang::exec(
    datagrid,
    model = model,
    !!!setNames(
      list(
        seq(
          min(mf[[var1]], na.rm = TRUE),
          max(mf[[var1]], na.rm = TRUE),
          length.out = 100
        ),
        v2_vals
      ),
      c(var1, var2)
    )
  )

  preds <- predictions(
    model,
    newdata = nd,
    re.form = re.form,
    type = "response"
  )

  # Issue 2 fix: only apply 10^ back-transform for lmer/lm, not glmer
  if (!fam_info$is_glm && response_scale) {
    preds$estimate <- 10^preds$estimate
    preds$conf.low <- 10^preds$conf.low
    preds$conf.high <- 10^preds$conf.high
  }

  # Issue 4 fix: y-axis label
  ylab <- if (fam_info$is_glm || response_scale) {
    expression("Biomass (mg/m"^3 * ")")
  } else {
    expression(log[10] * "(Biomass)")
  }

  if (!is_factor2) {
    preds[[var2]] <- round(preds[[var2]], 1)
  }
  preds[[var2]] <- factor(preds[[var2]])

  p <- ggplot(
    preds,
    aes(
      x = .data[[var1]],
      y = estimate,
      colour = .data[[var2]],
      fill = .data[[var2]]
    )
  ) +
    geom_ribbon(
      aes(ymin = conf.low, ymax = conf.high),
      alpha = 0.15,
      colour = NA
    ) +
    geom_line(linewidth = 0.8) +
    scale_colour_viridis_d(option = "plasma", end = 0.85) +
    scale_fill_viridis_d(option = "plasma", end = 0.85) +
    labs(
      x = var1,
      y = ylab,
      colour = var2,
      fill = var2,
      title = paste0(var1, " \u00d7 ", var2)
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 9, face = "bold")
    )

  if (var1 == "HarmTOD") {
    brks <- .harmonic_breaks_tod()
    p <- p + scale_x_continuous(breaks = brks$breaks, labels = brks$labels)
  }
  if (var1 == "HarmDOY") {
    brks <- .harmonic_breaks_doy()
    p <- p +
      scale_x_continuous(
        breaks = brks$breaks,
        labels = brks$labels,
        guide = guide_axis(angle = 45)
      )
  }

  p
}


# ── Helper: extract random-effect grouping variable names ─────────────────────
#
# Returns a character vector of the grouping variable names from the random-
# effect structure of a merMod (lmer / glmer) model.  Returns character(0) for
# plain lm / glm objects (which have no random effects).
#
# Examples:
#   (1 | Gear)              → "Gear"
#   (1 | DatasetID) + (1 | Gear) → c("DatasetID", "Gear")
#   (Depth | Gear)          → "Gear"   (random slope — grouping var is still Gear)
#
# Strategy: insight::find_random() returns a list with element $random
# containing the grouping variable names.  We use that in preference to
# parsing the formula ourselves, because insight handles edge cases
# (nested effects, crossed effects, etc.) correctly.

.extract_random_grouping_vars <- function(model) {
  if (!is(model, "merMod")) {
    return(character(0L))
  }
  rnd <- tryCatch(
    insight::find_random(model, flatten = TRUE),
    error = function(e) character(0L)
  )
  if (is.null(rnd) || length(rnd) == 0L) {
    return(character(0L))
  }
  unique(rnd)
}


# ── Helper: caterpillar plot of BLUPs for one random-effect grouping variable ──
#
# Extracts the conditional modes (BLUPs) and their standard errors for the
# intercept term of the named grouping variable, then plots them as a
# dot-and-whisker (caterpillar) chart sorted by BLUP value.
#
# For models with a log10-transformed response (lmer, response_scale = TRUE)
# the BLUPs are on the log10 scale; we label the axis accordingly.
# For glmer models the BLUPs are on the link scale (e.g. log); we label
# accordingly.
#
# @param model        A merMod object (lmer or glmer).
# @param grp_var      Character. Name of the grouping variable, e.g. "Gear".
# @param fam_info     List from .model_family_info().
# @param response_scale Logical. Only relevant for lmer; if TRUE the y-axis
#                     label says "BLUP (log₁₀ scale)".
#
# @return A ggplot object.

.plot_random_effect <- function(model, grp_var, fam_info, response_scale) {
  # Extract random-effect conditional modes and SEs via lme4::ranef()
  re_obj <- lme4::ranef(model, condVar = TRUE)

  if (!grp_var %in% names(re_obj)) {
    warning(
      "Random-effect grouping variable '",
      grp_var,
      "' not found in ranef(model). Skipping.",
      call. = FALSE
    )
    return(NULL)
  }

  re_df <- re_obj[[grp_var]]

  # condVar attribute: a 3-D array [n_terms, n_terms, n_levels]
  # We want the intercept term (first row/col) for each level.
  cond_var_arr <- attr(re_df, "postVar") # lme4 < 1.1-30 name
  if (is.null(cond_var_arr)) {
    cond_var_arr <- attr(re_df, "condVar") # lme4 >= 1.1-30 name
  }

  # Build a tidy data frame: level, blup, se, lower, upper
  # Use the "(Intercept)" column if present; fall back to the first column.
  intercept_col <- if ("(Intercept)" %in% colnames(re_df)) {
    "(Intercept)"
  } else {
    colnames(re_df)[1L]
  }

  blup_vals <- re_df[[intercept_col]]
  n_levels <- length(blup_vals)

  # Standard errors from the conditional variance array
  if (!is.null(cond_var_arr) && length(dim(cond_var_arr)) == 3L) {
    # Intercept is always the first term (index 1)
    se_vals <- sqrt(cond_var_arr[1L, 1L, ])
  } else {
    # Fallback: no SEs available — plot points only
    se_vals <- rep(NA_real_, n_levels)
  }

  # ── Assign level numbers ───────────────────────────────────────────────────
  # Determine the canonical ordering of levels:
  #   - factor  → use existing factor level order
  #   - character → sort alphabetically
  # Then assign integer indices 1, 2, 3, … in that canonical order.
  # The caterpillar is still sorted by BLUP value on the y-axis, but the
  # tick labels show the level number rather than the (potentially long) name.
  mf <- insight::get_data(model)
  grp_col <- mf[[grp_var]]

  canonical_order <- if (is.factor(grp_col)) {
    levels(grp_col) # respect existing factor levels
  } else {
    sort(unique(as.character(grp_col))) # alphabetical for character
  }

  # Map each level name → its canonical index
  level_index <- setNames(seq_along(canonical_order), canonical_order)

  plot_df <- data.frame(
    level = rownames(re_df),
    blup = blup_vals,
    se = se_vals,
    lower = blup_vals - 1.96 * se_vals,
    upper = blup_vals + 1.96 * se_vals,
    stringsAsFactors = FALSE
  )

  # Add the numeric index for each level
  plot_df$level_num <- level_index[plot_df$level]

  # Sort rows by BLUP value for the caterpillar shape
  plot_df <- plot_df[order(plot_df$blup), ]

  # y-axis: ordered factor of level numbers (in BLUP-sorted order)
  plot_df$level_num_fct <- factor(plot_df$level_num, levels = plot_df$level_num)

  # x-axis label: depends on model type
  xlab <- if (fam_info$is_glm) {
    paste0("BLUP (", fam_info$link, " scale)")
  } else {
    expression("BLUP (log"[10] * " scale)")
  }

  p <- ggplot(plot_df, aes(x = blup, y = level_num_fct)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    {
      if (all(!is.na(plot_df$se))) {
        geom_errorbar(
          aes(xmin = lower, xmax = upper),
          width = 0.3,
          colour = "steelblue",
          linewidth = 0.5,
          orientation = "y"
        )
      }
    } +
    geom_point(colour = "steelblue", size = 1.8) +
    labs(
      x = xlab,
      y = paste0(grp_var, " (level no.)"),
      title = paste0("Random: ", grp_var)
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(size = 9, face = "bold"),
      axis.text.y = element_text(size = 7),
      axis.ticks = element_line()
    )

  p
}


# ── Main function ──────────────────────────────────────────────────────────────

#' Plot marginal covariate effects for lm / lmer / glmer models
#'
#' Works correctly for:
#'   - `lmer` with a log10-transformed response (Gaussian, identity link on
#'     the log10 scale).  Set `response_scale = TRUE` to back-transform
#'     predictions to mg/m³ via 10^.
#'   - `glmer` with any family/link (e.g. Gamma + log).  `marginaleffects`
#'     returns predictions on the response scale when `type = "response"`,
#'     so no manual back-transform is applied.  The `response_scale` argument
#'     is ignored for glmer models.
#'
#' @param model        A fitted lm, lmer, or glmer object.
#' @param covariates   Named list mapping variable names to x-axis labels,
#'                     e.g. list(SST = "SST (°C)", Chl = "Chl-a (mg/m³)").
#'                     If NULL (default), all fixed-effect predictors are
#'                     auto-detected via insight::find_predictors() and their
#'                     raw variable names are used as axis labels.
#' @param interactions Optional list of 2-element character vectors naming
#'                     interaction pairs to plot as overlaid lines.
#'                     e.g. list(c("HarmDOY", "SST"), c("HarmTOD", "Depth"))
#' @param random_effects Character vector of random-effect grouping variable
#'                     names to plot as caterpillar charts of BLUPs, e.g.
#'                     c("Gear", "DatasetID").  If NULL (default), grouping
#'                     variables are auto-detected from the model via
#'                     insight::find_random().  Set to character(0) to suppress
#'                     all random-effect panels.
#' @param re.form      Passed to marginaleffects. NA (default) marginalises over
#'                     random effects (population-level prediction). NULL uses BLUPs.
#' @param response_scale Logical. For **lmer** models only: TRUE back-transforms
#'                     log10 predictions to mg/m³ via 10^. For **glmer** models
#'                     this argument is ignored — predictions are always on the
#'                     response scale. Default FALSE.
#' @param show_rug               Logical. TRUE adds a rug of observed data
#'                               locations along the x-axis of each continuous
#'                               covariate panel. Default TRUE.
#' @param show_partial_residuals Logical. TRUE overlays partial residuals on
#'                               each continuous covariate panel. Suppressed
#'                               automatically for glmer models. Default FALSE.
#' @param ncol         Columns in the patchwork layout. NULL = auto (max 3).
#'
#' @return A patchwork object. Save with ggplot2::ggsave().

fPlotBiomassLMERs <- function(
  model,
  covariates = NULL,
  interactions = NULL,
  random_effects = NULL,
  re.form = NA,
  response_scale = TRUE,
  show_rug = TRUE,
  show_partial_residuals = FALSE,
  ncol = NULL
) {
  # Detect model family / link once and pass to all helpers
  fam_info <- .model_family_info(model)

  # Auto-build covariates from model predictors if not supplied
  if (is.null(covariates)) {
    pred_names <- insight::find_predictors(model)$conditional
    covariates <- setNames(as.list(pred_names), pred_names)
  }

  # Auto-detect interaction pairs from the model formula if not supplied
  if (is.null(interactions)) {
    interactions <- .extract_interactions(model)
  }

  # Auto-detect random-effect grouping variables if not supplied.
  # NULL means "auto-detect"; character(0) means "plot none".
  if (is.null(random_effects)) {
    random_effects <- .extract_random_grouping_vars(model)
  }

  # 1. Individual covariate panels
  panels <- mapply(
    FUN = .plot_one_covariate,
    var = names(covariates),
    xlab = unname(covariates),
    MoreArgs = list(
      model = model,
      re.form = re.form,
      response_scale = response_scale,
      show_rug = show_rug,
      show_partial_residuals = show_partial_residuals,
      fam_info = fam_info
    ),
    SIMPLIFY = FALSE
  )

  # 2. Interaction panels
  if (!is.null(interactions) && length(interactions) > 0L) {
    int_panels <- lapply(interactions, function(vars) {
      .plot_interaction(
        model,
        vars,
        re.form = re.form,
        response_scale = response_scale,
        fam_info = fam_info
      )
    })
    panels <- c(panels, int_panels)
  }

  # 3. Random-effect caterpillar panels
  if (length(random_effects) > 0L) {
    re_panels <- lapply(random_effects, function(grp) {
      .plot_random_effect(model, grp, fam_info, response_scale)
    })
    # Drop any NULLs (e.g. grouping var not found in ranef)
    re_panels <- Filter(Negate(is.null), re_panels)
    panels <- c(panels, re_panels)
  }

  # 4. Assemble
  if (is.null(ncol)) {
    ncol <- min(3L, length(panels))
  }

  wrap_plots(panels, ncol = ncol) +
    plot_annotation(
      theme = theme(plot.background = element_rect(fill = "white", colour = NA))
    )
}
