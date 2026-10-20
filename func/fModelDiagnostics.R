# fModelDiagnostics.R
# Render a self-contained HTML model diagnostic report via Quarto.
#
# Sections:
#   1. Model Parameters  — full model_parameters() output (fixed effects,
#                          random effects variance, fit indices)
#   2. Model Diagnostics — check_model() panel via check_model_pw()
#   3. Covariate Effects — fPlotBiomassLMERs() marginal effect plots
#   4. Model Comparison  — compare_performance() table + radar plot
#
# The model is serialised to a temp RDS file and passed to the Quarto template
# as a parameter, so each section renders in a clean chunk with no device
# auto-print side effects.
#
# Author: Jason Everett
# Updated: October 2026
#
# Dependencies:
#   quarto (R package + CLI), withr, parameters, performance, see, patchwork
#   marginaleffects, ggplot2, rlang, insight   (via fPlotBiomassLMERs)
#
# Usage:
#   source("func/fModelDiagnostics.R")
#   fModelDiagnosticsHTML(
#     model          = m1,
#     file           = "Figures/m1_diagnostics.html",
#     compare_models = list(m4, m5, m6, m7, m8)   # optional
#   )

library(quarto)


#' Render a self-contained HTML model diagnostic report
#'
#' @param model                  A fitted lm, lmer, or glmer object.
#' @param file                   Output HTML path, e.g. "Figures/m1_diagnostics.html".
#' @param model_label            Optional subtitle label. Defaults to the formula.
#' @param compare_models         Optional list of fitted model objects to compare
#'                               against `model` via compare_performance(). When
#'                               non-NULL a Section 4 "Model Comparison" is added
#'                               to the report with a ranked table and radar plot.
#'                               Default NULL (section omitted).
#' @param covariates             Named list of variable → x-axis label pairs passed
#'                               to fPlotBiomassLMERs(). NULL = auto-detect.
#' @param interactions           Optional list of 2-element character vectors for
#'                               interaction panels, passed to fPlotBiomassLMERs().
#' @param random_effects         Character vector of random-effect grouping variable
#'                               names to plot as caterpillar charts of BLUPs, e.g.
#'                               c("Gear", "DatasetID"). NULL (default) = auto-detect
#'                               from the model. Use character(0) to suppress all
#'                               random-effect panels.
#' @param re.form                NA (default) marginalises over random effects;
#'                               NULL uses BLUPs.
#' @param response_scale         Logical. Back-transform to response scale? Default FALSE.
#' @param show_rug               Logical. Show rug on continuous covariate panels.
#' @param show_partial_residuals Logical. Overlay partial residuals on covariate panels.
#' @param ncol                   Columns in the covariate-effects patchwork. NULL = auto.
#' @param n_simulations          Number of simulations for the posterior predictive
#'                               check inside check_model(). Lower values are faster;
#'                               the performance package default is 250. Default here
#'                               is 100 to keep rendering time reasonable for large
#'                               mixed models.
#' @param check_model_checks     Character vector of diagnostic panels to include,
#'                               passed as the `check` argument to check_model().
#'                               NULL (default) includes all panels. Valid values:
#'                               "pp_check", "linearity", "homogeneity", "outliers",
#'                               "normality", "reqq" (random-effects QQ), "vif".
#'                               Example: c("pp_check", "normality", "reqq") for a
#'                               fast subset on very large models.
#' @param ci_method              Passed to model_parameters(). NULL = package default.
#'
#' @return Invisibly returns the output file path.
#'
#' @examples
#' fModelDiagnosticsHTML(
#'   model          = m1,
#'   file           = "Figures/m1_diagnostics.html",
#'   re.form        = NA,
#'   response_scale = FALSE
#' )

fModelDiagnosticsHTML <- function(
  model,
  file = NULL,
  model_label = NULL,
  compare_models = NULL,
  covariates = NULL,
  interactions = NULL,
  random_effects = NULL,
  re.form = NA,
  response_scale = FALSE,
  show_rug = TRUE,
  show_partial_residuals = FALSE,
  ncol = NULL,
  ci_method = NULL,
  n_simulations = 100L,
  check_model_checks = NULL
) {
  # ── Capture model name before any evaluation ─────────────────────────────────
  # deparse(substitute()) must be called first, before `model` is touched.
  model_name <- deparse(substitute(model))

  # ── Capture compare_models variable names ────────────────────────────────────
  # substitute() on the compare_models argument gives us the unevaluated call,
  # e.g. list(m4, m5, m6, m7, m8).  We deparse each element to get the symbol
  # names the user typed, then use them to name the list so compare_performance()
  # labels the output correctly.
  compare_models_call <- substitute(compare_models)
  if (!is.null(compare_models_call) && is.call(compare_models_call)) {
    # Elements 2..n of the call are the individual model arguments
    cmp_syms <- as.list(compare_models_call)[-1L]
    cmp_names <- vapply(cmp_syms, deparse, character(1))
    if (length(cmp_names) == length(compare_models)) {
      names(compare_models) <- cmp_names
    }
  }

  # ── Default label and file from model name ───────────────────────────────────
  if (is.null(model_label)) {
    model_label <- model_name
  }
  if (is.null(file) || !nzchar(file)) {
    file <- file.path("Figures", paste0(model_name, ".html"))
    message("No `file` supplied — writing to: ", file)
  }

  out_dir <- dirname(file)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  # ── Serialise model + data (+ comparison models) to temp RDS ────────────────
  # marginaleffects::predictions() looks up the original data frame by name in
  # the calling environment. In the Quarto session that environment is empty, so
  # we bundle the model's data alongside it and restore it in the template.
  tmp_rds <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp_rds), add = TRUE)
  model_data <- tryCatch(insight::get_data(model), error = function(e) NULL)
  saveRDS(
    list(model = model, data = model_data, compare_models = compare_models),
    tmp_rds
  )

  # ── Resolve absolute paths (Quarto runs in its own working directory) ────────
  template_path <- normalizePath(
    "func/model_diagnostics_template.qmd",
    mustWork = TRUE
  )
  utils_path <- normalizePath("ZooBiomassAus_0_Utils.R", mustWork = TRUE)
  func_path <- normalizePath("func/fPlotBiomassLMERs.R", mustWork = TRUE)
  out_abs <- normalizePath(file, mustWork = FALSE)

  # ── Build Quarto execute-params list ─────────────────────────────────────────
  # Note: Quarto cannot serialise character(0) as a YAML param — it becomes NULL.
  # We pass random_effects as-is; the template treats NULL as "auto-detect" and
  # fPlotBiomassLMERs() calls .extract_random_grouping_vars() internally.
  # compare_models is passed via the RDS bundle (not as a YAML param) because
  # Quarto cannot serialise arbitrary R objects in execute_params.
  execute_params <- list(
    rds_path = tmp_rds,
    model_label = model_label,
    utils_path = utils_path,
    func_path = func_path,
    covariates = covariates,
    interactions = interactions,
    random_effects = if (length(random_effects) == 0L) NULL else random_effects,
    re_form_na = is.na(re.form),
    response_scale = response_scale,
    show_rug = show_rug,
    show_partial_residuals = show_partial_residuals,
    ncol = ncol,
    ci_method = ci_method,
    has_comparison = !is.null(compare_models),
    n_simulations = n_simulations,
    check_model_checks = if (length(check_model_checks) == 0L) {
      NULL
    } else {
      check_model_checks
    }
  )

  # ── Render ───────────────────────────────────────────────────────────────────
  # Quarto places the output next to the .qmd by default.
  # We render there, then move to the requested destination.
  quarto::quarto_render(
    input = template_path,
    output_format = "html",
    output_file = basename(out_abs),
    execute_params = execute_params,
    quiet = FALSE
  )

  # The rendered file lands next to the template
  rendered <- file.path(dirname(template_path), basename(out_abs))

  if (!file.exists(rendered)) {
    stop(
      "Quarto rendered successfully but the output file was not found at: ",
      rendered
    )
  }

  file.rename(rendered, out_abs)

  message("Model diagnostic report written to: ", out_abs)
  invisible(file)
}
