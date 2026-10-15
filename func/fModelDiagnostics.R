# fModelDiagnostics.R
# Render a self-contained HTML model diagnostic report via Quarto.
#
# Sections:
#   1. Model Parameters  — full model_parameters() output (fixed effects,
#                          random effects variance, fit indices)
#   2. Model Diagnostics — check_model() panel via check_model_pw()
#   3. Covariate Effects — fPlotBiomassLMERs() marginal effect plots
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
#     model = m1,
#     file  = "Figures/m1_diagnostics.html"
#   )

library(quarto)


#' Render a self-contained HTML model diagnostic report
#'
#' @param model                  A fitted lm, lmer, or glmer object.
#' @param file                   Output HTML path, e.g. "Figures/m1_diagnostics.html".
#' @param model_label            Optional subtitle label. Defaults to the formula.
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
  covariates = NULL,
  interactions = NULL,
  random_effects = NULL,
  re.form = NA,
  response_scale = FALSE,
  show_rug = TRUE,
  show_partial_residuals = FALSE,
  ncol = NULL,
  ci_method = NULL
) {
  # ── Capture model name before any evaluation ─────────────────────────────────
  # deparse(substitute()) must be called first, before `model` is touched.
  model_name <- deparse(substitute(model))

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

  # ── Serialise model + data to temp RDS ──────────────────────────────────────
  # marginaleffects::predictions() looks up the original data frame by name in
  # the calling environment. In the Quarto session that environment is empty, so
  # we bundle the model's data alongside it and restore it in the template.
  tmp_rds <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp_rds), add = TRUE)
  model_data <- tryCatch(insight::get_data(model), error = function(e) NULL)
  saveRDS(list(model = model, data = model_data), tmp_rds)

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
    ci_method = ci_method
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
