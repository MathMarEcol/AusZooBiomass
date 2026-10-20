# Zooplankton Biomass — Statistical Models
# Jason Everett
# Updated: October 2026

# ── Libraries ──────────────────────────────────────────────────────────────────
library(tidyverse)
library(splines)
library(utils4mme) # fHarmonic()
library(marginaleffects) # covariate effect plots
library(performance) # r2(), check_model(), compare_performance()
library(see) # required backend for check_model() plots
library(parameters) # tidy coefficient tables
library(DHARMa) # simulation-based residuals for GLMMs
library(patchwork) # multi-panel figure assembly
library(glmmTMB)

source("ZooBiomassAus_0_Utils.R")
source("func/fPlotBiomassLMERs.R")
source("func/fModelDiagnostics.R")


# Load Data -----
# Prepared `dat` in ZooBiomassAus_1a_PrepareDatabase.R
dat <- read_rds(file = "Output/FinalBiomassData.rds")


# Models -----

# ## m01 LMER (Gaussian, log10 response)
# # Conditional R2: 0.818; Marginal R2: 0.574
# m01 <- glmmTMB(
#   log10_Biomass ~ Method +
#     log10(Chl) +
#     ns(Bathy, df = 2) +
#     harmonic(HarmDOY, k = 1) +
#     Depth +
#     Mesh +
#     (1 | Gear),
#   data = dat
# )

## Summary — HTML diagnostic report ----
# fModelDiagnosticsHTML(
#   model = m01,
#   re.form = NA,
#   response_scale = FALSE
# )

## m02 ------
# Conditional R2: 0.820; Marginal R2: 0.576
# m02 <- glmmTMB(
#   log10_Biomass ~ Method +
#     log10(Chl) +
#     ns(Bathy, df = 2) +
#     harmonic(HarmDOY, k = 1) +
#     Depth_decay +
#     Mesh +
#     harmonic(HarmTOD, k = 1) +
#     (1 | Gear),
#   data = dat
# )

# ## Summary — HTML diagnostic report ----
# fModelDiagnosticsHTML(
#   model = m02,
#   re.form = NA,
#   response_scale = FALSE,
# )

# ## m03 ------
# # Add intereaction between TOD and Depth
# # Conditional R2: 0.820; Marginal R2: 0.581
# # m02 still the best overall but this has the best marginal
# m03 <- glmmTMB(
#   log10_Biomass ~ Method +
#     log10(Chl) +
#     ns(Bathy, df = 2) +
#     harmonic(HarmDOY, k = 1) +
#     Mesh +
#     harmonic(HarmTOD, k = 1) * Depth_decay +
#     (1 | Gear),
#   data = dat
# )

# ## Summary — HTML diagnostic report ----
# fModelDiagnosticsHTML(
#   model = m03,
#   re.form = NA,
#   response_scale = FALSE,
# )

# ## m04 ------
# # Try a log link on m01
# # Conditional R2: 0.789; Marginal R2: 0.646
# # Best marginal.

# m04 <- glmmTMB(
#   Biomass ~ Method +
#     log10(Chl) +
#     ns(Bathy, df = 2) +
#     harmonic(HarmDOY, k = 1) +
#     Depth +
#     Mesh +
#     (1 | Gear),
#   family = Gamma(link = "log"),
#   data = dat
# )

# ## Summary — HTML diagnostic report ----
# fModelDiagnosticsHTML(
#   model = m04,
#   re.form = NA,
#   response_scale = FALSE,
# )

# ## m05 ------
# # m03 with gamma log link
# # Conditional R2: 0.792; Marginal R2: 0.648
# # Best marginal
# m05 <- glmmTMB(
#   Biomass ~ Method +
#     log10(Chl) +
#     ns(Bathy, df = 2) +
#     harmonic(HarmDOY, k = 1) +
#     Mesh +
#     harmonic(HarmTOD, k = 1) * Depth_decay +
#     (1 | Gear),
#   family = Gamma(link = "log"),
#   data = dat
# )

# ## Summary — HTML diagnostic report ----
# fModelDiagnosticsHTML(
#   model = m05,
#   re.form = NA,
#   response_scale = FALSE,
# )

# ## m06 ------
# # m05 + Latitude * DOY
# # CConditional R2: 0.809; Marginal R2: 0.632
# #

# m06 <- glmmTMB(
#   Biomass ~ Method +
#     log10(Chl) +
#     ns(Bathy, df = 2) +
#     Mesh +
#     harmonic(HarmTOD, k = 1) * Depth_decay +
#     ns(Latitude, 2) * harmonic(HarmDOY, k = 2) +
#     (1 | Gear),
#   family = Gamma(link = "log"),
#   data = dat,
#   control = glmmTMBControl(
#     optCtrl = list(iter.max = 1000, eval.max = 1000)
#   )
# )

# ## Summary — HTML diagnostic report ----
# fModelDiagnosticsHTML(
#   model = m06,
#   re.form = NA,
#   response_scale = FALSE,
# )

# ## m07 ------
# # m06 replacing Gear and Mesh with Gear_Mesh
# # Conditional R2: 0.854; Marginal R2: 0.758
# # Best Conditional and Marginal so far
# m07 <- glmmTMB(
#   Biomass ~ Method +
#     log10(Chl) +
#     ns(Bathy, df = 2) +
#     harmonic(HarmTOD, k = 1) * Depth_decay +
#     ns(Latitude, 2) * harmonic(HarmDOY, k = 2) +
#     (1 | Gear_Mesh),
#   family = Gamma(link = "log"),
#   data = dat,
#   control = glmmTMBControl(
#     optCtrl = list(iter.max = 1000, eval.max = 1000)
#   )
# )

# ## Summary — HTML diagnostic report ----
# fModelDiagnosticsHTML(
#   model = m07,
#   re.form = NA,
#   response_scale = FALSE,
# )

# ## m08 ------
# # m07 + (1 | Project)
# # Conditional R2: 0.881; Marginal R2: 0.692
# # Best Conditional but worse Marginal

# m08 <- glmmTMB(
#   Biomass ~ Method +
#     log(Chl) +
#     ns(Bathy, df = 2) +
#     harmonic(HarmTOD, k = 1) * Depth_decay +
#     ns(Latitude, 2) * harmonic(HarmDOY, k = 2) +
#     (1 | Gear_Mesh) +
#     (1 | Project),
#   family = Gamma(link = "log"),
#   data = dat,
#   control = glmmTMBControl(
#     optCtrl = list(iter.max = 1000, eval.max = 1000)
#   )
# )

# ## Summary — HTML diagnostic report ----
# fModelDiagnosticsHTML(
#   model = m08,
#   re.form = NA,
#   response_scale = FALSE,
#   compare_models = list(m04, m05, m06, m07, m08)
# )

# drop1(m08) - NOTHING
# drop1(m08, k = log(nobs(m08))) - NOTHING

## m09 ------
# m08 but remove gear_mesh
#
#

# m09 <- glmmTMB(
#   Biomass ~ Method +
#     log(Chl) +
#     ns(Bathy, df = 2) +
#     harmonic(HarmTOD, k = 1) * Depth_decay +
#     ns(Latitude, 2) * harmonic(HarmDOY, k = 2) +
#     (1 | Project),
#   family = Gamma(link = "log"),
#   data = dat,
#   control = glmmTMBControl(
#     optCtrl = list(iter.max = 1000, eval.max = 1000)
#   )
# )

# ## Summary — HTML diagnostic report ----
# fModelDiagnosticsHTML(
#   model = m09,
#   re.form = NA,
#   response_scale = FALSE,
#   compare_models = list(m04, m05, m06, m07, m08, m09)
# )

# m10 ------
# m8 with Tow
m10 <- glmmTMB(
  Biomass ~ Method +
    log(Chl) +
    ns(Bathy, df = 2) +
    harmonic(HarmTOD, k = 1) * Depth_decay +
    harmonic(HarmDOY, k = 2) * ns(Latitude, 2) +
    Tow +
    (1 | Gear_Mesh) +
    (1 | Project),
  family = Gamma(link = "log"),
  data = dat,
  control = glmmTMBControl(
    optCtrl = list(iter.max = 1000, eval.max = 1000),
    parallel = list(n = parallel::detectCores() - 5L, autopar = TRUE)
  )
)

write_rds(m10, "Output/m10.rds")

## Summary — HTML diagnostic report ----
fModelDiagnosticsHTML(
  model = m10,
  re.form = NA,
  response_scale = FALSE,
  compare_models = NULL
)

# m11 -----
# m10 without Gear_Mesh
m11 <- glmmTMB(
  Biomass ~ Method +
    log(Chl) +
    ns(Bathy, df = 2) +
    harmonic(HarmTOD, k = 1) * Depth_decay +
    harmonic(HarmDOY, k = 2) * ns(Latitude, 2) +
    Tow +
    (1 | Project),
  family = Gamma(link = "log"),
  data = dat,
  control = glmmTMBControl(
    optCtrl = list(iter.max = 1000, eval.max = 1000),
    parallel = list(n = parallel::detectCores() - 5L, autopar = TRUE)
  )
)
write_rds(m11, "Output/m11.rds")

## Summary — HTML diagnostic report ----
fModelDiagnosticsHTML(
  model = m11,
  re.form = NA,
  response_scale = FALSE,
  compare_models = NULL
)


# m12 - m10 without Lat DOY interaction
m12 <- glmmTMB(
  Biomass ~ Method +
    log(Chl) +
    ns(Bathy, df = 2) +
    harmonic(HarmTOD, k = 1) * Depth_decay +
    ns(Latitude, 2) +
    harmonic(HarmDOY, k = 2) +
    Tow +
    (1 | Gear_Mesh) +
    (1 | Project),
  family = Gamma(link = "log"),
  data = dat,
  control = glmmTMBControl(
    optCtrl = list(iter.max = 1000, eval.max = 1000),
    parallel = list(n = parallel::detectCores() - 5L, autopar = TRUE)
  )
)


## Summary — HTML diagnostic report ----
fModelDiagnosticsHTML(
  model = m12,
  re.form = NA,
  response_scale = FALSE,
  compare_models = NULL
)


# m13 ------
# m10 swapping Latitude for SST
m13 <- glmmTMB(
  Biomass ~ Method +
    log(Chl) +
    ns(Bathy, df = 2) +
    harmonic(HarmTOD, k = 1) * Depth_decay +
    harmonic(HarmDOY, k = 2) * ns(SST, 2) +
    Tow +
    (1 | Gear_Mesh) +
    (1 | Project),
  family = Gamma(link = "log"),
  data = dat,
  control = glmmTMBControl(
    optCtrl = list(iter.max = 1000, eval.max = 1000),
    parallel = list(n = parallel::detectCores() - 5L, autopar = TRUE)
  )
)

write_rds(m13, "Output/m13.rds")

## Summary — HTML diagnostic report ----
fModelDiagnosticsHTML(
  model = m13,
  re.form = NA,
  response_scale = FALSE,
  compare_models = list(m10, m11, m12, m13)
)

# # Change to Tweedie -----
1
# sv <- list(
#   beta = fixef(m11)$cond, # fixed-effect coefficients
#   theta = getME(m11, "theta") # random-effect log-SDs
#   # glmmTMB will initialise the Tweedie p (psi) at its default
# )

# m13 <- glmmTMB(
#   Biomass ~ Method +
#     log(Chl) +
#     ns(Bathy, df = 2) +
#     harmonic(HarmTOD, k = 1) * Depth_decay +
#     ns(Latitude, 2) * harmonic(HarmDOY, k = 2) +
#     Tow +
#     (1 | Project),
#   family = tweedie(link = "log"),
#   data = dat,
#   start = sv,
#   control = glmmTMBControl(
#     optCtrl = list(iter.max = 5000, eval.max = 5000),
#     parallel = list(n = parallel::detectCores() - 5L, autopar = TRUE)
#   )
# )

# write_rds(m13, "Output/m13.rds")

# ## Summary — HTML diagnostic report ----
# fModelDiagnosticsHTML(
#   model = m13,
#   re.form = NA,
#   response_scale = FALSE,
#   compare_models = list(m11, m13),
#   check_model_checks = c("reqq")
# )

# m9 with Tow

# TODO
# rename HTML report function
# Write some code for predition - mean and month
# Write some code to save the netcdf

# Why is 12/19 [covariate-effects] so slow?
