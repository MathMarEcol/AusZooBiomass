# Zooplankton Biomass — Statistical Models
# Jason Everett
# Updated: October 2026

# ── Libraries ──────────────────────────────────────────────────────────────────
library(tidyverse)
library(splines)
library(lme4)
library(lmerTest) # p-values for lmer fixed effects
library(utils4mme) # fHarmonic()
library(marginaleffects) # covariate effect plots
library(performance) # r2(), check_model(), compare_performance()
library(see) # required backend for check_model() plots
library(parameters) # tidy coefficient tables
library(DHARMa) # simulation-based residuals for GLMMs
library(patchwork) # multi-panel figure assembly

source("ZooBiomassAus_0_Utils.R")
source("func/fPlotBiomassLMERs.R")
source("func/fModelDiagnostics.R")


# Load Data -----
# Prepared `dat` in ZooBiomassAus_1a_PrepareDatabase.R
dat <- read_csv(file = "Output/FinalBiomassData.csv") %>%
  mutate(log10_Biomass = log10(Biomass))


# Models -----

## M1 LMER (Gaussian, log10 response)
# Conditional R2: 0.818; Marginal R2: 0.574
m1 <- lmer(
  log10_Biomass ~ Method +
    log10(Chl) +
    ns(Bathy, df = 2) +
    harmonic(HarmDOY, k = 1) +
    Depth +
    Mesh +
    (1 | Gear),
  data = dat
)

## Quick visreg covariate plots ----
# GREAT. Looks identical
library(visreg)
par(mfrow = c(2, 3))
visreg(m1, "Chl", scale = "linear", gg = FALSE)
visreg(m1, "Bathy", scale = "linear", gg = FALSE)
visreg(m1, "HarmDOY", scale = "linear", gg = FALSE)
visreg(m1, "Depth", scale = "linear", gg = FALSE)
visreg(m1, "Mesh", scale = "linear", gg = FALSE)
visreg(m1, "Method", scale = "linear", gg = FALSE)
par(mfrow = c(1, 1))


## Summary — HTML diagnostic report ----
fModelDiagnosticsHTML(
  model = m1,
  re.form = NA,
  response_scale = FALSE
)


## m2 ------
# Conditional R2: 0.820; Marginal R2: 0.576
m2 <- lmer(
  log10_Biomass ~ Method +
    log10(Chl) +
    ns(Bathy, df = 2) +
    harmonic(HarmDOY, k = 1) +
    exp(-Depth / 1000) +
    Mesh +
    harmonic(HarmTOD, k = 1) +
    (1 | Gear),
  data = dat
)


## Summary — HTML diagnostic report ----
fModelDiagnosticsHTML(
  model = m2,
  re.form = NA,
  response_scale = FALSE,
)


## m3 ------
# Add intereaction between TOD and Depth
# Conditional R2: 0.820; Marginal R2: 0.581
# m2 still the best overall but this has the best marginal
m3 <- lmer(
  log10_Biomass ~ Method +
    log10(Chl) +
    ns(Bathy, df = 2) +
    harmonic(HarmDOY, k = 1) +
    Mesh +
    harmonic(HarmTOD, k = 1) * exp(-Depth / 1000) +
    (1 | Gear),
  data = dat
)


## Summary — HTML diagnostic report ----
fModelDiagnosticsHTML(
  model = m3,
  re.form = NA,
  response_scale = FALSE,
)

## m4 ------
# Try a log link on m1
# Conditional R2: 0.789; Marginal R2: 0.646
# Best marginal.
m4 <- glmer(
  Biomass ~ Method +
    log10(Chl) +
    ns(Bathy, df = 2) +
    harmonic(HarmDOY, k = 1) +
    Depth +
    Mesh +
    (1 | Gear),
  family = Gamma(link = "log"),
  data = dat
)


## Summary — HTML diagnostic report ----
fModelDiagnosticsHTML(
  model = m4,
  re.form = NA,
  response_scale = FALSE,
)

## m5 ------
# m3 with gamma log link
# Conditional R2: 0.792; Marginal R2: 0.648
# Best marginal
m5 <- glmer(
  Biomass ~ Method +
    log10(Chl) +
    ns(Bathy, df = 2) +
    harmonic(HarmDOY, k = 1) +
    Mesh +
    harmonic(HarmTOD, k = 1) * exp(-Depth / 1000) +
    (1 | Gear),
  family = Gamma(link = "log"),
  data = dat
)


## Summary — HTML diagnostic report ----
fModelDiagnosticsHTML(
  model = m5,
  re.form = NA,
  response_scale = FALSE,
)


## m6 ------
# m5 + Latitude * DOY
# CConditional R2: 0.809; Marginal R2: 0.632
#

m6 <- glmer(
  Biomass ~ Method +
    log10(Chl) +
    ns(Bathy, df = 2) +
    Mesh +
    harmonic(HarmTOD, k = 1) * exp(-Depth / 1000) +
    ns(Latitude, 2) * harmonic(HarmDOY, k = 2) +
    (1 | Gear),
  family = Gamma(link = "log"),
  data = dat,
  control = glmerControl(
    calc.derivs = FALSE,
    optimizer = "nloptwrap",
    optCtrl = list(algorithm = "NLOPT_LN_BOBYQA")
  )
)


## Summary — HTML diagnostic report ----
fModelDiagnosticsHTML(
  model = m6,
  re.form = NA,
  response_scale = FALSE,
)

## m7 ------
# m6 replacing Gear and Mesh with Gear_Mesh
# Conditional R2: 0.854; Marginal R2: 0.758
# Best Conditional and Marginal so far

m7 <- glmer(
  Biomass ~ Method +
    log10(Chl) +
    ns(Bathy, df = 2) +
    harmonic(HarmTOD, k = 1) * exp(-Depth / 1000) +
    ns(Latitude, 2) * harmonic(HarmDOY, k = 2) +
    (1 | Gear_Mesh),
  family = Gamma(link = "log"),
  data = dat,
  control = glmerControl(
    calc.derivs = FALSE,
    optimizer = "nloptwrap",
    optCtrl = list(algorithm = "NLOPT_LN_BOBYQA")
  )
)


## Summary — HTML diagnostic report ----
fModelDiagnosticsHTML(
  model = m7,
  re.form = NA,
  response_scale = FALSE,
)


## m8 ------
# m7 + (1 | Project)
# Conditional R2: 0.881; Marginal R2: 0.692
# Best Conditional but worse Marginal

m8 <- glmer(
  Biomass ~ Method +
    log(Chl) +
    ns(Bathy, df = 2) +
    harmonic(HarmTOD, k = 1) * exp(-Depth / 1000) +
    ns(Latitude, 2) * harmonic(HarmDOY, k = 2) +
    (1 | Gear_Mesh) +
    (1 | Project),
  family = Gamma(link = "log"),
  data = dat,
  control = glmerControl(
    calc.derivs = FALSE,
    optimizer = "nloptwrap",
    optCtrl = list(algorithm = "NLOPT_LN_BOBYQA")
  )
)


## Summary — HTML diagnostic report ----
fModelDiagnosticsHTML(
  model = m8,
  re.form = NA,
  response_scale = FALSE,
)


compare_performance(m4, m5, m6, m7, m8, rank = TRUE)
plot(compare_performance(m4, m5, m6, m7, m8))

# ══════════════════════════════════════════════════════════════════════════════
#### EXAMPLE 1 — GLMM (Gamma, log link) ####
# Use when response is strictly positive and right-skewed (raw Biomass).
# DHARMa residuals are ESSENTIAL — standard Pearson/deviance residuals are
# not valid for GLMMs.
# ══════════════════════════════════════════════════════════════════════════════

# ══════════════════════════════════════════════════════════════════════════════
#### EXAMPLE 2 — LMER (Gaussian, log10 response) — Nets ####
# Use when response has been log10-transformed to meet normality assumptions.
# Random effect accounts for variance across gear type × mesh size combinations.
# ══════════════════════════════════════════════════════════════════════════════
