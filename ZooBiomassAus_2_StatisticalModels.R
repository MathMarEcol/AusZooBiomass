# Biomass analysis around Australia
# Ant & Jase 4/7/2019
# Data provided by Felicity 20/6/2019

##### PRELIMINARIES #####
library(tidyverse)
library(splines)
library(lme4)
library(MuMIn) # r2 for mixed models
library(lmerTest) # p-levels for mixed models
library(visreg)
library(utils4mme)
library(visreg)

source("func/fHarmonic.R") # Harmonic function
source("func/fPlotBiomassLM.R")
source("func/fPlotBiomassLMERs.R")

# dat <- readRDS("Data/BiomassData.rds") # n = 55,131
dat <- read_rds("../GlobalZoopBiomass/Output/dat_cleaned.rds") %>%
  sf::st_as_sf(
    coords = c("Longitude", "Latitude"),
    remove = FALSE,
    crs = "EPSG:4326"
  )

bndry <- sf::st_bbox(
  c(xmin = 100, xmax = 170, ymax = -8, ymin = -50),
  crs = sf::st_crs("EPSG:4326")
) %>%
  sf::st_as_sfc()

dat <- sf::st_filter(dat, bndry)

## Filter for Australia!

#### Spatial model for Biomass ####
# Use GLMM because can restrict the relationships a bit more where there are few data
# Note: Chl only has about 10 unique values but Chl_gm3 has more

system.time(
  m1_map <- glmer(
    Biomass ~ Method +
      Mesh +
      Tow +
      log(Chl) +
      ns(SST, df = 3) +
      ns(Bathy, df = 2) +
      harmonic(HarmTOD, k = 1) * ns(Depth, df = 3) +
      ns(Latitude, 2) * harmonic(HarmDOY, k = 2) +
      (1 | DatasetID) +
      (1 | Gear),
    family = Gamma(link = "log"),
    data = dat,
    control = glmerControl(
      calc.derivs = FALSE,
      optimizer = "nloptwrap",
      optCtrl = list(algorithm = "NLOPT_LN_BOBYQA")
    )
  )
  # harmonic(HarmDOY, k = 2)*ns(Thetao, df = 3)
)
#

summary(m1_map)

saveRDS(m1_map, "Output/m1_map.rds")

r.squaredGLMM(m1_map) # Use trigamma. 70.3% fixed effects; 84.7% total (with random effects)

# Check residuals
lattice::qqmath(m1_map) # Terrible without removing very high values, but still pretty poor
plot(m1_map, sqrt(abs(resid(.))) ~ fitted(.)) # Good fitted residuals, but some big ones at small values

# Check distribution of random effects

model <- readRDS("Output/m1_map.rds")

fPlotBiomassLMERs(m1_map, "Biomass", "log10")

# p_out <- PlotGLMM_Map(m1_map, FileName = "Plot_GLMM_m1_map")

dat <- dat %>%
  mutate(
    Time = (hour(DateLocal) * 60 + minute(DateLocal) + second(DateLocal) / 60) /
      (24 * 60),
    HarmTOD = Time * 2 * pi,
    HarmDOY = (DOY / 365) * 2 * pi, # Convert to radians
    Biomass = replace(Biomass, Biomass > 10000, 10000),
    Bathy = replace(Bathy, Bathy > 6000, 6000),
    Depth = replace(Depth, Depth > 500, 500),
    Mesh = replace(Mesh, Mesh > 400, 400),
    Mesh = as.numeric((Mesh)),
    Tow = as.character(Tow),
    SST = replace(SST, SST < 5, 5),
    SST = replace(SST, SST > 30, 30),
    # Tow = replace(Tow, Tow == "Unknown", "Horizontal"),
    BiomassMethod = as.factor(BiomassMethod),
    Project = as.character(Project)
  ) %>%
  droplevels()

# Calculate Gear_Mesh
dat <- dat %>% unite(Gear_Mes2, NetType, Mesh, remove = FALSE)

# Code OPC/LOPC
dat <- dat %>%
  mutate(
    LOPC_Type = case_when(
      grepl("SS", dat$Project, fixed = TRUE) ~ "OPC",
      grepl("IN", dat$Project, fixed = TRUE) ~ "LOPC"
    )
  )
dat$Gear_Mes2[dat$NetType2 == "LOPC"] <- dat$LOPC_Type[dat$NetType2 == "LOPC"]

# Create dfs for CPR, Net, LOPC
# Define each dataset
dat_Net <- dat %>% filter(NetType2 == "Net") %>% droplevels()
dat_CPR <- dat %>% filter(NetType2 == "CPR") %>% droplevels()
dat_LOPC <- dat %>% filter(NetType2 == "LOPC") %>% droplevels()

# Times are OK for CPR and LOPC, but lots of 00:00 time for nets, which need to have a random time
dat_Net$Time[dat_Net$Time == 0] <- runif(sum(dat_Net$Time == 0), 0, 1) # Random time between 0 and 24 hrs for those with Time = 0; n = 3419 00:00 time
dat_Net <- dat_Net %>% mutate(HarmTOD = Time * 2 * pi)

##### MODELS - NETS = LM (NO RANDOM EFFECT + NO GAMMA + Gear_Mes2) #####
m1 <- lmer(
  log10(Biomass) ~ BiomassMethod +
    log10(Chl) +
    fHarmonic(HarmTOD, k = 1) * exp(-Depth / 1000) +
    ns(Bathy, df = 3) +
    fHarmonic(HarmDOY, k = 1) * ns(SST, 3) +
    (1 | Gear_Mes2),
  data = dat_Net
)
summary(m1)
#
r.squaredGLMM(m1) # Fixed r2 = 71.4%, Total r2 = 87.4%
fPlotBiomassLMERs(m1, "Biomass", "log10")
dev.copy2pdf(file = "Figures/Biomass_Nets_Final.pdf", paper = "A4r")
saveRDS(m1, "ModelOutput/Biomass_Nets_Final.rds")

############################################################
##### MODELS - CPR = LM (NO RANDOM EFFECT + NO GAMMA + ALL SAME SURVEY SO NO RANDOM EFFECT AND NO DEPTH) #####
m1 <- lm(
  log10(Biomass) ~ log10(Chl) +
    fHarmonic(HarmTOD, k = 1) +
    ns(Bathy, df = 3) +
    fHarmonic(HarmDOY, k = 1) * ns(SST, 3),
  data = dat_CPR
)
summary(m1) # r2 = 19.5%
anova(m1)
fPlotBiomassLM(m1, "Biomass_lm1")
dev.copy2pdf(file = "Figures/Biomass_CPR_Final.pdf", paper = "A4r") #
# par(mfrow = c(2,2))
# plot(m1)
saveRDS(m1, "ModelOutput/Biomass_CPR_Final.rds")

##### MODELS - LOPC = LM (NO RANDOM EFFECT + NO GAMMA + FIXED EFFECT FOR OPC/LOPC) #####
m1 <- lm(
  log10(Biomass) ~ LOPC_Type +
    log10(Chl) +
    fHarmonic(HarmTOD, k = 1) * exp(-Depth / 1000) +
    ns(Bathy, df = 3) +
    fHarmonic(HarmDOY, k = 1) * ns(SST, 3),
  data = dat_LOPC
)
summary(m1) # r2 = 55.5%
anova(m1)
fPlotBiomassLM(m1, "Biomass_lm1")
dev.copy2pdf(file = "Figures/Biomass_LOPC_Final.pdf", paper = "A4r") #
# par(mfrow = c(2,2))
# plot(m1)
saveRDS(m1, "ModelOutput/Biomass_LOPC_Final.rds")

##### ***MODELS - ALL = LMER (RANDOM EFFECT + NO GAMMA) #####
m1 <- lmer(
  log10(Biomass) ~ BiomassMethod +
    log10(Chl) +
    fHarmonic(HarmTOD, k = 1) * exp(-Depth / 1000) +
    ns(Bathy, df = 3) +
    fHarmonic(HarmDOY, k = 1) * ns(SST, 3) +
    (1 | Gear_Mes2),
  data = dat
)
summary(m1)
r.squaredGLMM(m1) # Fixed r2 = 58.3%, Total r2 = 82.3%
fPlotBiomassLMERs(m1, "Biomass", "log10")
dev.copy2pdf(
  file = "Figures/Biomass_All_Final_UnequalWeighting.pdf",
  paper = "A4r"
) #
saveRDS(m1, "ModelOutput/Biomass_All_Final_UnequalWeighting.rds")

# Give Nets, CPR and LOPC equal weighting
Summary <- dat %>% group_by(NetType2) %>% summarise(N = n())
dat <- dat %>%
  mutate(
    WtVec = case_when(
      NetType2 == "Net" ~ 1 / Summary$N[3] * nrow(dat),
      NetType2 == "CPR" ~ 1 / Summary$N[1] * nrow(dat),
      NetType2 == "LOPC" ~ 1 / Summary$N[2] * nrow(dat)
    )
  )
Mn <- mean(dat$WtVec)
dat$WtVec <- dat$WtVec / Mn
sum(dat$WtVec) # Check - correct!

m1 <- lmer(
  log10(Biomass) ~ BiomassMethod +
    log10(Chl) +
    fHarmonic(HarmTOD, k = 1) * exp(-Depth / 1000) +
    ns(Bathy, df = 3) +
    fHarmonic(HarmDOY, k = 1) * ns(SST, 3) +
    (1 | Gear_Mes2),
  weights = WtVec,
  data = dat
)
summary(m1)
r.squaredGLMM(m1) # Fixed r2 = 59.5%, Total r2 = 82.4%
fPlotBiomassLMERs(m1, "Biomass", "log10")
dev.copy2pdf(
  file = "Figures/BiomassModel_All_Final_EqualWeighting.pdf",
  paper = "A4r"
) #
saveRDS(m1, "ModelOutput/Biomass_All_Final_EqualWeighting.rds")
