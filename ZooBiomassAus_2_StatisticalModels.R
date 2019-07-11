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
library(lubridate)

source("../ZoopAbundanceModels/fHarmonic.R") # Harmonic function
source("fPlotBiomassLM.R")
source("fPlotBiomassLMERs.R")

dat <- readRDS("BiomassData.rds") # n = 55,131
dat <- dat %>% 
  mutate(
    Time = (hour(DateLocal) * 60 + minute(DateLocal) + second(DateLocal) / 60) / (24 * 60),
    HarmTOD = Time * 2 * pi,
    HarmDOY = (DOY/365)*2*pi, # Convert to radians
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

# Create dfs for CPR, Net, LOPC
# Define each dataset
dat_Net <- dat %>% filter(NetType2 == "Net") %>% droplevels()
dat_CPR <- dat %>% filter(NetType2 == "CPR") %>% droplevels()
dat_LOPC <- dat %>% filter(NetType2 == "LOPC") %>% droplevels()
dat_LOPC <- dat_LOPC %>% mutate(LOPC_Type = case_when(grepl("SS", dat_LOPC$Project, fixed = TRUE) ~ "OPC",
                                                      grepl("IN", dat_LOPC$Project, fixed = TRUE) ~ "LOPC"))
                                
# Times are OK for CPR and LOPC, but lots of 00:00 time for nets, which need to have a random time
dat_Net$Time[dat_Net$Time == 0] <- runif(sum(dat_Net$Time == 0), 0, 1) # Random time between 0 and 24 hrs for those with Time = 0; n = 3419 00:00 time
dat_Net <- dat_Net %>% mutate(HarmTOD = Time * 2 * pi)

##### MODELS - NETS = LM (NO RANDOM EFFECT + NO GAMMA) #####
m1 <- lm(log10(Biomass) ~ BiomassMethod + Tow + Mesh + log10(Chl) + 
           fHarmonic(HarmTOD, k = 1) * exp(-Depth/1000) + ns(Bathy, df = 3) + 
           fHarmonic(HarmDOY, k = 1) * ns(SST, 3), 
         data = dat_Net)
# With Gear_Mes2
m1 <- lm(log10(Biomass) ~ BiomassMethod + Gear_Mes2 + log10(Chl) + 
           fHarmonic(HarmTOD, k = 1) * exp(-Depth/1000) + ns(Bathy, df = 3) + 
           fHarmonic(HarmDOY, k = 1) * ns(SST, 3), 
         data = dat_Net)
summary(m1) # r2 = 75.0%
anova(m1)
fPlotBiomassLM(m1, "Biomass_lm1")
dev.copy2pdf(file = "BiomassModel_Nets_GearMes2_Fixed.pdf", paper = "A4r") #
# par(mfrow = c(2,2))
# plot(m1)

##### ***MODELS - NETS = LMER (RANDOM EFFECT + NO GAMMA) #####
# LMER (random effect but no Gamma())
m1 <- lmer(log10(Biomass) ~ BiomassMethod + Tow + Mesh + log10(Chl) + 
             fHarmonic(HarmTOD, k = 1) * exp(-Depth/1000) + ns(Bathy, df = 3) + 
              fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + (1|Project), data = dat_Net)
summary(m1)
r.squaredGLMM(m1) # Fixed r2 = 64.6%, Total r2 = 92.1%
fPlotBiomassLMERs(m1, "Biomass", "log10")
dev.copy2pdf(file = "BiomassModel_Nets.pdf", paper = "A4r") #
saveRDS(m1, "LMER_Nets.rds")

x11(width = 8, height = 5, title = "Biomass Diagnostics")
par(mfrow = c(1,2))
plot(m1)
qqnorm(residuals(m1))
qqline(residuals(m1), col = "blue")
plot(ranef(m1)$Project, main = "Project")

##### ***MODELS - NETS = LMER (RANDOM EFFECT + NO GAMMA + GEAR_MESH) #####
m1 <- lmer(log10(Biomass) ~ BiomassMethod + log10(Chl) + 
             fHarmonic(HarmTOD, k = 1) * exp(-Depth/1000) + ns(Bathy, df = 3) + 
             fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + (1|Gear_Mes2), data = dat_Net)
summary(m1)
r.squaredGLMM(m1) # Fixed r2 = 64.6%, Total r2 = 92.1%
fPlotBiomassLMERs(m1, "Biomass", "log10")
dev.copy2pdf(file = "BiomassModel_Nets_GearMes2_Random.pdf", paper = "A4r") #
saveRDS(m1, "LMER_Nets.rds")


##### MODELS - NETS = GLM (NO RANDOM EFFECT + GAMMA) #####
m1 <- glm(Biomass ~ BiomassMethod + Tow + Mesh + log10(Chl) + 
            fHarmonic(HarmTOD, k = 1) * exp(-Depth/1000) + ns(Bathy, df = 3) + 
            fHarmonic(HarmDOY, k = 1) * ns(SST, 3), data = dat_Net,  maxit = 1000, family = Gamma(link = log))
summary(m1)
anova(m1, test = "F")
1-(m1$deviance/m1$null.deviance) # r2 = 69.3%

# par(mfrow = c(2,2))
# plot(m1) # Residuals actually look better for the log10 transformation rather than log-link
fPlotBiomassLM(m1, "Biomass") # Surface is bad

##### MODELS - NETS = GLMER (RANDOM EFFECT + GAMMA) #####
m1 <- glmer(Biomass ~ BiomassMethod + Tow + Mesh + log10(Chl) + 
              fHarmonic(HarmTOD, k = 1) * exp(-Depth/1000) + ns(Bathy, df = 3) + 
              fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + (1|Project), data = dat_Net, family = Gamma(link = log), nAGQ=0)
m1 <- glmer(Biomass ~ BiomassMethod + log10(Chl) + 
              fHarmonic(HarmTOD, k = 1) * exp(-Depth/1000) + ns(Bathy, df = 3) + 
              fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + (1|Gear_Mes2), data = dat_Net, family = Gamma(link = log), nAGQ=0)
summary(m1)
r.squaredGLMM(m1)
fPlotBiomassLMERs(m1, "Biomass", "log")

x11(width = 8, height = 5, title = "Biomass Diagnostics")
par(mfrow = c(1,2))
plot(m1)
qqnorm(residuals(m1))
qqline(residuals(m1), col = "blue")
plot(ranef(m1)$Project, main = "Project")

############################################################
##### MODELS - CPR = LM (NO RANDOM EFFECT + NO GAMMA) #####
m1 <- lm(log10(Biomass) ~ log10(Chl) + 
           fHarmonic(HarmTOD, k = 1) + ns(Bathy, df = 3) + 
           fHarmonic(HarmDOY, k = 1) * ns(SST, 3), 
         data = dat_CPR)
summary(m1) # r2 = 19.5%
anova(m1)
fPlotBiomassLM(m1, "Biomass_lm1")
# par(mfrow = c(2,2))
# plot(m1)

##### ***MODELS - CPR = LMER (RANDOM EFFECT + NO GAMMA) #####
# LMER (random effect but no Gamma())
m1 <- lmer(log10(Biomass) ~ log10(Chl) + 
             fHarmonic(HarmTOD, k = 1) + ns(Bathy, df = 3) + 
             fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + (1|Project), data = dat_CPR)
summary(m1)
r.squaredGLMM(m1) # Fixed r2 = 13.3%, Total r2 = 27.2%
fPlotBiomassLMERs(m1, "Biomass", "log10") # Surface is good
dev.copy2pdf(file = "BiomassModel_CPR.pdf", paper = "A4r") #
saveRDS(m1, "LMER_CPR.rds")

x11(width = 8, height = 5, title = "Biomass Diagnostics")
par(mfrow = c(1,2))
plot(m1)
qqnorm(residuals(m1))
qqline(residuals(m1), col = "blue")
plot(ranef(m1)$Project, main = "Project")

##### ***MODELS - CPR = LMER (RANDOM EFFECT + NO GAMMA + GEAR_MESH) #####
m1 <- lmer(log10(Biomass) ~ BiomassMethod + log10(Chl) + 
             fHarmonic(HarmTOD, k = 1) + ns(Bathy, df = 3) + 
             fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + (1|Gear_Mes2), data = dat_CPR)
summary(m1)
r.squaredGLMM(m1) # Fixed r2 = 64.6%, Total r2 = 92.1%
fPlotBiomassLMERs(m1, "Biomass", "log10")
dev.copy2pdf(file = "BiomassModel_Nets.pdf", paper = "A4r") #
saveRDS(m1, "LMER_Nets.rds")



##### MODELS - CPR = GLM (NO RANDOM EFFECT + GAMMA) #####
m1 <- glm(Biomass ~ log10(Chl) + 
            fHarmonic(HarmTOD, k = 1) + ns(Bathy, df = 3) + 
            fHarmonic(HarmDOY, k = 1) * ns(SST, 3), data = dat_CPR,  maxit = 1000, family = Gamma(link = log))
summary(m1)
anova(m1, test = "F")
1-(m1$deviance/m1$null.deviance) # r2 = 21.4%

# par(mfrow = c(2,2))
# plot(m1) # Residuals actually look better for the log10 transformation rather than log-link
fPlotBiomassLM(m1, "Biomass") # Surface is bad

##### MODELS - CPR = GLMER (RANDOM EFFECT + GAMMA) #####
m1 <- glmer(Biomass ~ log10(Chl) + 
              fHarmonic(HarmTOD, k = 1) + ns(Bathy, df = 3) + 
              fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + (1|Project), data = dat_CPR, family = Gamma(link = log), nAGQ=0)
summary(m1)
r.squaredGLMM(m1)
fPlotBiomassLMERs(m1, "Biomass", "log") # Surface bad

x11(width = 8, height = 5, title = "Biomass Diagnostics")
par(mfrow = c(1,2))
plot(m1)
qqnorm(residuals(m1))
qqline(residuals(m1), col = "blue")
plot(ranef(m1)$Project, main = "Project")

##### MODELS - LOPC = LM (NO RANDOM EFFECT + NO GAMMA) #####
m1 <- lm(log10(Biomass) ~ log10(Chl) + 
           fHarmonic(HarmTOD, k = 1) * exp(-Depth/1000) + ns(Bathy, df = 3) + 
           fHarmonic(HarmDOY, k = 1) * ns(SST, 3), 
         data = dat_LOPC)
m1 <- lm(log10(Biomass) ~ LOPC_Type + log10(Chl) + 
           fHarmonic(HarmTOD, k = 1) * exp(-Depth/1000) + ns(Bathy, df = 3) + 
           fHarmonic(HarmDOY, k = 1) * ns(SST, 3), 
         data = dat_LOPC)

summary(m1) # r2 = 48.7%
anova(m1)
fPlotBiomassLM(m1, "Biomass_lm1")
# par(mfrow = c(2,2))
# plot(m1)

##### ***MODELS - LOPC = LMER (RANDOM EFFECT + NO GAMMA) #####
m1 <- lmer(log10(Biomass) ~ log10(Chl) + 
             fHarmonic(HarmTOD, k = 1) * exp(-Depth/1000) + ns(Bathy, df = 3) + 
             fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + (1|Project), data = dat_LOPC)
summary(m1)
r.squaredGLMM(m1) # Fixed r2 = 43.4%, Total r2 = 93.3%
fPlotBiomassLMERs(m1, "Biomass", "log10")
dev.copy2pdf(file = "BiomassModel_LOPC.pdf", paper = "A4r") #

saveRDS(m1, "LMER_LOPC.rds")

# Try LOPC data w/o DOY because not enough seasonal data
m1 <- lmer(log10(Biomass) ~ log10(Chl) + 
             fHarmonic(HarmTOD, k = 1) * exp(-Depth/1000) + ns(Bathy, df = 3) + 
             ns(SST, 3) + (1|Project), data = dat_LOPC)
summary(m1)
r.squaredGLMM(m1) # Fixed r2 = 26.8%, Total r2 = 61.6%
fPlotBiomassLMERs(m1, "Biomass", "log10")
dev.copy2pdf(file = "BiomassModel_LOPC_NoSeasonality.pdf", paper = "A4r") #

# DVM quite variable and difficult to see
x11(width = 14, height = 8)
ggplot(data = dat_LOPC, aes(x = Time, y = log10(Biomass))) +
  geom_point() +
  geom_smooth(method = "loess") +
  facet_wrap(vars(Project), scales = "free_y", ncol = 4)

plot(m1) # Pretty reasonable
x11(width = 8, height = 5, title = "Biomass Diagnostics")
par(mfrow = c(1,2))
qqnorm(residuals(m1))
qqline(residuals(m1), col = "blue")
plot(ranef(m1)$Project, main = "Project")

##### MODELS - LOPC = GLM (NO RANDOM EFFECT + GAMMA) #####
m1 <- glm(Biomass ~ log10(Chl) + 
            fHarmonic(HarmTOD, k = 1) * exp(-Depth/1000) + ns(Bathy, df = 3) + 
            fHarmonic(HarmDOY, k = 1) * ns(SST, 3), data = dat_LOPC,  maxit = 1000, family = Gamma(link = log))
summary(m1)
# anova(m1, test = "F")
1-(m1$deviance/m1$null.deviance) # r2 = 47.9%

# par(mfrow = c(2,2))
# plot(m1) # Residuals actually look better for the log10 transformation rather than log-link
fPlotBiomassLM(m1, "Biomass") # Surface bad

##### MODELS - LOPC = GLMER (RANDOM EFFECT + GAMMA) #####
m1 <- glmer(Biomass ~ log10(Chl) + 
              fHarmonic(HarmTOD, k = 1) * exp(-Depth/1000) + ns(Bathy, df = 3) + 
              fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + (1|Project), data = dat_LOPC, family = Gamma(link = log), nAGQ=0)
summary(m1)
r.squaredGLMM(m1)
fPlotBiomassLMERs(m1, "Biomass", "log")

plot(m1) # Very poor
x11(width = 8, height = 5, title = "Biomass Diagnostics")
par(mfrow = c(1,2))
qqnorm(residuals(m1))
qqline(residuals(m1), col = "blue")
plot(ranef(m1)$Route, main = "Route")


##### ***MODELS - ALL = LMER (RANDOM EFFECT + NO GAMMA) #####
# NOTE: If try interactions with NetType2 (Net, CPR, LOPC) then model is rank deficient so drops coefficients, and visreg won't plot them
# Just combine them
m1 <- lmer(log10(Biomass) ~ BiomassMethod + Tow + Mesh + log10(Chl) + 
             fHarmonic(HarmTOD, k = 1) * exp(-Depth/1000) + ns(Bathy, df = 3) + 
             fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + (1|Project), data = dat)
summary(m1)
r.squaredGLMM(m1) # Fixed r2 = 67.9%, Total r2 = 86.6%
fPlotBiomassLMERs(m1, "Biomass", "log10")
dev.copy2pdf(file = "BiomassModel_All.pdf", paper = "A4r") #
saveRDS(m1, "LMER_All.rds")

# Give Nets, CPR and LOPC equal weighting
Summary <- dat %>% group_by(NetType2) %>% 
  summarise(N = n()) 
dat <- dat %>% 
  mutate(WtVec = case_when(NetType2 == "Net" ~ 1 / Summary$N[3] * nrow(dat),
                           NetType2 == "CPR" ~ 1 / Summary$N[1] * nrow(dat),
                           NetType2 == "LOPC" ~ 1 / Summary$N[2] * nrow(dat)))
Mn <- mean(dat$WtVec)
dat$WtVec <- dat$WtVec / Mn
sum(dat$WtVec) # Check - correct!

m1 <- lmer(log10(Biomass) ~ BiomassMethod + Tow + Mesh + log10(Chl) + 
             fHarmonic(HarmTOD, k = 1) * exp(-Depth/1000) + ns(Bathy, df = 3) + 
             fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + (1|Project), weights = WtVec, data = dat)
summary(m1)
r.squaredGLMM(m1) # Fixed r2 = 69.4%, Total r2 = 87.5%
fPlotBiomassLMERs(m1, "Biomass", "log10")
dev.copy2pdf(file = "BiomassModel_All_EqualWt.pdf", paper = "A4r") #

x11(width = 8, height = 5, title = "Biomass Diagnostics")
par(mfrow = c(1,2))
plot(m1)
qqnorm(residuals(m1))
qqline(residuals(m1), col = "blue")
plot(ranef(m1)$Project, main = "Project")
