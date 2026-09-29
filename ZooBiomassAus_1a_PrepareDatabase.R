# Biomass analysis around Australia
# Ant & Jase 1/7/2019
# Data provided by Felicity 20/6/2019
# Jase Updated 10/10/2019
# Jase Updated 21/09/2026

library(tidyverse)
# library(janitor) # install.packages("janitor")

### MOVE ALL THIS TO GLOBALZOOBIOMASS SO WE CREATE EVERYTING IN ONE PLACE.
# WE THEN FILTER FOR AUSTRALIA OR DIFFERENT DATATYPES WHEN WE GO TO USE
# THE DATA FOR STATISTICAL MODELLING

dat1 <- read_delim(
  "Data/BiomassDatabase20190620.csv",
  delim = "\t",
  na = "(null)"
) %>%
  filter(NET_TYPE != "Optical Plankton Counter (OPC-2T)") %>% # Remove because I re-add it again with more data below.
  mutate(PROJECT_ID = as.character(PROJECT_ID))

dat2 <- read_csv("Data/MNF_LOPC_Biomass.csv", na = "NaN") %>%
  mutate(
    SAMPLE_DATE_LOCAL = dmy_hms(SAMPLE_DATE_LOCAL),
    SAMPLE_DATE_UTC = dmy_hms(SAMPLE_DATE_UTC),
    RECORD_ID = as.character(RECORD_ID),
    SAMPLE_ID = as.character(SAMPLE_ID),
    NET_MESH_SIZE = 300
  ) %>%
  mutate(
    DAY_NIGHT = case_when(
      hour(SAMPLE_DATE_LOCAL) > 6 & hour(SAMPLE_DATE_LOCAL) < 18 ~ "day",
      hour(SAMPLE_DATE_LOCAL) <= 6 | hour(SAMPLE_DATE_LOCAL) >= 18 ~ "night"
    )
  )

dat3 <- read_rds("../GlobalZoopBiomass/Output/dat_cleaned.rds")


dat <- bind_rows(dat1, dat2)

dat <- janitor::clean_names(dat, case = "upper_camel")

dat <- dat %>%
  rename(
    DateLocal = SampleDateLocal,
    DateUTC = SampleDateUtc,
    Year = SampleYear,
    Month = SampleMonth,
    Day = SampleDay
  )

# Function which adds satellite SST, Chl and Bathy data to any data frame with "Latitude" and "Longitude" variables
dat <- fAddEnviro(dat, EnviroDir, bathy_file)

summary(dat$SST)
summary(dat$Chl)
summary(dat$Bathy)

dat <- dat %>%
  rename(
    Mesh = NetMeshSize,
    Project = ProjectId,
    SampleID = SampleId,
    RecordID = RecordId,
    Tow = TowType
  ) %>%
  mutate(
    Chl = replace(Chl, Chl > 10, 10),
    DOY = yday(DateLocal),
    Bathy = replace(Bathy, Bathy > -20, -20),
    Bathy = abs(Bathy), # Make Bathymetry positive
    MinDepth = replace(MinDepth, MinDepth == -999, NA),
    MaxDepth = replace(MaxDepth, MaxDepth == -999, NA),
    Depth = (MinDepth + MaxDepth) / 2, # There are a few where the depths are reversed due to the sampling design (upcast of MOCNESS)
    BiomassMethod = str_replace(BiomassMethod, "wet mass", "Wet"),
    BiomassMethod = str_replace(BiomassMethod, "dry mass", "Dry"),
    BiomassMethod = str_replace(
      BiomassMethod,
      "displacement volume",
      "Displacement"
    ),
    BiomassMethod = str_replace(BiomassMethod, "settled volume", "Settled"),
    BiomassMethod = str_replace(
      BiomassMethod,
      "Estimated Wet Weight from BioVolume",
      "EstimatedWetWeight"
    ),
    BiomassMethod = as.factor(BiomassMethod),
    Tow = str_replace(Tow, "horizontal tow", "Horizontal"),
    Tow = str_replace(Tow, "double oblique tow", "Oblique"),
    Tow = str_replace(Tow, "oblique tow", "Oblique"),
    Tow = str_replace(Tow, "vertical", "Vertical"),
    Tow = str_replace(Tow, "unknown", "Unknown"),
    Tow = as.factor(Tow),
    DayNight = str_to_title(DayNight),
    DayNight = as.factor(DayNight),
    Project = as.factor(Project)
  ) %>%
  filter(!is.na(SST)) %>%
  filter(!is.na(Chl)) %>% # Remove NAs
  filter(Biomass >= 0) %>% # Remove negative biomass
  filter(!is.na(Depth))

summary(dat) # There are zeros in the Biomass data
sum(dat$Biomass == 0) # 212 from the LOPC - delete

dat <- dat %>%
  mutate(
    NetType2 = case_when(
      NetType == "LOPC" ~ "LOPC",
      NetType == "Continuous plankton recorder" ~ "CPR",
      TRUE ~ "Net"
    ),
    NetType2 = as.factor(NetType2)
  ) %>%
  filter(Biomass > 0) %>% # n = 54,033
  mutate(
    Tow = as.character(Tow),
    Tow = replace(Tow, NetType == "Continuous plankton recorder", "Surface"),
    Tow = as.factor(Tow)
  )

# Sort out Project for CPRs
Routes <- read_csv("Data/cpr_route.csv")
dat <- left_join(dat, Routes, by = c("SampleID" = "SAMPLE_ID"))
rm(Routes)

# Move the CPR route to the Project column
dat$Project <- as.character(dat$Project)
dat$Project[dat$Project == "597"] <- dat$ROUTE[dat$Project == "597"]

dat <- dat %>%
  dplyr::select(-ROUTE) %>%
  mutate(Project = as.factor(Project))

saveRDS(dat, file = "Output/BiomassData.rds")

write_csv(dat, "Output/FinalBiomassData.csv")
