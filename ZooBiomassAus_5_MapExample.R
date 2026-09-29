# Biomass mapping around Australia
# Jase 22nd September 2026

##### PRELIMINARIES #####
library(tidyverse)
library(lme4)
library(terra)
library(sf)
library(splines) # needed so ns() is available at predict time
library(tidyterra) # for geom_spatraster
library(utils4mme)

# ── Load model ────────────────────────────────────────────────────────────────
# m1_map is a glmer (Gamma log-link) saved in ZooBiomassAus_2_StatisticalModels.R
# If you want to use the lmer (log10 Gaussian) models instead, swap the path:
#   "ModelOutput/Biomass_All_Final_EqualWeighting.rds"
model <- readRDS("Output/m1_map.rds")

# Confirm model class
stopifnot(inherits(model, c("lmerMod", "glmerMod")))
is_glmer <- inherits(model, "glmerMod") # TRUE  → Gamma/log-link (no back-transform needed beyond exp)
is_lmer <- inherits(model, "lmerMod") # TRUE  → log10 Gaussian (need 10^x back-transform)

# ── Load & filter data (used only to check variable names / levels) ───────────
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

# ── Load spatial rasters ──────────────────────────────────────────────────────

# SST (MODIS AQUA – pick a representative month)
sst <- terra::rast(
  "/Volumes/T9/MODIS_AQUA/SST/AQUA_MODIS.20210701_20210731.L3m.MO.SST.sst.4km.nc",
  subds = "sst"
) %>%
  terra::crop(bndry)
names(sst) <- "SST" # match model variable name exactly

# CHL (Copernicus)
chl <- terra::rast(
  "/Volumes/T9/Copernicus/OCEANCOLOUR_GLO_BGC_L4_MY_009_104/cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M_202603/2021/20210701-20210731_cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M.nc",
  subds = "CHL"
) %>%
  terra::resample(sst)
names(chl) <- "Chl" # match model variable name exactly

# Bathymetry (GEBCO tiles)
bathy_files <- list.files(
  "/Volumes/T9/gebco_2026_sub_ice_topo_geotiff/",
  "tif$",
  full.names = TRUE
)
ic <- terra::sprc(map(bathy_files, terra::rast))
rbathy <- terra::mosaic(ic, fun = "mean")
names(rbathy) <- "Bathy"
rbathy <- terra::resample(rbathy, sst)
rbathy <- terra::ifel(rbathy >= 0, NA, rbathy) # ocean only

# ── Latitude raster (needed by glmer: ns(Latitude, 2)) ───────────────────────
# Build a raster where every cell value equals its own latitude.
# This is the cleanest way to pass a spatially-varying Latitude to
# terra::predict() without relying on xy = TRUE (which is not a valid
# argument for terra::predict and causes "unused argument" errors).
lat_rast <- terra::init(sst, "y") # fills each cell with its y-coordinate
names(lat_rast) <- "Latitude"

# ── Stack rasters (only the variables that vary spatially) ────────────────────
# Variable names MUST match the column names used in the model formula.
rstack <- c(sst, chl, rbathy, lat_rast) # SST, Chl, Bathy, Latitude

# ── Constant (non-spatial) predictors ─────────────────────────────────────────
# Supply one row; terra::predict() will recycle it across all pixels.
# Adjust these to the "reference" conditions you want to map.
#
# glmer model formula:
#   Biomass ~ Method + Mesh + Tow + log(Chl) + ns(SST,3) + ns(Bathy,2) +
#             harmonic(HarmTOD,1) * ns(Depth,3) +
#             ns(Latitude,2) * harmonic(HarmDOY,2) +
#             (1|DatasetID) + (1|Gear)
#
# lmer model formula (Biomass_All_Final_EqualWeighting):
#   log10(Biomass) ~ BiomassMethod + log10(Chl) +
#                    fHarmonic(HarmTOD,1) * exp(-Depth/1000) +
#                    ns(Bathy,3) + fHarmonic(HarmDOY,1) * ns(SST,3) +
#                    (1|Gear_Mes2)

const_df <- data.frame(
  # ── shared / common ──────────────────────────────────────────────
  Depth = 10, # surface tow depth (m)
  HarmTOD = 1, # midday (12:00 → π radians)
  HarmDOY = pi, # mid-year (~July, π radians)
  Mesh = 100, # mesh size (µm) – typical plankton net
  Tow = "V", # tow type
  # ── glmer-specific ───────────────────────────────────────────────
  Method = "Wet", # BiomassMethod label used in glmer
  DatasetID = "IMOS", # random effect level (zeroed with re.form=NA)
  Gear = "1000", # random effect level (zeroed with re.form=NA)
  # ── lmer-specific ────────────────────────────────────────────────
  BiomassMethod = "Wet"
)

# ── Prediction function passed to terra::predict() ────────────────────────────
# terra::predict() calls fun(model, newdata) — it does NOT forward extra named
# arguments reliably. The solution is a closure: make_predict_fun() captures
# const_df and is_lmer in its enclosing environment and returns a two-argument
# function that terra can call without any extra args.

make_predict_fun <- function(const_df, is_lmer) {
  function(model, newdata) {
    # Merge per-pixel raster columns (SST, Chl, Bathy, Latitude) with constants
    nd <- cbind(newdata, const_df[rep(1, nrow(newdata)), , drop = FALSE])

    # Population-level prediction (re.form = NA → random effects = 0)
    preds <- predict(model, newdata = nd, re.form = NA, allow.new.levels = TRUE)

    # Back-transform to biomass units
    if (is_lmer) {
      return(10^preds) # lmer: log10(Biomass) → mg m⁻³
    } else {
      return(exp(preds)) # glmer Gamma(log): log scale → mg m⁻³
    }
  }
}

predict_fun <- make_predict_fun(const_df, is_lmer)

# ── Run spatial prediction ────────────────────────────────────────────────────
# terra::predict() passes a data.frame of raster values (one chunk at a time)
# to predict_fun(model, newdata). Latitude arrives as a raster column (lat_rast)
# so no xy argument is needed.

layer <- terra::predict(
  object = rstack,
  model = model,
  fun = predict_fun,
  na.rm = TRUE
)

names(layer) <- "Biomass_mg_m3"

# ── Save output ───────────────────────────────────────────────────────────────
terra::writeRaster(layer, "Output/Biomass_Map_Winter.tif", overwrite = TRUE)
saveRDS(layer, "Output/Biomass_Map_Winter.rds")

# ── Plot ──────────────────────────────────────────────────────────────────────
# Load Australian coastline
aus <- rnaturalearth::ne_countries()

imcra <- sf::st_read(
  "Data/Integrated_Marine_and_Coastal_Regionalisation_of_Australia_(IMCRA)_v4.0_-_Provincial_Bioregions.geojson"
)


ggplot() +
  tidyterra::geom_spatraster(data = layer, show.legend = FALSE) +
  scale_fill_viridis_c(
    option = "plasma",
    na.value = "transparent",
    transform = "log10",
    limits = c(10, 1000),
    oob = scales::squish
  ) +
  geom_sf(
    data = imcra,
    colour = "black",
    fill = NA,
    show.legend = FALSE,
    linewidth = 0.5
  ) +
  geom_sf(data = aus, fill = "grey50", colour = "grey80") +
  coord_sf(xlim = c(100, 170), ylim = c(-50, -8), expand = FALSE) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw()

ggsave(
  "Figures/Biomass_Map_Winter.png",
  dpi = 600,
  height = 6,
  width = 8
)

## NOW DO SUMMER

# SST (MODIS AQUA – pick a representative month)
sst <- terra::rast(
  "/Volumes/T9/MODIS_AQUA/SST/AQUA_MODIS.20210101_20210131.L3m.MO.SST.sst.4km.nc",
  subds = "sst"
) %>%
  terra::crop(bndry)
names(sst) <- "SST" # match model variable name exactly

# CHL (Copernicus)
chl <- terra::rast(
  "/Volumes/T9/Copernicus/OCEANCOLOUR_GLO_BGC_L4_MY_009_104/cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M_202603/2021/20210101-20210131_cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M.nc",
  subds = "CHL"
) %>%
  terra::resample(sst)
names(chl) <- "Chl" # match model variable name exactly

rstack <- c(sst, chl, rbathy, lat_rast) # SST, Chl, Bathy, Latitude

const_df <- data.frame(
  # ── shared / common ──────────────────────────────────────────────
  Depth = 10, # surface tow depth (m)
  HarmTOD = 1, # midday (12:00 → π radians)
  HarmDOY = 0.017, # mid-year (~July, π radians)
  Mesh = 100, # mesh size (µm) – typical plankton net
  Tow = "V", # tow type
  # ── glmer-specific ───────────────────────────────────────────────
  Method = "Wet", # BiomassMethod label used in glmer
  DatasetID = "IMOS", # random effect level (zeroed with re.form=NA)
  Gear = "1000", # random effect level (zeroed with re.form=NA)
  # ── lmer-specific ────────────────────────────────────────────────
  BiomassMethod = "Wet"
)


layer <- terra::predict(
  object = rstack,
  model = model,
  fun = predict_fun,
  na.rm = TRUE
)

names(layer) <- "Biomass_mg_m3"

# ── Save output ───────────────────────────────────────────────────────────────
terra::writeRaster(layer, "Output/Biomass_Map_Summer.tif", overwrite = TRUE)
saveRDS(layer, "Output/Biomass_Map_Summer.rds")


ggplot() +
  tidyterra::geom_spatraster(data = layer, show.legend = FALSE) +
  scale_fill_viridis_c(
    option = "plasma",
    na.value = "transparent",
    transform = "log10",
    limits = c(10, 1000),
    oob = scales::squish
  ) +
  geom_sf(
    data = imcra,
    colour = "black",
    fill = NA,
    show.legend = FALSE,
    linewidth = 0.5
  ) +
  geom_sf(data = aus, fill = "grey50", colour = "grey80") +
  coord_sf(xlim = c(100, 170), ylim = c(-50, -8), expand = FALSE) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw()

ggsave(
  "Figures/Biomass_Map_Summer.png",
  dpi = 600,
  height = 6,
  width = 8
)


# MOCK UP BAYSEIAN OUTPUT

ggplot() +
  tidyterra::geom_spatraster(data = layer, show.legend = FALSE) +
  scale_fill_viridis_c(
    option = "plasma",
    na.value = "transparent",
    transform = "log10",
    limits = c(1, 5000),
    oob = scales::squish
  ) +
  geom_sf(
    data = imcra,
    colour = "black",
    fill = NA,
    show.legend = FALSE,
    linewidth = 0.5
  ) +
  geom_sf(data = aus, fill = "grey50", colour = "grey80") +
  coord_sf(xlim = c(100, 170), ylim = c(-50, -8), expand = FALSE) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw()

ggsave(
  "Figures/Biomass_Map_Future.png",
  dpi = 600,
  height = 6,
  width = 8
)
