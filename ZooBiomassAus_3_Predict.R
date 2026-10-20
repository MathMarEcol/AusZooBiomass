# Predict m13 across spatial grid using terra SpatRasters ----
# Chl: CMEMS OCEANCOLOUR_GLO_BGC_L4_MY_009_104 climatological monthly means
# SST: MODIS-Aqua L3m 4 km monthly composites, climatological monthly means
# Bathy: GEBCO 2026, resampled to Chl grid

library(terra)
library(glmmTMB)
library(splines)
library(utils4mme) # harmonic()
library(lubridate)
library(tidyverse)
library(ncdf4)
library(tidyterra)

reprocess <- FALSE

# Month definitions ----
month_names <- c(
  "01January",
  "02February",
  "03March",
  "04April",
  "05May",
  "06June",
  "07July",
  "08August",
  "09September",
  "10October",
  "11November",
  "12December"
)

# Study region boundary ----
bndry <- sf::st_bbox(
  c(xmin = 90, xmax = 179, ymax = 0, ymin = -55),
  crs = sf::st_crs("EPSG:4326")
) %>%
  sf::st_as_sfc()


if (isTRUE(reprocess)) {
  # Load Chl Data ----
  # Source: CMEMS OCEANCOLOUR_GLO_BGC_L4_MY_009_104, 4 km monthly composites
  # Variable: CHL (mg m-3)
  # Strategy: climatological monthly mean — average all available years per calendar month

  chl_base <- "/Volumes/T9/Copernicus/OCEANCOLOUR_GLO_BGC_L4_MY_009_104/cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M_202603"

  chl_year_range <- 2003:2015 # adjust as needed

  chl_files <- list.files(
    chl_base,
    pattern = "\\.nc$",
    full.names = TRUE,
    recursive = TRUE
  )
  chl_years <- as.integer(substr(basename(chl_files), 1, 4)) # chars 1-4 of YYYYMMDD... = YYYY
  chl_files <- chl_files[chl_years %in% chl_year_range]
  chl_months <- as.integer(substr(basename(chl_files), 5, 6)) # chars 5-6 of YYYYMMDD... = MM

  chl_clim <- vector("list", 12)

  for (m in 1:12) {
    files_m <- chl_files[chl_months == m]

    r_m <- terra::rast(
      lapply(files_m, function(f) {
        terra::rast(f, subds = "CHL") %>%
          terra::crop(bndry)
      })
    )

    chl_clim[[m]] <- terra::app(r_m, mean, na.rm = TRUE)
  }

  # 12-layer SpatRaster: one climatological mean per calendar month
  rchl <- terra::rast(chl_clim)
  names(rchl) <- gsub("\\d", "", month_names)
  terra::writeRaster(rchl, "Output/rchl_clim.tif", overwrite = TRUE)

  # Load SST Data ----
  # Source: MODIS-Aqua L3m 4 km monthly composites (NASA OceanColor)
  # Variable: sst (degrees C)
  # Files: flat directory, named AQUA_MODIS.YYYYMMDD_YYYYMMDD.L3m.MO.SST.sst.4km.nc
  # Strategy: climatological monthly mean — average all available years per calendar month
  # Note: year range matched to chl_year_range for consistency

  sst_base <- "/Volumes/T9/MODIS_AQUA/SST"

  # Exclude macOS metadata files (._AQUA_MODIS...) with ^AQUA_MODIS pattern
  sst_files <- list.files(
    sst_base,
    pattern = "^AQUA_MODIS.*\\.nc$",
    full.names = TRUE
  )
  sst_years <- as.integer(substr(basename(sst_files), 12, 15)) # chars 12-15 = YYYY
  sst_files <- sst_files[sst_years %in% chl_year_range]
  sst_months <- as.integer(substr(basename(sst_files), 16, 17)) # chars 16-17 = MM

  sst_clim <- vector("list", 12)

  for (m in 1:12) {
    files_m <- sst_files[sst_months == m]

    r_m <- terra::rast(
      lapply(files_m, function(f) {
        terra::rast(f, subds = "sst") %>%
          terra::crop(bndry)
      })
    )

    sst_clim[[m]] <- terra::app(r_m, mean, na.rm = TRUE)
  }

  # 12-layer SpatRaster: one climatological mean per calendar month
  rsst <- terra::rast(sst_clim)
  names(rsst) <- gsub("\\d", "", month_names)
  terra::writeRaster(rsst, "Output/rsst_clim.tif", overwrite = TRUE)

  # Load Bathy Data ----
  # Source: GEBCO 2026 sub-ice topo GeoTIFF tiles
  # Mosaicked and land-masked here; resampled to Chl grid inside the month loop
  # (resample must happen inside the loop because Chl defines the target grid)

  bathy_files <- list.files(
    "/Volumes/T9/gebco_2026_sub_ice_topo_geotiff/",
    "tif$",
    full.names = TRUE
  )

  ic <- terra::sprc(map(bathy_files, terra::rast))
  rbathy <- terra::mosaic(ic, fun = "mean")
  rbathy <- terra::crop(rbathy, bndry)
  rbathy <- terra::ifel(rbathy >= 0, NA, rbathy)
  terra::set.names(rbathy, "Bathy")
  terra::writeRaster(rbathy, "Output/rbathy.tif", overwrite = TRUE)
} else {
  rsst <- terra::rast("Output/rsst_clim.tif")
  rchl <- terra::rast("Output/rchl_clim.tif")
  rbathy <- terra::rast("Output/rbathy.tif")
}

# Load model and training data ----
m13 <- readRDS("Output/m13.rds")
dat <- readRDS("Output/FinalBiomassData.rds")


# Fixed constants ----
# Held constant across all grid cells and months.
# Method and Tow are stored as character here; they are converted to factors
# with the correct levels inside pred_fun.  We pull levels from the fitted
# model frame (not from dat) because dat may not have factor levels stored
# correctly after being saved/loaded from .rds.
# Note: Tow levels in the model are abbreviated ("H", "O", "V").
fixed_method <- "Wet"
fixed_tow <- "O" # "O" = Oblique in the model's Tow factor
fixed_harm_tod <- 0 # Midnight
fixed_depth_decay <- exp(-10 / 1000) # 10 m depth

# Pull factor levels from the fitted model frame ----
# model.frame(m13) is the only reliable source of the exact levels glmmTMB
# used at fit time — dat loaded from .rds may have lost factor levels.
mf <- model.frame(m13)
method_levels <- levels(mf$Method)
tow_levels <- levels(mf$Tow)

# Predict wrapper ----
# terra::predict() extracts raster values into a plain data.frame and calls
# fun(model, newdata, ...).  We rebuild newdata here so that:
#   - Method and Tow are proper factors with the correct levels from the model
#   - re.form = NA suppresses random effects (population-level prediction)
#   - type = "response" back-transforms from the log link
pred_fun <- function(model, newdata, ...) {
  newdata$Method <- factor(fixed_method, levels = method_levels)
  newdata$Tow <- factor(fixed_tow, levels = tow_levels)
  newdata$HarmTOD <- fixed_harm_tod
  newdata$Depth_decay <- fixed_depth_decay
  predict(model, newdata = newdata, re.form = NA, type = "response")
}

# Month loop: build raster stack and predict ----
pred_list <- vector("list", length(month_names))

for (i in seq_along(month_names)) {
  # 1. Select this month's climatological Chl layer
  Chl <- rchl[[i]]
  names(Chl) <- "Chl"

  # 2. Resample Bathy to match Chl grid (GEBCO ~500 m -> Chl 4 km)
  Bathy <- terra::resample(rbathy, Chl, method = "bilinear")
  names(Bathy) <- "Bathy"

  # 3. DOY for the 15th of this month (reference year 2010)
  doy <- lubridate::yday(
    lubridate::parse_date_time(
      paste0("15 ", gsub("\\d", "", month_names[i]), " 2010"),
      "dmy"
    )
  )
  harm_doy <- (doy / 365) * 2 * pi

  # 4. Select this month's climatological SST layer, resampled to Chl grid
  SST <- terra::resample(rsst[[i]], Chl, method = "bilinear")
  names(SST) <- "SST"

  # 5. HarmDOY: scalar this month, broadcast to every cell as a raster layer
  #    (must be a layer so terra extracts it as a column in newdata)
  HarmDOY_r <- terra::setValues(Chl, harm_doy)
  names(HarmDOY_r) <- "HarmDOY"

  # 6. Stack and predict via pred_fun wrapper
  r_stack <- c(Chl, Bathy, SST, HarmDOY_r)

  pred_list[[i]] <- terra::predict(r_stack, m13, fun = pred_fun)
}

# Stack and summarise ----
pred_stack <- terra::rast(pred_list)
names(pred_stack) <- gsub("\\d", "", month_names) # "January", "February", ...

pred_annual <- terra::app(pred_stack, mean, na.rm = TRUE)


# Attach time dimension ----
month_dates <- lubridate::as_date(
  lubridate::parse_date_time(
    paste0("15 ", gsub("\\d", "", month_names), " 2010"),
    "dmy"
  )
)
terra::time(pred_stack) <- month_dates


write_rds(pred_stack, "Output/Biomass_m13_monthly.rds")
write_rds(pred_annual, "Output/Biomass_m13_annual.rds")


# Write monthly stack ----
terra::writeCDF(
  pred_stack,
  filename = "Output/Biomass_m13_monthly.nc",
  varname = "zooplankton_biomass",
  longname = "Zooplankton wet biomass predicted by m13 (glmmTMB Gamma-log)",
  unit = "mg m-3",
  zname = "time",
  compression = 5,
  overwrite = TRUE
)

# Write annual mean ----
terra::writeCDF(
  pred_annual,
  filename = "Output/Biomass_m13_annual.nc",
  varname = "zooplankton_biomass_annual",
  longname = "Annual mean zooplankton wet biomass predicted by m13 (glmmTMB Gamma-log)",
  unit = "mg m-3",
  compression = 5,
  overwrite = TRUE
)


# Append global attributes ----
nc <- ncdf4::nc_open("Output/Biomass_m13_monthly.nc", write = TRUE)
ncdf4::ncatt_put(
  nc,
  0,
  "title",
  "Australian zooplankton biomass — m13 predictions"
)

ncdf4::ncatt_put(
  nc,
  0,
  "author",
  "Jason Everett (Everdat Solutions Pty Ltd & The University of Queensland) and Claire Davies (CSIRO)"
)
ncdf4::ncatt_put(
  nc,
  0,
  "model",
  "m13: glmmTMB Gamma(log), fixed effects only"
)
ncdf4::ncatt_put(
  nc,
  0,
  "fixed_predictors",
  "Method=Wet, Tow=Oblique, HarmTOD=0, Depth=10m"
)
ncdf4::ncatt_put(nc, 0, "creation_date", as.character(Sys.Date()))
ncdf4::ncatt_put(nc, 0, "Conventions", "CF-1.8")
ncdf4::nc_close(nc)


## Plot output ----

aus <- rnaturalearth::ne_countries(scale = "medium") %>%
  sf::st_make_valid() %>%
  sf::st_crop(bndry)

(gg_ann <- ggplot() +
  geom_spatraster(data = pred_annual, aes(fill = mean)) +
  geom_sf(data = aus, fill = "grey50", colour = "grey80") +
  theme_bw() +
  coord_sf(xlim = c(100, 170), ylim = c(-50, -8), expand = FALSE) +
  cmocean::scale_fill_cmocean(name = "thermal", transform = "log10"))

ggsave(
  "Figures/pred_annual.png",
  plot = gg_ann,
  width = 20,
  height = 12,
  dpi = 1000
)

ggsave(
  "Figures/pred_annual.pdf",
  plot = gg_ann,
  width = 20,
  height = 12
)


(gg_seas <- ggplot() +
  geom_spatraster(data = pred_stack) +
  geom_sf(data = aus, fill = "grey50", colour = "grey80") +
  theme_bw() +
  coord_sf(xlim = c(100, 170), ylim = c(-50, -8), expand = FALSE) +
  facet_wrap(~lyr) +
  cmocean::scale_fill_cmocean(name = "thermal", transform = "log10"))

ggsave(
  "Figures/pred_season.png",
  plot = gg_seas,
  width = 20,
  height = 12,
  dpi = 1000
)

ggsave(
  "Figures/pred_season.pdf",
  plot = gg_seas,
  width = 20,
  height = 12
)
