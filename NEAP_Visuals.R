library(tidyverse)
library(tidyterra)
library(sf)


dat1 <- read_delim(
  "Data/BiomassDatabase20190620.csv",
  delim = "\t",
  na = "(null)"
) %>%
  filter(NET_TYPE != "Optical Plankton Counter (OPC-2T)") %>% # Remove because I re-add it again with more data below.
  select(LATITUDE, LONGITUDE) %>%
  rename(Latitude = LATITUDE, Longitude = LONGITUDE)


dat2 <- read_csv("Data/MNF_LOPC_Biomass.csv", na = "NaN") %>%
  select(LATITUDE, LONGITUDE) %>%
  rename(Latitude = LATITUDE, Longitude = LONGITUDE)

dat3 <- read_rds("../GlobalZoopBiomass/Output/dat_cleaned.rds") %>%
  select(Latitude, Longitude)


dat <- bind_rows(dat1, dat2, dat3) %>%
  distinct(Latitude, Longitude) %>%
  sf::st_as_sf(
    coords = c("Longitude", "Latitude"),
    remove = FALSE,
    crs = "EPSG:4326"
  )

aus <- rnaturalearth::ne_countries(country = "Australia")

imcra <- sf::st_read(
  "Data/Integrated_Marine_and_Coastal_Regionalisation_of_Australia_(IMCRA)_v4.0_-_Provincial_Bioregions.geojson"
)

# bndry <- imcra %>%
#   sf::st_buffer(dist = 10) %>%
#   summarise(geometry = sf::st_union(geometry)) %>%
#   sf::st_make_valid() %>%
#   sf::st_bbox() %>%
#   sf::st_as_sfc()

bndry <- sf::st_bbox(
  c(xmin = 100, xmax = 170, ymax = -8, ymin = -50),
  crs = sf::st_crs("EPSG:4326")
) %>%
  sf::st_as_sfc()

dat_aus <- sf::st_filter(dat, bndry)

(gg <- ggplot() +
  geom_sf(data = aus, fill = "grey50", colour = "grey80") +
  geom_sf(data = imcra, colour = "red", fill = NA, show.legend = FALSE) +
  geom_sf(data = dat_aus, colour = "blue", size = 0.1) +
  theme_bw() +
  coord_sf(xlim = c(100, 170), ylim = c(-50, -8), expand = FALSE))

ggsave(
  "Figures/SamplingLocations.png",
  plot = gg,
  dpi = 600,
  height = 6,
  width = 8
)


sst <- terra::rast(
  "/Volumes/T9/MODIS_AQUA/SST/AQUA_MODIS.20060501_20060531.L3m.MO.SST.sst.4km.nc",
  subds = "sst"
) %>%
  terra::crop(bndry)

(gg_sst <- ggplot() +
  geom_spatraster(data = sst, aes(fill = sst), show.legend = FALSE) +
  geom_sf(
    data = imcra,
    colour = "black",
    fill = NA,
    show.legend = FALSE,
    linewidth = 0.5
  ) +
  geom_sf(data = aus, fill = "grey50", colour = "grey80") +
  theme_bw() +
  coord_sf(xlim = c(100, 170), ylim = c(-50, -8), expand = FALSE) +
  cmocean::scale_fill_cmocean(name = "thermal"))


ggsave(
  "Figures/SST.png",
  plot = gg_sst,
  dpi = 600,
  height = 6,
  width = 8
)


chl <- terra::rast(
  "/Volumes/T9/Copernicus/OCEANCOLOUR_GLO_BGC_L4_MY_009_104/cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M_202603/2021/20210501-20210531_cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M.nc",
  subds = "CHL"
) %>%
  terra::crop(bndry)

(gg_chl <- ggplot() +
  geom_spatraster(data = chl, aes(fill = CHL), show.legend = FALSE) +
  geom_sf(
    data = imcra,
    colour = "black",
    fill = NA,
    show.legend = FALSE,
    linewidth = 0.5
  ) +
  geom_sf(data = aus, fill = "grey50", colour = "grey80") +
  theme_bw() +
  coord_sf(xlim = c(100, 170), ylim = c(-50, -8), expand = FALSE) +
  cmocean::scale_fill_cmocean(name = "algae", transform = "log10"))


ggsave(
  "Figures/Chl.png",
  plot = gg_chl,
  dpi = 600,
  height = 6,
  width = 8
)


bathy_files <- list.files(
  "/Volumes/T9/gebco_2026_sub_ice_topo_geotiff/",
  "tif$",
  full.names = TRUE
)

# Load bathy data
ic <- terra::sprc(map(bathy_files, terra::rast))

rbathy <- terra::mosaic(ic, fun = "mean")
terra::set.names(rbathy, "Bathy") # rename layer

rbathy <- rbathy %>%
  terra::crop(bndry)

rbathy <- terra::ifel(rbathy >= 0, NA, rbathy)

(gg_bathy <- ggplot() +
  geom_spatraster(data = rbathy, aes(fill = Bathy), show.legend = FALSE) +
  geom_sf(
    data = imcra,
    colour = "black",
    fill = NA,
    show.legend = FALSE,
    ,
    linewidth = 0.5
  ) +
  geom_sf(data = aus, fill = "grey50", colour = "grey80") +
  theme_bw() +
  coord_sf(xlim = c(100, 170), ylim = c(-50, -8), expand = FALSE) +
  cmocean::scale_fill_cmocean(name = "deep", direction = -1))


ggsave(
  "Figures/Bathy.png",
  plot = gg_bathy,
  dpi = 600,
  height = 6,
  width = 8
)
