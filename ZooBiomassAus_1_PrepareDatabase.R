# Biomass analysis around Australia
# Ant & Jase 1/7/2019
# Data provided by Felicity 20/6/2019
# Jase Updated 10/10/2019
# Jase Updated 21/09/2026

library(tidyverse)


# Load the data: n = 257,474 ----
dat <- read_rds("../GlobalZoopBiomass/Output/dat_cleaned.rds") %>%
  sf::st_as_sf(
    coords = c("Longitude", "Latitude"),
    remove = FALSE,
    crs = "EPSG:4326"
  )

bndry <- sf::st_bbox(
  c(xmin = 80, xmax = 179, ymax = 0, ymin = -60),
  crs = sf::st_crs("EPSG:4326")
) %>%
  sf::st_as_sfc()

## Filter for Australia! (n = 51,412)
dat <- sf::st_filter(dat, bndry)
# 50K

ggplot(data = dat, aes(x = Biomass)) +
  geom_histogram() +
  facet_wrap(vars(Method), scales = "free") +
  scale_y_log10()


## Check LOPC again

dat <- dat %>%
  mutate(
    Chl = if_else(Chl > 5, 5, Chl),
    Depth = if_else(Depth > 500, 500, Depth),
    Bathy = if_else(Bathy > 6000, 6000, Bathy),
    Biomass = if_else(Biomass > 2000, 2000, Biomass) # Max 2g m-3
  ) %>%
  filter(!Project %in% c("1099", "2", "21"))

# Check the data ----

## Biomass methods
unique(dat$Method)
# "Dry"     "CPR"     "Wet"     "Disp"    "Settled" "LOPC"

unique(dat$Units)
# mg/m3            ml/m3

# dat_ln <- dat %>% filter(Units == "ml/m3 log normal")

# Check combinations
dat %>%
  sf::st_drop_geometry() %>%
  dplyr::count(Method, Units)

# Check Tow Type
unique(dat$Tow)

hist(log10(dat$Biomass))
hist(dat$TimeLocal)
hist(dat$Mesh)
hist(log10(dat$Volume))

hist(dat$Depth)

# Check environmental data

#### SST
summary(dat$SST) # No NAs
hist(dat$SST)
ggplot(data = dat, aes(x = SST, y = Latitude)) +
  geom_point()


#### Chl
summary(dat$Chl) # 1 NA - Code needs a rewrite. Not worth the time for 1 point. Just remove below. It looks to be on land anyway.
hist(dat$Chl)
hist(log10(dat$Chl))


ggplot(data = dat, aes(x = Chl, y = Latitude)) +
  geom_point()

ggplot(
  data = dat,
  aes(x = log10(Chl), y = lubridate::month(Date))
) +
  geom_point()
# TODO JDE - Not convinced that the Chl in the southern hemisphere is correct, but it looks good in the North.

dat <- dat %>%
  filter(!is.na(Chl)) # Remove NAs

summary(dat$Chl)

ggplot() +
  geom_point(
    data = dat,
    aes(x = Longitude, y = Latitude, fill = log10(Chl), colour = log10(Chl))
  ) +
  cmocean::scale_fill_cmocean(
    name = "algae",
    aesthetics = c("colour", "fill"),
    limits = c(-2, 0),
    oob = scales::squish
  ) +
  geom_sf(
    data = rnaturalearth::ne_countries(returnclass = "sf"),
    fill = "grey50"
  ) + #world coastline
  theme_bw() +
  coord_sf(
    xlim = sf::st_bbox(bndry)[c("xmin", "xmax")],
    ylim = sf::st_bbox(bndry)[c("ymin", "ymax")]
  )


### Bathymetry

summary(dat$Bathy)
hist(dat$Bathy)

# Plot the data to check
ggplot() +
  geom_point(
    data = dat,
    aes(x = Longitude, y = Latitude, colour = log10(Bathy))
  ) +
  cmocean::scale_fill_cmocean(name = "deep", aesthetics = c("colour", "fill")) +
  geom_sf(
    data = rnaturalearth::ne_countries(returnclass = "sf"),
    fill = "grey50"
  ) + #world coastline
  theme_bw() +
  coord_sf(
    xlim = sf::st_bbox(bndry)[c("xmin", "xmax")],
    ylim = sf::st_bbox(bndry)[c("ymin", "ymax")]
  )


## Plot Longhurst Provinces to check
ggplot() +
  geom_point(
    data = dat,
    aes(x = Longitude, y = Latitude, colour = Longhurst_ProvCode),
    show.legend = FALSE
  ) +
  geom_sf(
    data = rnaturalearth::ne_countries(returnclass = "sf"),
    fill = "grey50"
  ) + #world coastline
  theme_bw() +
  coord_sf(
    xlim = sf::st_bbox(bndry)[c("xmin", "xmax")],
    ylim = sf::st_bbox(bndry)[c("ymin", "ymax")]
  )


##

write_rds(
  x = dat %>% sf::st_drop_geometry(),
  file = "Output/FinalBiomassData.rds"
)
