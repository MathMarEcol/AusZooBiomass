# Biomass analysis around Australia
# Ant & Jase 1/7/2019
# Data provided by Felicity 20/6/2019

library(tidyverse)
library(lubridate)
library(janitor) # install.packages("janitor")

library(maptools)
data(wrld_simpl)

dat <- readRDS("BiomassData.rds")

library(sf)
shapef <- "../CommonResources/WorldBorders/TM_WORLD_BORDERS-0.3.shp"

wb_sf <- read_sf(shapef)
wb_sp <- as(wb_sf, "Spatial")
wb_df <- fortify(wb_sp)  
wb_df$group <- as.numeric(wb_df$group)
wb_df$hole <- as.numeric(wb_df$hole)

# Define themes for plots 
theme_opts <- list(theme(panel.grid.minor = element_blank(),
                         panel.grid.major = element_blank(),
                         panel.background = element_blank(),
                         plot.background = element_rect(fill="white"),
                         panel.border = element_blank(),
                         axis.line = element_blank(),
                         axis.text.x = element_blank(),
                         axis.text.y = element_blank(),
                         axis.ticks = element_blank(),
                         axis.title.x = element_blank(),
                         axis.title.y = element_blank(),
                         plot.title = element_text(size=22)))

latlonCRS <- CRS("+proj=longlat +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m")

##################################
## Plot data distribution

gg <- ggplot() + 
  geom_polygon(data=wb_df, mapping=aes(x=long, y=lat, group = group)) +
  geom_point(data = dat, aes(x = Longitude, y = Latitude, colour=Project), size=1) +
  theme_bw() +
  theme_opts +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  coord_map(projection = "mercator") +
  scale_x_continuous(expand = c(0, 0), limits = c(100, 160)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-55, 1))

ggsave("Figures/BiomassDataDistribution_Project.png", width = 297, height = 210, 
       units = "mm", dpi=300)



gg <- ggplot() + 
  geom_polygon(data=wb_df, mapping=aes(x=long, y=lat, group = group)) +
  geom_point(data = dat, aes(x = Longitude, y = Latitude, colour=NetType2), size=1) +
  theme_bw() +
  theme_opts +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  coord_map(projection = "mercator") +
  scale_x_continuous(expand = c(0, 0), limits = c(100, 160)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-60, 1))

ggsave("Figures/BiomassDataDistribution_GearType.png", width = 297, height = 210, 
       units = "mm", dpi=300)


ggfw <- gg +
  facet_wrap(vars(BiomassMethod), ncol = 3)

ggsave("Figures/BiomassDataDistribution_FacetWrap.png", width = 297, height = 210, 
       units = "mm", dpi=300)



