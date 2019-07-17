library(tidyverse)
library(raster)
library(rasterImage)
library(sf)
library(ggpubr)
library(splines)
library(gganimate)

source("fHarmonic.R")
source("fPlotAusLMER.R")
source("fBuildRasterBrick.R")
source("fWorldBorders.R")


# Set the max values to be removed in rasterbrick
maxm <- tibble(Bathy = 6000, SST = 28, Chl = 10)
res <- "one"
####################################################
NthMonth <- c("01January", "02February","03March","04April","05May","06June",
              "07July","08August","09September","10October","11November","12December")

df <- data.frame(HarmTOD = 1, BiomassMethod = "Wet", 
                 Project = "1090", Depth = 10, Gear_Mes2 = "Bongo net_150")

# These are zeroed out in the final model, but they need to be included (in df) initially
# excl <- c("s(Project)") 

##############################################################################

# Predict abundance using the predict function in raster and plot
tframe <- "month"



stat_brick <- readRDS("Biomass_Layer_Joint.rds")


data <- projectRaster(stat_brick, crs = CRS, over = TRUE)

# Load the World Borders
wb_df <- fWorldBorders()

# Convert data to dataframe
data_df <- as.data.frame(rasterToPoints(data))




data_df_long <- gather(data_df, key = "Month", value = "Biomass", January:December)
data_df_long <- data_df_long %>% 
  mutate(Month = as.factor(Month), ordered==TRUE)

data_df_long <- data_df_long %>% 
  mutate(Month2 = match(data_df_long$Month, month.name))

gg <- ggplot() +
  geom_tile(data = data_df_long, aes(x = x, y = y, fill = Biomass))  +
  # geom_polygon(fill = "white") +
  geom_polygon(data = wb_df, aes(x = long, y = lat, group = group)) +
  scale_fill_gradientn(name =  expression(paste("mg m"^-3)), 
                       colours = rev(col),
                       limits = c(0.8, 3),
                       # breaks = cticks,
                       # labels = clabel,
                       position = "bottom",
                       na.value = "grey50") +
  coord_quickmap() +
  theme_opts +
  guides(fill = guide_colourbar(barwidth = 15, barheight = 1), label.theme = element_text(size=20)) +
  scale_alpha(range = c(-0, 0.5)) +
  scale_x_continuous(expand = c(0, 0), limits = c(100, 170)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-45, 1)) +
  geom_text(data = data_df_long, aes(label = Month), x = 135, y = -25, size = 14, color = "white", hjust = "middle") +
  transition_time(Month2)

animate(gg, nframes = 12, fps = 0.5)
anim_save("Figures/MonthlyBiomass.gif")




