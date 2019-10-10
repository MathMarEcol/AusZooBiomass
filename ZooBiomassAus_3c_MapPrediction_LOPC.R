# https://stats.stackexchange.com/questions/131106/predicting-with-random-effects-in-mgcv-gam

library(tidyverse)
library(raster)
library(rasterImage)
library(sf)
library(ggpubr)
library(splines)

source("../ZoopAbundanceModels/fHarmonic.R")
source("func/fPlotAusLMER.R")
source("../ZoopAbundanceModels/fBuildRasterBrick.R")
source("../ZoopAbundanceModels/fWorldBorders.R")

# Set the max values to be removed in rasterbrick
maxm <- tibble(Bathy = 6000, SST = 28, Chl = 10)

####################################################
NthMonth <- c("01January", "02February","03March","04April","05May","06June",
              "07July","08August","09September","10October","11November","12December")

df <- data.frame(HarmTOD = 1, LOPC_Type = "LOPC", Depth = 10)

# These are zeroed out in the final model, but they need to be included (in df) initially
# excl <- c("s(Project)") 
##############################################################################

# Predict abundance using the predict function in raster and plot
tframe <- "annual"

quick <- 1
CRS <- CRS("+proj=longlat +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84")

myplots <- list() # Initialise plot list

# Define themes for plots 
theme_opts <- list(theme(panel.grid.minor = element_blank(),
                         panel.grid.major = element_blank(),
                         panel.background = element_blank(),
                         plot.background = element_rect(fill="white"),
                         plot.title = element_text(hjust = 0.5, size = rel(0.5)),
                         panel.border = element_blank(),
                         axis.line = element_blank(),
                         axis.text.x = element_blank(),
                         axis.text.y = element_blank(),
                         axis.ticks = element_blank(),
                         axis.title.x = element_blank(),
                         axis.title.y = element_blank(),
                         legend.position = "bottom",
                         legend.text = element_text(size = rel(0.5))))

mdl <- readRDS("Biomass_LOPC_Final.rds")

for (i in 1:length(NthMonth)) {
  
  # The data brick of satellite data, and DOY
  data_brick <- fBuildRasterBrick(NthMonth[i], maxm,"one")
  
  layer <- predict(data_brick, mdl, 
                   const = df, re.form=NULL)
  
  names(layer) <- paste0(gsub("\\d", "", NthMonth[i]))
  
  if (i == 1){
    stat_brick <- brick(layer)} else{
      stat_brick <- addLayer(stat_brick,layer)
    }
  
}

saveRDS(stat_brick, file = "Biomass_Layer_LOPC.rds")

if (tframe == "annual") {
  stat_brick  <- calc(stat_brick, fun = mean)
  n <- 1
}
if (tframe=="month"){
  n <- 12
}

data <- projectRaster(stat_brick, crs = CRS, over = TRUE)

# Load the World Borders in the Moll projection
wb_df <- fWorldBorders()


for (i in 1:n){
  # Convert data to dataframe
  data_df <- as.data.frame(rasterToPoints(data))
  
  ## Colormap options
  col <- rev(rasterImage::colorPalette(n = 14, type = "jet.colors")) # using rasterImage pckge
  
  # ## Do the plotting
  gg <- ggplot() +
    geom_tile(data = data_df, aes_string(x = "x", y = "y", fill = names(data_df[i+2]))) +
    geom_polygon(fill = "white") +
    geom_polygon(data = wb_df, aes(x = long, y = lat, group = group)) +
    scale_fill_gradientn(name =  expression(paste("mg m"^-3)), colours = rev(col),
                         #limits = clim, breaks = cticks,
                         #labels = clabel, 
                         position = "bottom",
                         na.value = "grey50") +
    {if (quick==1) coord_quickmap() else coord_map(projection = "mollweide") } +
    theme_opts +
    scale_alpha(range = c(-0, 0.5)) +
    scale_x_continuous(expand = c(0, 0), limits = c(100, 170)) +
    scale_y_continuous(expand = c(0, 0), limits = c(-60, 1))
  
  if (tframe == "month"){
    gg <- gg +  ggtitle(paste0(gsub("\\d", "", NthMonth[i])," - ",taxa)) + 
      guides(fill = guide_colorbar(barheight = 0.5, barwidth = 10))
  }
  
  if (tframe == "annual"){
    gg <- gg +  ggtitle("LOPC") + 
      guides(fill = guide_colorbar(barheight = 2, barwidth = 20)) + 
      theme(legend.text = element_text(size = rel(2)),
            legend.title = element_text(size = rel(2)),
            plot.title = element_text(size = rel(2)))
  }
  
  myplots[[i]] <- gg
}


if (tframe=="month"){
  graphics.off()
  x11(width = 15, height = 15)
  ggarrange(plotlist=myplots, ncol = 3, nrow = 4)
  ggsave(paste0("Figures/GlobalMaps/GlobalModel_All_",taxa,".png"), dpi=300, scale=1,
         width = 210, height = 297, units = "mm")
} 

if (tframe=="annual"){
  graphics.off()
  x11(width = 15, height = 15)
  ggsave(paste0("Figures/AusBiomass_Annual_LOPC.png"), dpi=300, scale=1,
         width = 210, height = 297, units = "mm")
} 

