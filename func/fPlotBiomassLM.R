fPlotBiomassLM <- function (m1, Name) {
  Terms <- as.character(m1$terms)[3] # Terms from the model so we can print blank if n.s.
  
  x11(width = 15, height = 6)
  par(mfrow = c(2,5), mar = c(4,4,2,2))
  
  # if(grepl("BiomassMethod", Terms, fixed = TRUE)) {
  #   visreg(m1, "BiomassMethod", rug = FALSE, scale = "response", xlab = "Method")} else {
  #     plot.new()
  #   }
  
  # if(grepl("Tow", Terms, fixed = TRUE)) {
  #   visreg(m1, "Tow", rug = FALSE, scale = "response", xlab = "Tow")} else {
  #     plot.new()
  #   }
  # 
  # if(grepl("Mesh", Terms, fixed = TRUE)) {
  #   visreg(m1, "Mesh", scale = "response", xlab = "Mesh  (microns)")} else {
  #     plot.new()
  #   }
  
  if(grepl("Chl", Terms, fixed = TRUE)) {
    visreg(m1, "Chl", scale = "response", xlab = "Chl-a (mg/m3)")} else {
      plot.new()
    }
  if(grepl("Bathy", Terms, fixed = TRUE)) {
    visreg(m1, "Bathy", scale = "response", xlab = "Bathy (m)")} else {
      plot.new()
    }
  
  if(grepl("HarmTOD", Terms, fixed = TRUE)) {
    visreg(m1, "HarmTOD", rug = FALSE, scale = "response", xlab = "Time of day")} else {
      plot.new()
    }
  
  # if(grepl("Depth", Terms, fixed = TRUE)) {
  #   visreg(m1, "Depth", scale = "response", xlab = "Depth")} else {
  #     plot.new()
  #   }
  # 
  # if(grepl("exp(-Depth/1000) * fHarmonic(HarmTOD, k = 1)", Terms, fixed = TRUE) | grepl("fHarmonic(HarmTOD, k = 1) * exp(-Depth/1000)", Terms, fixed = TRUE)) {
  #   visreg(m1, "HarmTOD", by = "Depth", breaks = c(0, 50, 100, 300, 500), xlab = "Time of day", ylab = "log10(Abundance)", 
  #          type = "contrast", scale = "response", overlay = TRUE, rug = 0, strip.names = c("Depth=0","Depth=50", "Depth=100", "Depth=300", "Depth=500"))} else {
  #            plot.new()
  #          }
  
  if(grepl("HarmDOY", Terms, fixed = TRUE)) {
    visreg(m1, "HarmDOY", scale = "response", xlab = "Day of year")} else {
      plot.new()
    }
  
  if(grepl("SST", Terms, fixed = TRUE)) {
    visreg(m1, "SST", scale = "response", xlab = "SST (ºC)")} else {
      plot.new()
    }
  
  if(grepl("fHarmonic(HarmDOY, k = 1) * ns(SST, 3)", Terms, fixed = TRUE)) {
    visreg2d(m1, yvar = "HarmDOY", xvar = "SST", scale = "response", plot.type = "persp", theta = 45, phi = 10, r = 100, 
             ticktype = "detailed", xlab = "\nSST (ºC)", ylab = "\nDay of Year", zlab = "\nlog10(Biomass)", color = "deepskyblue2")} else {
               plot.new()
             }
  
  dev.print(pdf, paste0("Figures/", Name, ".pdf"))
}