PlotGLMM <- function(model, FileName) {
  # Note: in visreg(), plot = T, gg = T, turns it into a ggplot object,
  # but plot = F, gg = F turns the object into a visreg object :-(.
  # Problem is turning it into a visreg object if you don't have points on gets y-axis scale wrong. Easier to build ggplot from scratch.
  # Put in a loop maybe sometime the different layers in Method, but could be tricky. # for(i in 1:NlevelBiomMeth){
  # Ymin <- c(-5, -5, -5, -5, -5, -5, -7)
  # Ymax <- c(3, 10, 5, 5, 5, 3, 3)
  Ymin <- c(-5, -5, -5, -5, -5, -8)
  Ymax <- c(3, 10, 4, 4, 4, 4)

  # 1. YEAR
  Tempory <- visreg(model, xvar = "Year", gg = TRUE)
  p1 <- ggplot() +
    geom_point(
      data = Tempory$data,
      aes(x = x, y = y),
      alpha = 0.01,
      stroke = NA
    ) +
    geom_polygon(
      data = Tempory$layers[[1]]$data,
      aes(x = x, y = y),
      colour = NA,
      fill = "slateblue1",
      alpha = 0.5
    ) +
    geom_line(
      data = Tempory$layers[[3]]$data,
      aes(x = x, y = y),
      colour = "blue",
      linewidth = 1
    ) +
    xlab(Tempory$labels$x) +
    ylab("log(Biomass)") +
    ylim(Ymin[1], Ymax[1])

  # 2. METHOD
  Tempory <- visreg(model, xvar = "Method", gg = TRUE)
  Labels <- ggplot_build(Tempory)$layout$panel_params[[1]]$x$get_labels() # Method labels
  Posns <- ggplot_build(Tempory)$layout$panel_params[[1]]$x$break_positions() # Method label positions

  p2 <- ggplot() +

    geom_point(
      data = Tempory$layers[[7]]$data,
      aes(x = visregGGX, y = visregGGY),
      alpha = 0.01,
      stroke = NA
    ) +

    geom_polygon(
      data = Tempory$layers[[1]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = NA,
      fill = "slateblue1",
      alpha = 0.5
    ) +
    geom_polygon(
      data = Tempory$layers[[2]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = NA,
      fill = "slateblue1",
      alpha = 0.5
    ) +
    geom_polygon(
      data = Tempory$layers[[3]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = NA,
      fill = "slateblue1",
      alpha = 0.5
    ) +
    geom_polygon(
      data = Tempory$layers[[4]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = NA,
      fill = "slateblue1",
      alpha = 0.5
    ) +
    geom_polygon(
      data = Tempory$layers[[5]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = NA,
      fill = "slateblue1",
      alpha = 0.5
    ) +
    geom_polygon(
      data = Tempory$layers[[6]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = NA,
      fill = "slateblue1",
      alpha = 0.5
    ) +

    geom_line(
      data = Tempory$layers[[8]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = "blue",
      linewidth = 1
    ) +
    geom_line(
      data = Tempory$layers[[9]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = "blue",
      linewidth = 1
    ) +
    geom_line(
      data = Tempory$layers[[10]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = "blue",
      linewidth = 1
    ) +
    geom_line(
      data = Tempory$layers[[11]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = "blue",
      linewidth = 1
    ) +
    geom_line(
      data = Tempory$layers[[12]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = "blue",
      linewidth = 1
    ) +
    geom_line(
      data = Tempory$layers[[13]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = "blue",
      linewidth = 1
    ) +

    # Remove ticks and numbers
    scale_x_continuous(breaks = Posns, labels = Labels) +
    theme(axis.text.x = element_text(angle = 90)) +
    xlab("Method") +
    ylab("log(Biomass)") +
    ylim(Ymin[2], Ymax[2])

  # 3. MESH
  Tempory <- visreg(model, xvar = "Mesh", gg = TRUE)
  p3 <- ggplot() +
    geom_point(
      data = Tempory$data,
      aes(x = x, y = y),
      alpha = 0.01,
      stroke = NA
    ) +
    geom_polygon(
      data = Tempory$layers[[1]]$data,
      aes(x = x, y = y),
      colour = NA,
      fill = "slateblue1",
      alpha = 0.5
    ) +
    geom_line(
      data = Tempory$layers[[3]]$data,
      aes(x = x, y = y),
      colour = "blue",
      linewidth = 1
    ) +
    xlab(Tempory$labels$x) +
    ylab("log(Biomass)") +
    ylim(Ymin[3], Ymax[3])

  # 4. TOW
  Tempory <- visreg(model, xvar = "Tow", gg = TRUE)
  Labels <- ggplot_build(Tempory)$layout$panel_params[[1]]$x$get_labels() # Tow labels
  Posns <- ggplot_build(Tempory)$layout$panel_params[[1]]$x$break_positions() # Tow label positions

  p4 <- ggplot() +

    geom_point(
      data = Tempory$layers[[4]]$data,
      aes(x = visregGGX, y = visregGGY),
      alpha = 0.01,
      stroke = NA
    ) +

    geom_polygon(
      data = Tempory$layers[[1]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = NA,
      fill = "slateblue1",
      alpha = 0.5
    ) +
    geom_polygon(
      data = Tempory$layers[[2]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = NA,
      fill = "slateblue1",
      alpha = 0.5
    ) +
    geom_polygon(
      data = Tempory$layers[[3]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = NA,
      fill = "slateblue1",
      alpha = 0.5
    ) +

    geom_line(
      data = Tempory$layers[[5]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = "blue",
      linewidth = 1
    ) +
    geom_line(
      data = Tempory$layers[[6]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = "blue",
      linewidth = 1
    ) +
    geom_line(
      data = Tempory$layers[[7]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = "blue",
      linewidth = 1
    ) +

    # Remove ticks and numbers
    scale_x_continuous(breaks = Posns, labels = Labels) +
    xlab(Tempory$labels$x) +
    ylab("log(Biomass)") +
    ylim(Ymin[4], Ymax[4])

  # 5. DEPTH
  Tempory <- visreg(model, xvar = "Depth", gg = TRUE)
  p5 <- ggplot() +
    geom_point(
      data = Tempory$data,
      aes(x = x, y = y),
      alpha = 0.01,
      stroke = NA
    ) +
    geom_polygon(
      data = Tempory$layers[[1]]$data,
      aes(x = x, y = y),
      colour = NA,
      fill = "slateblue1",
      alpha = 0.5
    ) +
    geom_line(
      data = Tempory$layers[[3]]$data,
      aes(x = x, y = y),
      colour = "blue",
      linewidth = 1
    ) +
    xlab(Tempory$labels$x) +
    ylab("log(Biomass)") +
    ylim(Ymin[5], Ymax[5])

  # Tempory <- visreg(model, xvar = "HarmTOD", gg = TRUE)
  # p6 <- ggplot() +
  #   geom_point(data = Tempory$data, aes(x = x, y = y), alpha = 0.01, stroke = NA) +
  #   geom_polygon(data = Tempory$layers[[1]]$data, aes(x = x, y = y), colour = NA, fill = "slateblue1", alpha = 0.5) +
  #   geom_line(data = Tempory$layers[[3]]$data, aes(x = x, y = y), colour = "blue", linewidth = 1) +
  #   xlab(Tempory$labels$x) +
  #   ylab("log(Biomass)") +
  #   ylim(Ymin[6], Ymax[6])

  # 6. HarmDOY BY LATITUDE
  Tempory <- visreg(
    model,
    xvar = "HarmDOY",
    by = "Latitude",
    breaks = c(-60, -30, 0, 30, 60),
    overlay = TRUE,
    gg = TRUE
  )
  p6 <- ggplot() +
    geom_point(
      data = Tempory$data,
      aes(x = x, y = y, colour = Latitude),
      alpha = 0.01,
      stroke = NA
    ) +
    geom_polygon(
      data = Tempory$layers[[1]]$data,
      aes(x = x, y = y, fill = Latitude),
      colour = NA,
      alpha = 0.5
    ) +
    geom_line(
      data = Tempory$layers[[3]]$data,
      aes(x = x, y = y, colour = Latitude),
      linewidth = 1
    ) +
    xlab("Day of year") +
    ylab("log(Biomass)") +
    ylim(Ymin[6], Ymax[6])

  Plots <- p1 + p2 + p3 + p4 + p5 + p6
  ggsave(
    paste0("Figures/", FileName, ".png"),
    Plots,
    width = 12,
    height = 8,
    dpi = 300
  )
  return(Plots)
}


PlotGLMM_Map <- function(model, FileName) {
  # Biomass ~ Method + Mesh + Tow + log(Chl_gm3) + ns(Thetao, df = 3) +
  # ns(Bathy, df = 3) +
  # harmonic(HarmTOD, k = 1)*ns(Depth, df = 3) +
  # ns(Latitude, 2) * harmonic(HarmDOY, k = 2) +
  # (1|DatasetID) + (1|Gear),

  Ymin <- c(-5, -5, -5, -5, -5, -6, -5, -5)
  Ymax <- c(10, 3, 3, 3, 2, 3, 3, 0)
  # p5, p6, p7, p8 done for Ylims

  # 1. METHOD
  Tempory <- visreg(model, xvar = "Method", gg = TRUE)
  Labels <- ggplot_build(Tempory)$layout$panel_params[[1]]$x$get_labels() # Method labels
  Posns <- ggplot_build(Tempory)$layout$panel_params[[1]]$x$break_positions() # Method label positions

  p1 <- ggplot() +

    geom_point(
      data = Tempory$layers[[7]]$data,
      aes(x = visregGGX, y = visregGGY),
      alpha = 0.01,
      stroke = NA
    ) +

    geom_polygon(
      data = Tempory$layers[[1]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = NA,
      fill = "slateblue1",
      alpha = 0.5
    ) +
    geom_polygon(
      data = Tempory$layers[[2]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = NA,
      fill = "slateblue1",
      alpha = 0.5
    ) +
    geom_polygon(
      data = Tempory$layers[[3]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = NA,
      fill = "slateblue1",
      alpha = 0.5
    ) +
    geom_polygon(
      data = Tempory$layers[[4]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = NA,
      fill = "slateblue1",
      alpha = 0.5
    ) +
    geom_polygon(
      data = Tempory$layers[[5]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = NA,
      fill = "slateblue1",
      alpha = 0.5
    ) +
    geom_polygon(
      data = Tempory$layers[[6]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = NA,
      fill = "slateblue1",
      alpha = 0.5
    ) +

    geom_line(
      data = Tempory$layers[[8]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = "blue",
      linewidth = 1
    ) +
    geom_line(
      data = Tempory$layers[[9]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = "blue",
      linewidth = 1
    ) +
    geom_line(
      data = Tempory$layers[[10]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = "blue",
      linewidth = 1
    ) +
    geom_line(
      data = Tempory$layers[[11]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = "blue",
      linewidth = 1
    ) +
    geom_line(
      data = Tempory$layers[[12]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = "blue",
      linewidth = 1
    ) +
    geom_line(
      data = Tempory$layers[[13]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = "blue",
      linewidth = 1
    ) +

    # Remove ticks and numbers
    scale_x_continuous(breaks = Posns, labels = Labels) +
    theme(
      axis.text.x = element_text(angle = 45, margin = margin(t = 20, l = -50)),
      axis.title.x = element_text(margin = margin(t = -20))
    ) +
    xlab("Method") +
    ylab("log(Biomass)") +
    ylim(Ymin[1], Ymax[1])

  # 2. MESH
  Tempory <- visreg(model, xvar = "Mesh", gg = TRUE)
  p2 <- ggplot() +
    geom_point(
      data = Tempory$data,
      aes(x = x, y = y),
      alpha = 0.01,
      stroke = NA
    ) +
    geom_polygon(
      data = Tempory$layers[[1]]$data,
      aes(x = x, y = y),
      colour = NA,
      fill = "slateblue1",
      alpha = 0.5
    ) +
    geom_line(
      data = Tempory$layers[[3]]$data,
      aes(x = x, y = y),
      colour = "blue",
      linewidth = 1
    ) +
    xlab(Tempory$labels$x) +
    ylab("log(Biomass)") +
    ylim(Ymin[2], Ymax[2])

  # 3. TOW
  Tempory <- visreg(model, xvar = "Tow", gg = TRUE)
  Labels <- ggplot_build(Tempory)$layout$panel_params[[1]]$x$get_labels() # Tow labels
  Posns <- ggplot_build(Tempory)$layout$panel_params[[1]]$x$break_positions() # Tow label positions

  p3 <- ggplot() +

    geom_point(
      data = Tempory$layers[[4]]$data,
      aes(x = visregGGX, y = visregGGY),
      alpha = 0.01,
      stroke = NA
    ) +

    geom_polygon(
      data = Tempory$layers[[1]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = NA,
      fill = "slateblue1",
      alpha = 0.5
    ) +
    geom_polygon(
      data = Tempory$layers[[2]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = NA,
      fill = "slateblue1",
      alpha = 0.5
    ) +
    geom_polygon(
      data = Tempory$layers[[3]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = NA,
      fill = "slateblue1",
      alpha = 0.5
    ) +

    geom_line(
      data = Tempory$layers[[5]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = "blue",
      linewidth = 1
    ) +
    geom_line(
      data = Tempory$layers[[6]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = "blue",
      linewidth = 1
    ) +
    geom_line(
      data = Tempory$layers[[7]]$data,
      aes(x = visregGGX, y = visregGGY),
      colour = "blue",
      linewidth = 1
    ) +

    # Remove ticks and numbers
    scale_x_continuous(breaks = Posns, labels = Labels) +
    xlab(Tempory$labels$x) +
    ylab("log(Biomass)") +
    ylim(Ymin[3], Ymax[3])

  # 4. CHL_GM3
  Tempory <- visreg(model, xvar = "Chl_gm3", gg = TRUE)
  p4 <- ggplot() +
    geom_point(
      data = Tempory$data,
      aes(x = x, y = y),
      alpha = 0.01,
      stroke = NA
    ) +
    geom_polygon(
      data = Tempory$layers[[1]]$data,
      aes(x = x, y = y),
      colour = NA,
      fill = "slateblue1",
      alpha = 0.5
    ) +
    geom_line(
      data = Tempory$layers[[3]]$data,
      aes(x = x, y = y),
      colour = "blue",
      linewidth = 1
    ) +
    xlab(Tempory$labels$x) +
    ylab("log(Biomass)") +
    ylim(Ymin[4], Ymax[4])

  # 5. THETAO
  Tempory <- visreg(model, xvar = "Thetao", gg = TRUE)
  p5 <- ggplot() +
    geom_point(
      data = Tempory$data,
      aes(x = x, y = y),
      alpha = 0.01,
      stroke = NA
    ) +
    geom_polygon(
      data = Tempory$layers[[1]]$data,
      aes(x = x, y = y),
      colour = NA,
      fill = "slateblue1",
      alpha = 0.5
    ) +
    geom_line(
      data = Tempory$layers[[3]]$data,
      aes(x = x, y = y),
      colour = "blue",
      linewidth = 1
    ) +
    xlab(Tempory$labels$x) +
    ylab("log(Biomass)") +
    ylim(Ymin[5], Ymax[5])

  # 6. BATHY
  Tempory <- visreg(model, xvar = "Bathy", gg = TRUE)
  p6 <- ggplot() +
    geom_point(
      data = Tempory$data,
      aes(x = x, y = y),
      alpha = 0.01,
      stroke = NA
    ) +
    geom_polygon(
      data = Tempory$layers[[1]]$data,
      aes(x = x, y = y),
      colour = NA,
      fill = "slateblue1",
      alpha = 0.5
    ) +
    geom_line(
      data = Tempory$layers[[3]]$data,
      aes(x = x, y = y),
      colour = "blue",
      linewidth = 1
    ) +
    xlab("Bathymetry") +
    ylab("log(Biomass)") +
    ylim(Ymin[6], Ymax[6])

  # 7. HarmTOD BY DEPTH
  Tempory <- visreg(
    model,
    xvar = "HarmTOD",
    by = "Depth",
    breaks = c(0, 100, 300, 1000),
    overlay = TRUE,
    gg = TRUE
  )
  p7 <- ggplot() +
    geom_point(
      data = Tempory$data,
      aes(x = x, y = y, colour = Depth),
      alpha = 0.01,
      stroke = NA
    ) +
    geom_polygon(
      data = Tempory$layers[[1]]$data,
      aes(x = x, y = y, fill = Depth),
      colour = NA,
      alpha = 0.5
    ) +
    geom_line(
      data = Tempory$layers[[3]]$data,
      aes(x = x, y = y, colour = Depth),
      linewidth = 1
    ) +
    xlab("Time of day") +
    ylab("log(Biomass)") +
    ylim(Ymin[7], Ymax[7]) +
    theme(legend.position = "bottom")

  # 8. HarmDOY BY LATITUDE
  Tempory <- visreg(
    model,
    xvar = "HarmDOY",
    by = "Latitude",
    breaks = c(-60, -30, 0, 30, 60),
    overlay = TRUE,
    gg = TRUE
  )
  p8 <- ggplot() +
    geom_point(
      data = Tempory$data,
      aes(x = x, y = y, colour = Latitude),
      alpha = 0.01,
      stroke = NA
    ) +
    geom_polygon(
      data = Tempory$layers[[1]]$data,
      aes(x = x, y = y, fill = Latitude),
      colour = NA,
      alpha = 0.5
    ) +
    geom_line(
      data = Tempory$layers[[3]]$data,
      aes(x = x, y = y, colour = Latitude),
      linewidth = 1
    ) +
    xlab("Day of year") +
    ylab("log(Biomass)") +
    ylim(Ymin[8], Ymax[8]) +
    theme(legend.position = "bottom")

  Plots <- p1 +
    p2 +
    patchwork::plot_spacer() +
    patchwork::plot_spacer() +
    p3 +
    p4 +
    patchwork::plot_spacer() +
    patchwork::plot_spacer() +
    p5 +
    p6 +
    patchwork::plot_spacer() +
    patchwork::plot_spacer() +
    p7 +
    p8 +
    patchwork::plot_spacer() +
    patchwork::plot_spacer() +
    patchwork::guide_area()
  patchwork::plot_layout(ncol = 4, guides = "collect")

  ggsave(
    paste0("Figures/", FileName, ".png"),
    Plots,
    width = 12,
    height = 8,
    dpi = 600
  )

  out <- list(p1, p2, p3, p4, p5, p6, p7, p8)

  return(out)
}
