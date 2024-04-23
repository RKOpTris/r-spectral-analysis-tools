plot_peak_polygons <- function(plot_polygons, col = "#FFFF0077"){
  for(i in 1:length(plot_polygons)){
    poly_dat <- plot_polygons[[i]]
    polygon(poly_dat$x, poly_dat$y, col = col)
  }
}