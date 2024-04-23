plot_peak_areas <- function(group, wn_hi = 1800, wn_lo = 950){
  group_spectra <- get_mean_spectra(spectra, metadata$sample, smoothed_derivative = T, w = 7)
  group_spectra <- group_spectra[group_spectra$group == group, ]
  my_peaks <- find_peaks(group_spectra$mean, group_spectra$wavenumber, wn_hi = wn_hi, wn_lo = wn_lo, ndowns = 3, type = "trough")
  plot_spectra_base(group_spectra, wn_hi = wn_hi, wn_lo = wn_lo, col_vector = "#000000", point_size = 0.5)
  #points(my_peaks$peak_wavenumber, my_peaks$abs, col = "red")
  #points(my_peaks$start_wavenumber, my_peaks$start_absorbance, col = "green", pch = 3)
  #points(my_peaks$end_wavenumber, my_peaks$end_absorbance, col = "blue", pch = 4)
  abline(v = c(1605, 1514, 1377, 1169), col = "darkgrey")
  text(my_peaks$peak_wavenumber, my_peaks$nudge, labels = round(my_peaks$peak_wavenumber, 0), cex = 0.5)
  
  my_polygons <- get_peak_areas(group_spectra$mean, group_spectra$wavenumber, peak_data = my_peaks)
  plot_polygons <- lapply(my_polygons, "[[", 2)
  polygon_areas <- sapply(my_polygons, "[[", 1)
  
  plot_peak_polygons(plot_polygons)
}