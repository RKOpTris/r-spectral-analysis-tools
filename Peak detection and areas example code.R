metadata <- readRDS("Pinaceae_Cupressaceae_metadata.RDS")
spectra <- readRDS("Pinaceae_Cupressaceae_spectra.RDS")
metadata <- metadata$metadata
spectra <- spectra$spectra


tidy_spectra_sub <- get_mean_spectra(spectra, metadata$sample, smoothed_derivative = T, w = 7) %>% 
  filter(group == "PINUS-COULTERI")
wn_hi <- 1800
wn_lo <- 950

my_peaks <- find_peaks(tidy_spectra_sub$mean, tidy_spectra_sub$wavenumber, wn_hi = wn_hi, wn_lo = wn_lo, ndowns = 3, type = "trough")
plot_spectra_base(tidy_spectra_sub, wn_hi = wn_hi, wn_lo = wn_lo, col_vector = "#000000", point_size = 0.5)
points(my_peaks$peak_wavenumber, my_peaks$abs, col = "red")
points(my_peaks$start_wavenumber, my_peaks$start_absorbance, col = "green", pch = 3)
points(my_peaks$end_wavenumber, my_peaks$end_absorbance, col = "blue", pch = 4)
abline(v = c(1605, 1514, 1377, 1169), col = "darkgrey")
text(my_peaks$peak_wavenumber, my_peaks$nudge, labels = round(my_peaks$peak_wavenumber, 0), cex = 0.5)

get_peak_areas <- function(x, wns, peak_data, trancated_areas = F){
  lapply(1:nrow(peak_data), function(n){
    section <- peak_data %>% select(peak_start, peak_end, peak_ind) %>% slice(n) %>% as.numeric()
    peak_dat <- x[section[1]:section[2]]
    return_dat <- seq(x[section[2]], x[section[1]], length.out = length(section[1]:section[2]))
    #### this will only work for derivative rather than raw spectrum because negative values
    #### could possibly fix by converting to absolute values and then seeing if return exceeds orig values
    #### or by refering to my_peaks to check for trough = derivative or peak = raw spectrum
    return_dat <- rev(return_dat)
    check_outside_boundary <- return_dat < peak_dat
    return_dat[check_outside_boundary] <- peak_dat[check_outside_boundary]
    return_dat <- rev(return_dat)
    peak_wns <- wns[section[1]:section[2]]
    return_wns <- rev(wns[section[1]:section[2]])
    poly_dat <- data.frame(x = c(peak_wns, return_wns), y = c(peak_dat, return_dat))
    area <- area::polygon_area(as.matrix(poly_dat))
    list(polygon_area = area, polygon_data = poly_dat)
  })
}

my_polygons <- get_peak_areas(tidy_spectra_sub$mean, tidy_spectra_sub$wavenumber, my_peaks)
plot_polygons <- lapply(my_polygons, "[[", 2)
polygon_areas <- sapply(my_polygons, "[[", 1)

plot_peak_polygons <- function(plot_polygons, col = "#FFFF0077"){
  for(i in 1:length(plot_polygons)){
    poly_dat <- plot_polygons[[i]]
    polygon(poly_dat$x, poly_dat$y, col = col)
  }
}

plot_peak_polygons(plot_polygons)
