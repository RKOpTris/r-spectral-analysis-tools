# loading libraries and functions
library(dplyr)
source("get_mean_spectra.R")
source("strip_non_numbers.R")
source("find_peaks.R")
source("nearest.R")
source("plot_spectra_base.R")
source("inRange.R")
source("dec_to_hex.R")
source("internal_axis_ticks.R")
source("get_peak_areas.R")
source("plot_peak_polygons.R")

# loading data (these data are already scatter-corrected and scaled)
metadata <- readRDS("Pinaceae_Cupressaceae_metadata.RDS")
spectra <- readRDS("Pinaceae_Cupressaceae_spectra.RDS")
metadata <- metadata$metadata
spectra <- spectra$spectra

# subset data to one sample
tidy_spectra_sub <- get_mean_spectra(spectra, metadata$sample, smoothed_derivative = T, w = 7) %>% 
  filter(group == "PINUS-COULTERI")

# set wavenumber range
wn_hi <- 1800
wn_lo <- 950

# detect peaks
my_peaks <- find_peaks(tidy_spectra_sub$mean, tidy_spectra_sub$wavenumber, wn_hi = wn_hi, wn_lo = wn_lo, ndowns = 3, type = "trough")

# plot the sample mean with the peak numbers
plot_spectra_base(tidy_spectra_sub, wn_hi = wn_hi, wn_lo = wn_lo, col_vector = "#000000", point_size = 0.5)
points(my_peaks$peak_wavenumber, my_peaks$abs, col = "red")
points(my_peaks$start_wavenumber, my_peaks$start_absorbance, col = "green", pch = 3)
points(my_peaks$end_wavenumber, my_peaks$end_absorbance, col = "blue", pch = 4)
text(my_peaks$peak_wavenumber, my_peaks$nudge, labels = round(my_peaks$peak_wavenumber, 0), cex = 0.5)

# find the peak polygons
my_polygons <- get_peak_areas(tidy_spectra_sub$mean, tidy_spectra_sub$wavenumber, my_peaks)
plot_polygons <- lapply(my_polygons, "[[", 2)
polygon_areas <- sapply(my_polygons, "[[", 1)

plot_polygons[[1]] # the first polygon
polygon_areas[[1]] # the area of the first polygon

# visualise the peak area polygons in a new plot
plot_spectra_base(tidy_spectra_sub, wn_hi = wn_hi, wn_lo = wn_lo, col_vector = "#000000", point_size = 0.5)
text(my_peaks$peak_wavenumber, my_peaks$nudge, labels = round(my_peaks$peak_wavenumber, 0), cex = 0.5)
plot_peak_polygons(plot_polygons)
