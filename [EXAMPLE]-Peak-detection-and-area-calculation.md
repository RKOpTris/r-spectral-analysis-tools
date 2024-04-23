\[EXAMPLE\] Peak detection and area calculation of a derivative spectrum
================
RKOpTris
2024-04-23

## Loading libraries

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
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
source("plot_peak_areas.R")
```

## Loading data (these data are already scatter-corrected and scaled)

``` r
metadata <- readRDS("Pinaceae_Cupressaceae_metadata.RDS")
spectra <- readRDS("Pinaceae_Cupressaceae_spectra.RDS")
metadata <- metadata$metadata
spectra <- spectra$spectra
```

## Subset data to one sample

``` r
tidy_spectra_sub <- get_mean_spectra(spectra, metadata$sample, smoothed_derivative = T, w = 7) %>% 
  filter(group == "PINUS-COULTERI")
```

    ## prospectr::savitzkyGolay is using w = 7, p = 2, m = 2

    ## `summarise()` has grouped output by 'group'. You can override using the
    ## `.groups` argument.

## Set wavenumber range

``` r
wn_hi <- 1800
wn_lo <- 950
```

## Detect peaks

``` r
my_peaks <- find_peaks(tidy_spectra_sub$mean, tidy_spectra_sub$wavenumber, wn_hi = wn_hi, wn_lo = wn_lo, ndowns = 3, type = "trough")
```

## Plot the sample mean with the peak numbers

``` r
plot_spectra_base(tidy_spectra_sub, wn_hi = wn_hi, wn_lo = wn_lo, col_vector = "#000000", point_size = 0.5)
points(my_peaks$peak_wavenumber, my_peaks$abs, col = "red")
points(my_peaks$start_wavenumber, my_peaks$start_absorbance, col = "green", pch = 3)
points(my_peaks$end_wavenumber, my_peaks$end_absorbance, col = "blue", pch = 4)
text(my_peaks$peak_wavenumber, my_peaks$nudge, labels = round(my_peaks$peak_wavenumber, 0), cex = 0.5)
```

![](%5BEXAMPLE%5D-Peak-detection-and-area-calculation_files/figure-gfm/s05-1.png)<!-- -->

## Find the peak polygons

``` r
my_polygons <- get_peak_areas(tidy_spectra_sub$mean, tidy_spectra_sub$wavenumber, my_peaks)
plot_polygons <- lapply(my_polygons, "[[", 2)
polygon_areas <- sapply(my_polygons, "[[", 1)

plot_polygons[[1]] # the first polygon
```

    ##         x             y
    ## 1   970.0  1.395011e-03
    ## 2   971.9  1.186164e-03
    ## 3   973.9  5.682777e-04
    ## 4   975.8 -3.977810e-04
    ## 5   977.7 -1.466349e-03
    ## 6   979.7 -2.296607e-03
    ## 7   981.6 -2.628448e-03
    ## 8   983.5 -2.420816e-03
    ## 9   985.4 -1.850866e-03
    ## 10  987.4 -1.195732e-03
    ## 11  989.3 -6.793150e-04
    ## 12  991.2 -3.701893e-04
    ## 13  993.2 -2.052236e-04
    ## 14  995.1 -1.075811e-04
    ## 15  997.0 -3.209418e-05
    ## 16  998.9  7.341836e-05
    ## 17 1000.9  2.292168e-04
    ## 18 1002.8  3.391409e-04
    ## 19 1002.8  3.391409e-04
    ## 20 1000.9  4.012509e-04
    ## 21  998.9  4.633609e-04
    ## 22  997.0  5.254709e-04
    ## 23  995.1  5.875810e-04
    ## 24  993.2  6.496910e-04
    ## 25  991.2  7.118010e-04
    ## 26  989.3  7.739110e-04
    ## 27  987.4  8.360211e-04
    ## 28  985.4  8.981311e-04
    ## 29  983.5  9.602411e-04
    ## 30  981.6  1.022351e-03
    ## 31  979.7  1.084461e-03
    ## 32  977.7  1.146571e-03
    ## 33  975.8  1.208681e-03
    ## 34  973.9  1.270791e-03
    ## 35  971.9  1.332901e-03
    ## 36  970.0  1.395011e-03

``` r
polygon_areas[[1]] # the area of the first polygon
```

    ## [1] 0.04909371

## Visualise the peak area polygons in a new plot

``` r
plot_spectra_base(tidy_spectra_sub, wn_hi = wn_hi, wn_lo = wn_lo, col_vector = "#000000", point_size = 0.5)
text(my_peaks$peak_wavenumber, my_peaks$nudge, labels = round(my_peaks$peak_wavenumber, 0), cex = 0.5)
plot_peak_polygons(plot_polygons)
```

![](%5BEXAMPLE%5D-Peak-detection-and-area-calculation_files/figure-gfm/s07-1.png)<!-- -->

## Do this all with one line for another sample

``` r
plot_peak_areas("CRYPTOMERIA-JAPONICA")
```

    ## prospectr::savitzkyGolay is using w = 7, p = 2, m = 2

    ## `summarise()` has grouped output by 'group'. You can override using the
    ## `.groups` argument.

![](%5BEXAMPLE%5D-Peak-detection-and-area-calculation_files/figure-gfm/s08-1.png)<!-- -->
