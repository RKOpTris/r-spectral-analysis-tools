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