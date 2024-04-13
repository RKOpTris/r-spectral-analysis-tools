### calculate the mean spectra by a grouping vector, giving a tidy output
### assumes sample column to average with other columns being wavenumbers
### smoothed_derivative (logical), and will return smoothed derivatives of the spectra according to w, p and m
### sd_upper_lower (logical) returns mean + sd and mean - sd for plotting sd ribbons
### sd (logical) returns the sd which may not be needed if upper and lower are already calculated
### w (numeric) smoothing window size
### p (numeric) smoothing polynomial
### m (numeric) degree of differentiation

get_mean_spectra <- function(spectra, grouping_vector, smoothed_derivative = F, sd_upper_lower = T, sd = F, w = 3, p = 2, m = 2, ...){
  wavenumber_cols <- strip_non_numbers(names(spectra)) %>% sapply(is.numeric)
  spectra <- spectra[wavenumber_cols]
  if(smoothed_derivative){
    message(paste0("prospectr::savitzkyGolay is using w = ", w, ", p = ", p, ", m = ", m))
    spectra <- spectra %>% prospectr::savitzkyGolay(w = w, p = p, m = m, ...) %>% data.frame()
  }
  spectra <- spectra %>% dplyr::mutate(group = grouping_vector) %>%
    tidyr::gather("wavenumber", "absorbance", -group) %>%
    dplyr::mutate(wavenumber = strip_non_numbers(wavenumber)) %>% 
    dplyr::group_by(group, wavenumber)
  if(!sd_upper_lower){
    spectra <- spectra %>% dplyr::summarise(mean = mean(absorbance), sd = sd(absorbance))
  } else {
    spectra <- spectra %>% dplyr::summarise(mean = mean(absorbance), sd = sd(absorbance), upper = mean + sd, lower = mean - sd)
  }
  if(!sd){
    spectra <- spectra %>% select(-sd)
  } 
  spectra %>% ungroup()
}

# get_mean_spectra(spectra, metadata$sample, sd = T)
# get_mean_spectra(spectra, metadata$sample, smoothed_derivative = T, w = 9)

