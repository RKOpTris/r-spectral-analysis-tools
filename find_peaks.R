find_peaks <- function(x, wns, type = "both", return = "both", wn_hi = NULL, wn_lo = NULL, nudge = 0.025, nups = 5, ndowns = 1, ...){
  nudge <- diff(range(x)) * nudge
  if(!is.null(wn_hi)){
    wn_hi <- nearest(wns, wn_hi)$index
  } else {
    wn_hi <- length(wns)
  }
  if(!is.null(wn_lo)){
    wn_lo <- nearest(wns, wn_lo)$index
  } else {
    wn_lo <- 1
  }
  x <- x[wn_lo:wn_hi]
  wns <- wns[wn_lo:wn_hi]
  #need to implement return
  if(type %in% c("both", "peak")){
    spec_diff_peaks <- pracma::findpeaks(x, nups = nups, ...) %>% data.frame()
    names(spec_diff_peaks) <- c("abs", "peak_ind", "peak_start", "peak_end")
    spec_diff_peaks$peak_wavenumber <- wns[spec_diff_peaks$peak_ind]
    spec_diff_peaks$start_wavenumber <- wns[spec_diff_peaks$peak_start]
    spec_diff_peaks$end_wavenumber <- wns[spec_diff_peaks$peak_end]
    spec_diff_peaks$start_absorbance <- x[spec_diff_peaks$peak_start]
    spec_diff_peaks$end_absorbance <- x[spec_diff_peaks$peak_end]
    spec_diff_peaks$peak_type <- "peak"
    spec_diff_peaks$nudge <- spec_diff_peaks$abs + nudge
    peaks <- spec_diff_peaks
  }
  if(type %in% c("both", "trough")){
    spec_diff_peaks <- pracma::findpeaks(-1 * x, nups = ndowns, ...) %>% data.frame()
    names(spec_diff_peaks) <- c("abs", "peak_ind", "peak_start", "peak_end")
    spec_diff_peaks$peak_wavenumber <- wns[spec_diff_peaks$peak_ind]
    spec_diff_peaks$start_wavenumber <- wns[spec_diff_peaks$peak_start]
    spec_diff_peaks$end_wavenumber <- wns[spec_diff_peaks$peak_end]
    spec_diff_peaks$start_absorbance <- x[spec_diff_peaks$peak_start]
    spec_diff_peaks$end_absorbance <- x[spec_diff_peaks$peak_end]
    spec_diff_peaks$peak_type <- "trough"
    spec_diff_peaks$nudge <- (spec_diff_peaks$abs + nudge) * -1
    spec_diff_peaks$abs <- spec_diff_peaks$abs * -1
  }
  if(type == "both"){
    bind_rows(peaks, spec_diff_peaks) %>% arrange(peak_wavenumber)
  } else {
    spec_diff_peaks
  }
}


