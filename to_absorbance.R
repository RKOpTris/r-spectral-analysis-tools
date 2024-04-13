# convert FTIR transmission (%) measurements to absorbance using "tidy" spectra with a "wavenumber" column

to_absorbance <- function(df, ...){
  wavenumbers <- df[["wavenumber"]]
  df$wavenumber <- NULL
  df <- lapply(df, trans_to_abs, ...)
  data.frame(wavenumber = wavenumbers, df)
}