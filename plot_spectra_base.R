# using base plot, quickly plot spectra from a tidy spectra object
# group is the name of the column to subset the data by, the default is group as this is created by default using get_mean_spectra()
# include_groups can be a vector of indices for which groups to include
# wn_hi is the maximum wavenumber to plot
# wn_lo is the minimum wavenumber to plot
# col_vector is a vector of colours used to colour the groups, the length of this vector should be the same as the number of groups
# alpha is to alter the opacity of the plotted spectra
# overlay allows the spectra to be overplotted on other data using par(new = T)
# ... send further parameters to lines()

plot_spectra_base <- function(tidy_spectra, group = "group", include_groups = "all", wn_hi = NULL, wn_lo = NULL, col_vector, alpha = 1, overlay = F, ...){
  range_x <- rev(range(tidy_spectra$wavenumber))
  if(!is.null(wn_hi)){
    range_x[1] <- wn_hi
  }
  if(!is.null(wn_lo)){
    range_x[2] <- wn_lo
  }
  plot_deriv_subset <- tidy_spectra %>% filter(wavenumber <= range_x[1] & wavenumber >= range_x[2])
  if(include_groups != "all" & is.numeric(include_groups)){
    groups <- unique(plot_deriv_subset[[group]])
    groups <- groups[include_groups]
    plot_deriv_subset <- plot_deriv_subset[plot_deriv_subset$group %in% groups, ]
  }
  if(overlay){
    par(new = T)
  }
  plot(mean ~ wavenumber, plot_deriv_subset[plot_deriv_subset[[group]] == unique(plot_deriv_subset[[group]])[1], ], 
       type = "n",
       xlim = range_x,
       ylim = inRange(range(plot_deriv_subset$mean), c(-0.1, 1.0)),
       las = 1,
       axes = F,
       xlab = expression(paste("Wavenumber (cm"^"-1"*")")),
       ylab = "",
       cex.lab = 2)
  abline(h = 0, col = "lightgrey")
  
  col_vector <- paste0(tolower(col_vector), dec_to_hex(alpha, 1))
  
  for(i in 1:(length(unique(plot_deriv_subset[[group]])))){
    lines(mean ~ wavenumber, plot_deriv_subset[plot_deriv_subset[[group]] == unique(plot_deriv_subset[[group]])[i], ],
          col = col_vector[i],
          ...)
  }
  internal_axis_ticks(1, 1)
  mtext(2, 1.5, at = inRange(plot_deriv_subset$mean, 0.5), text = "Relative absorbance", cex = 2)
}
