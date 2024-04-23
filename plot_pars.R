plot_pars <- function(r = 1, c = 1, s = 5.1, w = 4.1, n = 1.1, e = 1.1, pty = "m", reset = F){
  if(reset){
    mfrow()
    mar()
    par()
  } else {
    mfrow(r = r, c = c)
    mar(s = s, w = w, n = n, e = e)
    par(pty = pty)
  }
}