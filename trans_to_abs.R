trans_to_abs <- function(vector, infinite_as = NA){
  out <- sapply(vector, function(x){ -log10(x / 100) })
  out[is.infinite(out)] <- infinite_as
  out
}