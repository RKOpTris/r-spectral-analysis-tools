#convert a decimal value to a hexidecimal value
dec_to_hex <- function(x, xmax){
  as.hexmode(round(x / xmax * 255, 0))
}
