#' @description Find index and value of the nearest numeric value in a numeric or character vector
#'
#' @param character/numeric
#'
#' @return data.frame
#'
#' @examples
#' nearest(paste0("sample", 1, ".", seq(2, 10, 2)), 1.4)
#' nearest(seq(2.1, 10.1, 2), c(4, 10))

nearest <- function(vector, val){
  stopifnot(is.numeric(vector) | is.character(vector))
  if(is.numeric(vector)){
    ind <- sapply(val, function(x){which.min(abs(x - vector))})
  } else {
    ind <- sapply(val, function(x){which.min(abs(x - strip_non_numbers(vector)))})
  }
  data.frame(
    value_in = val,
    nearest = vector[ind],
    index = ind,
    stringsAsFactors = F
  )
}
