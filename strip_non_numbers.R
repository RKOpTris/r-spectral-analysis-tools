# remove non-number characters from a vector, useful for column names, sample numbers, etc.
# BUT ONLY WHEN THERE IS ONE COMPLETE NUMBER
strip_non_numbers <- function(vector){
  gsub("[^0-9.-]", "", vector) %>% as.numeric()
}

# (samples <- paste0("sample", 1:5))
# samples %>% strip_non_numbers
#
# not useful for:
# (samples <- paste0("9999sample", 1:5))
# samples %>% strip_non_numbers()
