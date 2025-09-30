INT <- function(data, pheno, k = 0.375){
  missing_index <- which(is.na(data[[pheno]]))
  
  if(length(missing_index) == 0){
    # All values are complete: use the whole data
    n <- nrow(data)
    r <- rank(data[[pheno]])
    data$int <- qnorm((r - k) / (n - 2 * k + 1))
  } else {
    # Some missing values exist: compute for complete cases only
    data.complete <- data[-missing_index, ]
    n <- nrow(data.complete)
    r <- rank(data.complete[[pheno]])
    data$int <- NA  # initialize
    data$int[-missing_index] <- qnorm((r - k) / (n - 2 * k + 1))
  }
  
  return(data)
}
