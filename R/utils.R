#' Search for values v in all columns of dtf
#'
#' @param v A character vector
#' @param dtf A data.frame
#'
#' @return The name of the column in dtf that contains the most
#' matches with the entries in v
#'
find_match_col <- function(v, dtf) {
  stopifnot(is.data.frame(dtf))
  map(dtf, ~sum(tolower(v) %in% tolower(.x))) %>%
    which.max() %>%
    names()
}


#' Compute the Euclidean distances between  the columns of two
#' matrices
#'
#'
matrix_euclidean <- function(X, Y) {
  stopifnot(nrow(X) == nrow(Y))
  co <- matrix(0, nrow = ncol(X), ncol = ncol(Y))
  dimnames(co) = list(colnames(X), colnames(Y))
  for (i in 1:ncol(X)) {
    for (j in 1:ncol(Y)) {
      # i = 1; j = 1
      co[i,j] <- sqrt(sum((X[, i] - Y[, j])^2))
    }
  }
  return(co)
}
