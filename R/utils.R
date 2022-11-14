#' 
#'
#' @export
get_out_dir <- function(...) UseMethod('get_out_dir')


#' 
#' 
#' @param obj DISTINCT object
#'
#' @export
get_out_dir.DISTINCT <- function(obj) {
  file.path(Sys.getenv('img_dir'),
    with(c(obj, unlist(obj$hyperparams)),
      glue::glue('DISTINCT-{obj$hyperparams$tier1}')))
}


#' 
#'
#' @export
get_plot_id <- function(...) UseMethod('get_plot_id')


#' Extract a plot ID from a DISTINCT object
#' 
#' @param obj DISTINCT object
#'
#' @export
get_plot_id.DISTINCT <- function(obj) {
  if (is.null(obj$plot_id)) {
    out <- with(c(obj, obj$hyperparams$tier1, as.list(obj$hyperparams$tier2)),
      glue::glue('{make_flag(ref_experiment)}\\
        {make_flag(log_gamma)}{make_flag(ltb)}'))
  } else {
    out <- obj$plot_id
  }
  return(out)
}


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


#' Compute the Euclidean distances between the columns of two
#' matrices. The two matrices need to have the same number of rows
#' ('features')
#'
#' @param X feature matrix
#' @param Y feature matrix
#'
#' @return \code{matrix} with dimensions [X_n, Y_n]
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


#' Resimulate from a probability mass function and return metadata
#' (sample annotation) corresponding to simulated classes
#'
#' @param X Probability mass function
#' @param sa Sample annotation
#' @param N_resamples Number of resamples per draw
#' @param N_draws Number of draws
resim_pmf <- function(X, sa, N_resamples = 1e3, N_draws = 1e2) {
  stopifnot(rowSums(X) - 1 <= 1e-6)
  # stopifnot(nrow(sa) == nrow(XR))
  if (length(N_resamples) == 1) {
    N_resamples <- rep(N_resamples, nrow(X))
  }
  if (length(N_draws) == 1) {
    N_draws <- rep(N_draws, nrow(X))
  }
  purrr::map(1:nrow(X), function(i) {
    XR <- stats::rmultinom(N_resamples[i], N_draws[i], prob = X[i, ])
    XR <- XR / colSums(XR)
    # colSums(XR)
    pred_sa <- t(XR) %*% sa
    return(pred_sa)
  })
}
