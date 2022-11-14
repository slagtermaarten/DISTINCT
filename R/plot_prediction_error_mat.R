#' 
#'
#' @export
plot_prediction_error_mat.DISTINCT <- function(obj, ...) {
  nrs <- as_nrs(obj)
  nrs$agg_Nhood_error()
  dots <- list(...)
  if (is.null(dots$out_dir)) {
    dots$out_dir <- get_out_dir(obj)
    dir.create(dots$out_dir, showWarnings = FALSE)
  }
  if (is.null(dots$plot_id)) {
    dots$plot_id <- get_plot_id(obj)
  }

  do.call(plot_prediction_error_mat, c(list('obj' = nrs), dots))
}
