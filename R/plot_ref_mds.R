#' 
#'
#' @export
plot_ref_mds <- function(...) UseMethod('plot_ref_mds')


#' 
#'
#' @export
plot_ref_mds.DISTINCT <- function(
  obj, ltb = obj$hyperparams$tier2$ltb, ...) {

  nrs <- as_nrs(obj, ltb = ltb)
  dots <- list(...)
  if (is.null(dots$out_dir)) {
    dots$out_dir <- get_out_dir(obj)
    dir.create(dots$out_dir, showWarnings = FALSE)
  }
  if (is.null(dots$plot_id)) {
    dots$plot_id <- get_plot_id(obj)
  }

  do.call(plot_ref_mds, c(list('obj' = nrs), dots))
}


