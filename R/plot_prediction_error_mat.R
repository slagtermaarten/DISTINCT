plot_prediction_error_mat <- function(...) 
  UseMethod('plot_prediction_error_mat')


#' 
#'
#' @export
plot_prediction_error_mat.DISTINCT <- function(obj, ...) {
  nrs <- as_nrs(obj)
  nrs$compute_error()
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


#' 
#'
#' @export
plot_prediction_error_mat.NhoodRefSim <- function(obj,
  plot_id = '', out_dir = Sys.getenv('img_dir')) {

  stopifnot(!is.null(obj$agg_error))
  EM <- obj$agg_error$EM
  IZ <- obj$agg_error$IZ %>%
    { .[, colnames(.) != 'duration'] }
  IZ_n <- IZ
  IZ_n[IZ == F] <- 1.0
  IZ_n[IZ == T] <- 0.0
  colnames(IZ_n) <- paste0(colnames(IZ_n), '_bin')


  obj <- SummarizedExperiment::SummarizedExperiment(
    assays = list(mean = cbind(EM, IZ_n)),
    rowData = extract_sa(
      as.data.frame(SummarizedExperiment::colData(obj$sce)),
      meta_fields = c('condition_i', 'stim_group', 'duration',
        'frozen')) %>%
    set_rownames(NULL) %>%
    dplyr::distinct() %>%
    dplyr::arrange(condition_i) %>%
    dplyr::select(-condition_i) %>%
    set_rownames(rownames(.)) %>%
    { . }
  )

  plot_prediction_error_mat(obj, plot_id = plot_id, out_dir = out_dir)
}
