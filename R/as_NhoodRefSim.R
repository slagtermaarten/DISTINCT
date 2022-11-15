#' 
#'
#' @export
as_nrs <- function(...) UseMethod('as_NhoodRefSim')
as_NhoodRefSim <- function(...) UseMethod('as_NhoodRefSim')


#' Convert a DISTINCT object to a nrs (NhoodRefSim sample) object 
#'
#' @export
as_NhoodRefSim.DISTINCT <- function(
  obj,
  ltb = obj$hyperparams$tier2$ltb,
  primary_ref_samples = 
    obj$primary_ref_samples %||%
    eval(formals(compute_error_estimates)$primary_ref_samples),
  exclusion_affinity_bms = 
    obj$exclusion_affinity_bms %||%
    eval(formals(compute_error_estimates)$exclusion_affinity_bms)) {

  nrs <- NhoodRefSim$new(
    dtf = obj$dtf,
    sce = obj$sce,
    query = SummarizedExperiment::colData(obj$sce)$exp[1],
    primary_ref_samples = primary_ref_samples,
    ref_factors = obj$ref_factors,
    exclusion_affinity_bms = exclusion_affinity_bms,
    ref_sa =
      obj$ref_so@meta.data %>%
      dplyr::select(any_of(c('sample_name')), stim_group,
        duration, condition_name, matches('conc|dilution')) %>%
      # { set_rownames(., tolower(.$sample_name)) } %>%
      order_duration() %>%
      dplyr::select(-stim_group) %>%
      { . },
    ref_experiment = obj$ref_so@meta.data$experiment[1],
    ltb = ltb
  )

  return(nrs)
}
