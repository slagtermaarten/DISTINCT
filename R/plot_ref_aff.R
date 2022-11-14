#' 
#'
#' @export
plot_ref_aff <- function(...) UseMethod('plot_ref_aff')


#' 
#'
#' @export
plot_ref_aff.DISTINCT <- function(obj, ...) {
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

  do.call(plot_ref_aff, c(list('obj' = nrs), dots))
}


#' 
#'
#' @export
plot_ref_aff.NhoodRefSim <- function(obj, out_dir, plot_id) {
  ## This is the apparent affinity for two samples that are not
  ## distinguishable
  # max_tol_aff <-
  #   aff_M['100 ng/ml IFNy - 24h', '100 ng/ml IFNy 10 ng/ml TNFa - 24h']
  # aff_M[which(abs(aff_M - 12.42) < 0.01, arr.ind = T)]
  # rownames(aff_M)[t(which(abs(aff_M - 12.42) < 0.01, arr.ind = T))]
  # log10(12.42 + 1)

  # p_aff_M <- exp(-self$ltb * as.matrix(dist(self$ref_M)))
  p_aff_M <- as.matrix(dist(obj$ref_M))
  p_aff_M[!is.finite(p_aff_M)] <- 0
  diag(p_aff_M) <- 0
  # p_aff_M <- log10(p_aff_M)
  # max(p_aff_M)
  # p_aff_M <-
  # var(p_aff_M['0.1 ng/ml TNFa - 12h', ])
  tree <- gen_clust_object(p_aff_M)
  tree_split <- gen_tree_split(tree,
    cluster_h = quantile(as.hclust(tree)$height, .9))

  idxs <- which(apply(p_aff_M, 2, var) == 0)
  if (length(idxs) > 0) {
    message('Dropping: ',
      paste(colnames(p_aff_M)[idxs], collapse = ', '))
    p_aff_M <- p_aff_M[-idxs, -idxs]
  }

  # p_aff_M[p_aff_M < max_tol_aff] <- 0
  # p_aff_M <- log10(p_aff_M + 1)
  HM <- gen_HM(
    p_aff_M,
    name = 'Sample affinity',
    cluster_rows = tree,
    cluster_columns = tree,
    row_split = max(c(tree_split, 2)),
    column_split = max(c(tree_split, 2)),
    gap = unit(1, 'mm'),
    show_row_names = T,
    show_column_names = T
  )

  o_fn <-
    file.path(out_dir, glue::glue('ref_aff_unfiltered_HM\\
        {prepend_hyphen(plot_id)}.pdf'))
  print_plot_eval(draw(HM, heatmap_legend_side = 'top'),
    width = 17, height = 17, filename = o_fn)

}
