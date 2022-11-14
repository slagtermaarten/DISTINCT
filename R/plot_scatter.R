#' 
#'
#' @export
plot_Nhood_milo <- function(...)
  UseMethod('plot_Nhood_milo')


#' 
#'
#' @export
plot_Nhood_milo.DISTINCT <- function(
  obj, 
  o_fn = file.path(get_out_dir(obj),
    glue::glue('Nhood_milo{get_plot_id(obj)}.pdf'))) {

  ## The first 64 entries in dtf are for the reference samples, the
  ## remaining are for the single cell Nhoods. Let's verify this
  stopifnot(unlist(nhoodIndex(obj$sce)) == obj$dtf$rn[65:nrow(obj$dtf)])

  dtf <- obj$dtf

  test_non_NA <-
    dplyr::select(dtf, matches('UMAP')) %>%
    is.na() %>%
    all()

  if (!'UMAP1' %in% colnames(dtf) || all(is.na(dtf$UMAP1))) {
    dtf <- add_umap(dtf, column_selector = matches('CF|NLPC'))
  }
  dtf <- order_concentration(dtf)
  dtf$tnf_conc <- factor(dtf$tnf_conc)

  set_coord_ranges <- gen_coord_ranges_fun(dtf)

  p0 <- plot_dim_reduc(
    dplyr::filter(dtf, experiment %in% obj$ref_experiment),
    ref_experiments = obj$ref_experiment,
    legend_ncol = 1,
    coord_regex = '^UMAP_*(1|2)$'
  ) %>% set_coord_ranges()
  p0 <- p0 + theme_cyto_inf()
  # print_plot_eval(print(p0),
  #   width = 17.4, height = 25, filename = o_fn)

  # class(dtf$tnf_conc)
  p1 <- plot_dim_reduc(
    dplyr::filter(dtf, experiment %in% obj$ref_experiment),
    legend_ncol = 2,
    ref_experiments = obj$ref_experiment,
    coord_regex = '^UMAP_*(1|2)$',
    colour_var = 'tnf_conc'
  ) %>% set_coord_ranges()
  p1 <- p1 + theme_cyto_inf()
  # print_plot_eval(print(p1),
  #   width = 17.4, height = 25, filename = o_fn)

  Nhood_dtf <- dplyr::filter(dtf, !is.na(N))

  milo_plots <- purrr::imap(obj$NH_test_DA, function(.x, .y) {
   out <- plotNhoodGraphDA(obj$sce, .x,
     layout = as.matrix(Nhood_dtf[, c('UMAP1', 'UMAP2')]),
     alpha=0.25) +
     ggtitle(.y) +
     theme_cyto_inf() +
     xlab('') + ylab('') +
     guides(size = 'none', edge_width = 'none')
    set_coord_ranges(out)
  })

  library(patchwork)
  pe <- quote(print((
      p0 + theme(legend.position = 'none') +
      p1 + theme(legend.position = 'right', 
        legend.direction = 'vertical')
    ) / patchwork::wrap_plots(milo_plots) + 
  plot_layout(heights = c(.3, .7))))

  dir.create(dirname(o_fn))
  print_plot_eval(eval(pe),
    width = 17.4, height = 25,
    # width = 17.4, height = 20,
    # width = knitr::opts_chunk$get('fig.width')*2.54,
    # height = knitr::opts_chunk$get('fig.height')*2.54,
    filename = o_fn
  )

  return(o_fn)
}
