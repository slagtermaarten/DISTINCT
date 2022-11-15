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


#' 
#'
#' @export
plot_ref_mds.NhoodRefSim <- function(
  obj, out_dir = get_out_dir(obj), plot_id = get_plot_id(obj),
  include_query = F) {

  if (!include_query) {
    fit <- cmdscale(dist(obj$ref_M))
    sa <- obj$ref_sa %>%
      recover_stim_group()
    exp_v <- rep(obj$ref_experiment, nrow(obj$ref_M))
  } else {
    fit <- cmdscale(dist(rbind(obj$ref_M, obj$query_M)))
    ref_sa <- obj$ref_sa %>%
      recover_stim_group()

    query_sa <-
      obj$dtf %>%
      dplyr::filter(experiment == obj$query) %>%
      dplyr::select(any_of(colnames(ref_sa)), any_of(c('N')),
        matches('weight|_N'))
    # query_sa <-
    #   SummarizedExperiment::colData(obj$sce) %>%
    #   as.data.frame() %>%
    #   set_rownames(NULL) %>%
    #   extract_sa(meta_fields = c(colnames(ref_sa), 'condition_i')) %>%
    #   { . }

    exp_v <- c(
      rep(obj$ref_experiment, nrow(obj$ref_M)),
      rep(obj$query, nrow(obj$query_M))
    )

    sa <- harmonize_bind_rows(ref_sa, query_sa)
  }

  p_dat <-
    tibble(
      CMD1 = fit[, 1], CMD2 = fit[, 2],
      experiment = exp_v
    ) %>%
    dplyr::bind_cols(sa) %>%
    dplyr::mutate(tnf_conc = factor(tnf_conc)) %>%
    dplyr::mutate(ifn_conc = factor(ifn_conc))

  p1 <-
    p_dat %>%
    dplyr::filter(experiment == obj$ref_experiment) %>%
    plot_dim_reduc(
      coord_regex = '^CMD',
      colour_var = 'ifn_conc'
    ) +
    ggtitle('IFNy') +
    theme(
      legend.direction = 'vertical',
      legend.position = 'right'
    )

  p2 <-
    p_dat %>%
    dplyr::filter(experiment == obj$ref_experiment) %>%
    plot_dim_reduc(
      coord_regex = '^CMD',
      colour_var = 'tnf_conc'
    ) +
    theme(
      legend.direction = 'vertical',
      legend.position = 'right'
    ) +
    ggtitle('TNFa') +
    guides(shape = 'none')

  if (!include_query) {
    print_plot_eval(
      print(p1 + p2 + plot_layout(guides = 'auto')),
      width = 17.4, height = 10,
      filename = file.path(out_dir,
        glue::glue('mds{prepend_hyphen(plot_id)}.pdf')))
  } else {
    p3 <-
      p_dat %>%
      # dplyr::filter(experiment != obj$ref_experiment) %>%
      plot_dim_reduc(
        coord_regex = '^CMD',
        filter_samples = !experiment %in% c('6434', '5029-6434'),
        ordering_code =
          rlang::expr({fg_dtf <- dplyr::arrange(fg_dtf, N)}),
        colour_var = 'N'
      ) +
      theme(
        legend.direction = 'horizontal',
        legend.position = 'bottom'
      ) +
      guides(color = guide_colourbar()) +
      scale_colour_viridis_c() +
      ggtitle('Nhood size') +
      guides(shape = 'none')
    p31 <-
      p_dat %>%
      # dplyr::filter(experiment != obj$ref_experiment) %>%
      plot_dim_reduc(
        coord_regex = '^CMD',
        filter_samples = !experiment %in% c('6434', '5029-6434'),
        ordering_code =
          rlang::expr({fg_dtf <- dplyr::arrange(fg_dtf, Nhood_bandwidth_weights)}),
        colour_var = 'Nhood_bandwidth_weights'
      ) +
      theme(
        legend.direction = 'horizontal',
        legend.position = 'bottom'
      ) +
      guides(color = guide_colourbar()) +
      scale_colour_viridis_c() +
      ggtitle('Bandwidth weight') +
      guides(shape = 'none')
    p4 <-
      p_dat %>%
      plot_dim_reduc(
        coord_regex = '^CMD',
        filter_samples = !experiment %in% c('6434', '5029-6434'),
        ordering_code =
          rlang::expr({fg_dtf <- dplyr::arrange(fg_dtf, Nhood_importance_weights)}),
        colour_var = 'Nhood_importance_weights'
      ) +
      theme(
        legend.direction = 'horizontal',
        legend.position = 'bottom'
      ) +
      guides(color = guide_colourbar()) +
      scale_colour_viridis_c() +
      ggtitle('Total weight') +
      guides(shape = 'none')
    print_plot_eval(
      print((p1 + p2) / (p3 + p31 + p4) + plot_layout(guides = 'auto')),
      width = 17.4, height = 20,
      filename = file.path(out_dir,
        glue::glue('mds{prepend_hyphen(plot_id)}.pdf')))
  }

}


