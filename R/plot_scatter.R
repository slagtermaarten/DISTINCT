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


#' 
#'
#' @export
plot_Nhood_scatter_pie <- function(...)
  UseMethod('plot_Nhood_scatter_pie')


#' 
#'
#' @export
plot_Nhood_scatter_pie.DISTINCT <- function(
  obj, 
  separate_duration = TRUE,
  o_fn = file.path(get_out_dir(obj),
    glue::glue('scatter_pie{get_plot_id(obj)}.pdf'))) {

  dtf <- obj$dtf
  ref_experiment = ds$dtf$exp[1]

  test_non_NA <-
    dplyr::select(dtf, matches('UMAP')) %>%
    is.na() %>%
    all()

  if (!'UMAP1' %in% colnames(dtf) || all(is.na(dtf$UMAP1))) {
    dtf <- add_umap(dtf, column_selector = matches('CF|NLPC'))
  }
  dtf <- order_concentration(dtf)
  dtf$tnf_conc <- factor(dtf$tnf_conc)
  dtf$ifn_conc <- factor(dtf$ifn_conc)

  set_coord_ranges <- gen_coord_ranges_fun(dtf)

  if (!separate_duration) {
    p0 <- plot_dim_reduc(
      dplyr::filter(dtf, experiment %in% ref_experiment),
      ref_experiments = ref_experiment,
      legend_ncol = 1,
      coord_regex = '^UMAP_*(1|2)$'
    ) %>% set_coord_ranges()
  } else {
    p0 <- plot_dim_reduc(
      dplyr::filter(dtf, experiment %in% ref_experiment),
      ref_experiments = ref_experiment,
      legend_ncol = 1,
      colour_var = 'ifn_conc',
      coord_regex = '^UMAP_*(1|2)$'
    ) %>% set_coord_ranges()
  }
  p0 <- p0 + theme_cyto_inf()
  # print_plot_eval(print(p0),
  #   width = 17.4, height = 25, filename = o_fn)

  # class(dtf$tnf_conc)
  p1 <- plot_dim_reduc(
    dplyr::filter(dtf, experiment %in% ref_experiment),
    legend_ncol = 2,
    ref_experiments = ref_experiment,
    coord_regex = '^UMAP_*(1|2)$',
    colour_var = 'tnf_conc'
  ) %>% set_coord_ranges()
  p1 <- p1 + theme_cyto_inf()
  # print_plot_eval(print(p1),
  #   width = 17.4, height = 25, filename = o_fn)

  Nhood_dtf <- dplyr::filter(dtf, !is.na(N))
  CN_levs <- obj$query_sa$condition_name
  # CN_levs <-
  #   get_obj('q_sc_so')@meta.data %>%
  #   # dtf %>%
  #   extract_query_dtf() %>%
  #   order_condition_name() %>%
  #   pull(condition_name) %>%
  #   levels
  # CN_levs
  idxs <- which(stringr::str_detect(colnames(Nhood_dtf), 'CN'))
  colnames(Nhood_dtf)[idxs] <- as.character(CN_levs)

  if (!separate_duration) {
    p2 <- ggplot(Nhood_dtf) +
      scatterpie::geom_scatterpie(
        data = Nhood_dtf,
        mapping = aes(x = UMAP1, y = UMAP2),
        cols = CN_levs) +
      theme(legend.position = 'right', legend.direction = 'vertical')
    p2 <- set_coord_ranges(p2)
  } else {
    l_condition_names <- colnames(Nhood_dtf) %>%
      stringr::str_subset('\\d+h$')
    durations <- map_chr(l_condition_names, 
      ~strsplit(.x, ' - ')[[1]][[2]])
    max_radius = 0
    for (ld in unique(durations)) {
      v <- rowSums(Nhood_dtf[, l_condition_names[durations == ld]])
      Nhood_dtf[[glue::glue('radius_{ld}')]] <- v
      max_radius <- max(c(v, max_radius))
    }
    for (ld in unique(durations)) {
      Nhood_dtf[[glue::glue('radius_{ld}')]] <- 
        Nhood_dtf[[glue::glue('radius_{ld}')]] / (3 * max_radius)
    }
    p2_plots <- map(unique(durations), function(ld) {
      p2 <- ggplot(Nhood_dtf) +
        scatterpie::geom_scatterpie(
          data = Nhood_dtf,
          mapping = aes(x = UMAP1, y = UMAP2, 
            r = .data[[glue::glue('radius_{ld}')]]),
          color = NA,
          cols = l_condition_names[durations == ld]
        ) +
        theme(legend.position = 'top', legend.direction =
          'vertical') +
        scale_colour_discrete(name = 'Single-cell condition name') + 
        coord_fixed() +
        scatterpie::geom_scatterpie_legend(
          Nhood_dtf[[glue::glue('radius_{ld}')]], 
          x=unname(quantile(Nhood_dtf$UMAP1, .5)), 
          y=quantile(Nhood_dtf$UMAP2, .5))
      p2 <- set_coord_ranges(p2)
      return(p2)
    })
    p2 <- wrap_plots(p2_plots)
  }

  library(patchwork)
  ref_bg_col <- "grey95"
  query_bg_col <- "white"
  ref_bg_col <- "white"
  ref_theme <- theme(legend.position = 'top', 
          legend.background = element_rect(fill = "transparent", 
            colour = "transparent"),
          plot.background = element_rect(fill = ref_bg_col),
          legend.direction = 'vertical', legend.box = 'horizontal')
  bulk_panels <- (p0 + ref_theme) | (p1 + ref_theme)
        # plot_layout(heights = c(.1, .9), guides = 'auto')
  t_grob <- patchwork::wrap_elements(
    grobTree(
      rectGrob(gp=gpar(fill=ref_bg_col, col= NA)), 
      textGrob('Bulk reference')
    )
  )
  qt_grob <- patchwork::wrap_elements(
    grobTree(
      rectGrob(gp=gpar(fill=query_bg_col, col= NA)), 
      textGrob('Single-cell query')
    )
  )
  pe <- quote(print((t_grob / bulk_panels / qt_grob / p2 + 
        plot_layout(heights = c(0.05, .3, 0.05, .7)))))
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
