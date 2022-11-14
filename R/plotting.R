#' Plot a dimension reduction of integrated data
#'
#' @param dtf Data.frame of latent variables + meta.data + UMAP
#' coordinates (one sample per row)
#'
#'
plot_dim_reduc <- function(
  dtf = umap_emb,
  colour_var = 'stim_group',
  use_stim_group_cols = TRUE,
  shape_var = 'duration',
  coord_regex = 'UMAP',
  legend_ncol = 1,
  zoom_range = 1 == 1,
  plot_range = NULL,
  zoom_name = NULL,
  filter_samples = 1 == 1,
  ordering_code = NULL,
  filter_name = NULL,
  ref_experiments = c('4910', '5029', '6434', '6623', '5029-6434'),
  point_alpha = .5,
  plot_ref_background = T,
  width = 17.4,
  height = 20,
  annotation_code = NULL,
  fn_app = '',
  print_to_file = F) {

  library(rlang)
  library(dplyr)
  library(ggplot2)

  if (maartenutils::null_dat(dtf)) return(NULL)

  if (!'experiment' %in% colnames(dtf)) {
    if ('exp' %in% colnames(dtf)) {
      dtf$experiment <- dtf$exp
    } else {
      stop('Experiment column not present')
    }
  }
  experiments <- gen_exp_string(unique(dtf$experiment))

  query_experiments <-
    setdiff(unique(dtf$experiment), ref_experiments)
  # library(rlang)

  # dtf <- denorm_regressors(dtf)
  dtf <- order_duration(dtf)

  fg_dtf <- 
    dtf %>%
    dplyr::filter(!!enquo(filter_samples)) %>%
    dplyr::arrange(
      across(any_of(c('experiment', 'duration',
            'ifn_conc', 'tnf_conc', 'sn_dilution')))
    )

  if (!is.null(ordering_code)) {
    eval(ordering_code)
  }

  stopifnot(is.factor(fg_dtf$duration))

  if (nrow(fg_dtf) == nrow(dtf)) {
    filter_name <- NULL
  }

  if (!is.null(shape_var) && 
      (!shape_var %in% colnames(fg_dtf) ||
      length(unique(fg_dtf[[shape_var]])) <= 1))
    shape_var <- NULL

  coord_vars <- stringr::str_subset(colnames(fg_dtf), coord_regex)
  # is.na(dtf[, coord_vars])

  zoom_samples <-
    dplyr::filter(fg_dtf, !!enquo(zoom_range))
  if (!is.null(plot_range)) {
    plot_range <- zoom_samples %>%
      dplyr::summarize(across(any_of(coord_vars), range))
  }
  if (nrow(zoom_samples) == nrow(fg_dtf)) {
    zoom_name <- NULL
  }

  ## Setup background dat for each query/non-ref experiment
  ref_background_dat <-
    intersect(query_experiments, fg_dtf$experiment) %>%
    purrr::map_dfr(function(e) {
      dplyr::filter(dtf,
        experiment %in% ref_experiments & experiment != e) %>%
        dplyr::mutate(experiment = e)
    }) %>%
    dplyr::arrange(across(matches('duration')))

  # fg_dtf$stim_group <- droplevels(fg_dtf$stim_group)
  extract_col_var <- function(x) stringr::str_replace(x, '_simp', '')
  aes_vec <- list(coord_vars[1], coord_vars[2],
      extract_col_var(colour_var), shape_var)

  if (aes_vec %>%
    purrr::discard(is.null) %>%
    map(length) %>%
    { all(. != 1) }) {
  }

  p <- ggplot(fg_dtf,
    aes_string(
      x = coord_vars[1],
      y = coord_vars[2],
      colour = extract_col_var(colour_var),
      shape = shape_var
    )
  )

  if (plot_ref_background &&
      !maartenutils::null_dat(ref_background_dat)) {
    ## 'Fool' ggplot facet_wrap by duplicating reference data labelled
    ## as query data

    p <- p +
      geom_path(
        data = ref_background_dat,
        mapping = aes(group = stim_group),
        colour = 'grey90', alpha = .5, show.legend = F) +
      geom_point(
        data = ref_background_dat,
        colour = 'grey90', alpha = .5)
  }

  if (plot_ref_background &&
      !maartenutils::null_dat(query_background_dat)) {
    query_background_dat <-
      intersect(query_experiments, fg_dtf$experiment) %>%
      purrr::map_dfr(function(e) {
      dplyr::filter(dtf, experiment %in% ref_experiments) %>%
        dplyr::mutate(experiment = e)
      }) %>%
      dplyr::arrange(across(matches('duration')))

    p <- p +
      geom_path(
        data = query_background_dat, aes(group = stim_group),
        colour = 'grey90', alpha = .2, show.legend = F) +
      geom_point(
        data = query_background_dat,
        colour = 'grey90', alpha = .2)
  }

  line_dat <- fg_dtf %>%
    dplyr::filter(experiment %in% ref_experiments) %>%
    dplyr::arrange(experiment,
      across(matches('duration|conc|dilution')))
  if (!null_dat(line_dat)) {
    p <- p +
      geom_path(
        data = line_dat, aes(group = stim_group),
        alpha = point_alpha, show.legend = F
      ) +
      geom_point(data = line_dat, alpha = .2, size = 4)
  }

  p <- p +
    geom_point(alpha = point_alpha) +
    theme_cyto_inf() +
    theme_tabula_rasa +
    # ggtitle(colour_var) +
    coord_cartesian(
      xlim = plot_range[[coord_vars[1]]],
      ylim = plot_range[[coord_vars[2]]]
    ) +
    guides(colour = guide_legend(ncol = legend_ncol))

  if (length(unique(fg_dtf$experiment)) > 1) {
    p <- p + facet_wrap(~experiment)
  }

  if (F) {
    print_plot_eval(
      {
        print(p)
      },
      width = 17.4, height = 15,
      filename = file.path(Sys.getenv('img_dir'),
        glue::glue('test.png')))
  }

  if (!is.null(shape_var) && shape_var == 'duration') {
    p <- p + scale_shape_duration(fg_dtf)
  }

  col_name <- switch(extract_col_var(colour_var),
    'stim_group' = 'Stimulus',
    'ifn_conc' = 'IFNy concentration',
    'tnf_conc' = 'TNFa concentration',
    NULL)
  if (!is.null(colour_var)) {
    if (use_stim_group_cols && colour_var == 'stim_group') {
      p <- p + scale_colour_stim_group(fg_dtf)
    } else if (colour_var == 'stim_group_simp' ||
      is.factor(fg_dtf[[colour_var]])) {
      p <- p + scale_colour_discrete(name = col_name)
    } else if (grepl(colour_var, 'MAE|error')) {
      library(viridis)
      p <- p +
        scale_colour_viridis(name = colour_var) +
        guides(colour = guide_colourbar())
    } else {
      p <- p +
        scale_colour_gradient2(name = colour_var, mid = 'grey90')
    }
  }

  if (!is.null(annotation_code)) {
    if (is_expression(annotation_code)) {
      eval(annotation_code)
    } else if (is.list(annotation_code)) {
      for (ac in annotation_code) {
        if (is_expression(ac)) {
          eval(ac)
        } else {
          rlang::halt(
            paste('Do not know what to do with ac of class: ',
              class(ac)))
        }
      }
    }
  }

  if (print_to_file) {
    print_umap(
      p = p, dtf = dtf, 
      coord_regex = coord_regex,
      experiments = experiments, 
      colour_var = colour_var,
      filter_name = filter_name, 
      zoom_name = zoom_name,
      fn_app = fn_app, 
      width = width, 
      height = height
    )
  } else {
    return(p)
  }
}



