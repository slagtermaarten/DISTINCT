post_error_table <- function(error_tables) {
  error_tables <-
    error_tables %>%
    # { .[!stringr::str_detect(names(.), 'kpca_(comb|6434)')] } %>%
    # { .[!stringr::str_detect(names(.), 'kpca')] } %>%
    purrr::discard(is.null) %>%
    purrr::imap(function(x, on) {
      out <- x %>%
        dplyr::filter(!is.na(duration_IZ)) %>%
        dplyr::mutate(global_error_min =
          min(mean_error[!is.na(duration_IZ)])) %>%
        dplyr::mutate(integration_method = ifelse(
            stringr::str_detect(on, 'kpca'),
            'kpca', 'transact'))
      if (is.null(out$ref_weights)) {
        out$ref_weights <- 'none'
      }
      return(out)
    })
  return(error_tables)
}


filtering_funs <- list(
  'basic' = function(obj, min_frac_Nhoods_retained = .20) {
    tryCatch({
      obj %>%
        as_tibble() %>%
        # dplyr::filter(tnf_conc_IZ <= .25) %>%
        # dplyr::filter(frac_conditions_retained >= .99) %>%
        dplyr::filter(N_non_zero_SC_conditions ==
          max(N_non_zero_SC_conditions, na.rm = T)) %>%
        dplyr::filter(
          frac_Nhoods_retained >= min_frac_Nhoods_retained) %>%
        { . }
    }, error = function(e) { browser() })
  }
)


if (F) {
  ## It's late at night and I don't feel like calculus
  ## The larger g, the larger the relative penalty for low frac_Nhoods
  map(auto_name(c(.01, .1, 10)), ~exp(-.x) / exp(-.x * .1))
}


sort_error_tables <- function(
  error_tables,
  filtering_fun = filtering_funs[[1]],
  rank_fun = rank_funs[['TNFa_duration']]) {

  if ('gamma_D' %in% colnames(error_tables) && 
      !'ltb' %in% colnames(error_tables)) {
    error_tables <- dplyr::rename(error_tables, ltb = gamma_D)
  }

  error_tables <- filtering_fun(error_tables)
  scores <- rank_fun(error_tables)
  entry_ord <- order(scores, decreasing = FALSE)

  if (FALSE) {
    scores[entry_ord]
    cor(entry_ord, error_tables$mean_error[entry_ord], 
      method = 'spearman')
    cor(entry_ord, error_tables$mean_error, method = 'spearman')
    cor(1:length(entry_ord), error_tables$mean_error[entry_ord], 
      method = 'spearman')
    scores[rev(entry_ord)]
    error_tables[rev(entry_ord), ]
    error_tables[entry_ord[1:5], ]
  }

  error_tables[entry_ord, ]
}


#' Reduce a large set of benchmarking results to 1 setting per
#' combination of tier 1 parameters
#' 
#' @param error_tables A list of benchmarking results
#' @param include_non_hl Also include non-highlighted tier 2 setings.
#' This setting will add he 'hl' column to the function's output.
#' @param filtering_fun Function to pre-filter benchmarking results
#' with
#' @param rank_fun Function to rank benchmarking results with
#'
#' @return An ordered tibble of benchmarking results of class
#' 'error_tables_summary'
#'
#' @export
summarize_error_tables <- function(
  error_tables,
  include_non_hl = FALSE,
  mod_names = FALSE,
  filtering_fun = filtering_funs[[1]],
  rank_fun = rank_funs[[1]]) {

  ## Aggregate a bunch of individual benchmarks (each for one set of
  ## tier 1 params)
  bm_res_flat <-
    purrr::imap_dfr(error_tables, function(.x, .y) {
      dplyr::select(.x, everything()) %>%
        dplyr::mutate(name = .y) %>%
        { . }
    })

  if (mod_names) {
    ## Turn ref_weights into a tier 1 param by including it in 'name'
    bm_res_flat$name <- paste0(bm_res_flat$name, '_',
      bm_res_flat$ref_weights)
  }

  ## Sort the aggregate and reduce to one row per tier1 setting
  out <-
    sort_error_tables(bm_res_flat, rank_fun = rank_fun) %>%
    dplyr::distinct(name, .keep_all = TRUE) %>%
    dplyr::mutate(hl_index = 1:n())

  if (include_non_hl) {
    bm_res_flat <- dplyr::anti_join(bm_res_flat, out)

    out <-
      dplyr::bind_rows(
        dplyr::mutate(out, hl = 'Highlight setting'),
        dplyr::mutate(bm_res_flat, hl = 'none')
      )
  }

  ## Assume at least one setting retains all settings
  out$frac_conditions_retained <-
    out$N_non_zero_SC_conditions /
    max(out$N_non_zero_SC_conditions, na.rm = T)
  out$frac_conditions_retained[is.na(out$frac_conditions_retained)] <- 0

  class(out) <- c('error_tables_summary', class(out))
  return(out)
}


#' Sort a set of tier1 objects by a ranking function and extract the
#' i-th item, along with the error metrics for this item
#'
#' @export
extract_best <- function(
  error_tables = NULL,
  error_tables_sum = NULL,
  filtering_fun = filtering_funs[[1]],
  rank_fun = rank_funs[[1]],
  tier1_rank = 1L) {

  stopifnot(xor(is.null(error_tables), is.null(error_tables_sum)))

  if (is.null(error_tables_sum) && !is.null(error_tables)) {
    error_tables_sum <- summarize_error_tables(
      error_tables,
      rank_fun = rank_fun,
      filtering_fun = filtering_fun
    )
  }

  tier1_rank <- min(tier1_rank, nrow(error_tables_sum))

  optimal_tier1 <- error_tables_sum$name[tier1_rank]

  ## What are the associated tier2 settings?
  optimal_tier2 <-
    error_tables_sum %>%
    dplyr::slice(tier1_rank) %>%
    # dplyr::right_join(tibble(name = optimal_tier1), by = 'name') %>%
    dplyr::select(any_of(c('kernel', 'log_gamma', 'ltb'))) %>%
    # as.list()
    { . }

  agg_error <-
    error_tables_sum %>%
    dplyr::slice(tier1_rank) %>%
    # dplyr::right_join(tibble(name = optimal_tier1), by = 'name') %>%
    dplyr::select(matches('error|_IZ$|conc')) %>%
    dplyr::select(-any_of('optimal_error')) %>%
    dplyr::select(-matches('compute_error|obj|global')) %>%
    { . }

  list(
    'tier1' = optimal_tier1,
    'tier2' = optimal_tier2,
    'agg_error' = agg_error
  )
}


#' Assess DISTINCT results, class instantiator
#'
#' This class copies all related objects into one 'basket'/class so
#' that we don't have to worry about mismatched objects down the line
#'
#' @export
DISTINCT <- function(
  error_tables = NULL,
  error_tables_sum = NULL,
  rank_fun = rank_funs[['TNFa_duration']],
  filtering_fun = filtering_funs[[1]],
  lookup_hyperparams_f = lookup_tier1,
  exclusion_affinity_bms = 
    attr(error_tables[[1]], 'exclusion_affinity_bms') %||% 
    eval(formals(compute_error_estimates)$exclusion_affinity_bms),
  primary_ref_samples = 
    attr(error_tables[[1]], 'primary_ref_samples') %||% 
    eval(formals(compute_error_estimates)$primary_ref_samples),
  # settings_grid = exp6369_grid,
  settings_grid = NH_t_grid,

  error_level_computation = 'conditions',
  ref_weights = 'none',
  tier1_rank = 1L,
  plot_id = NULL) {

  stopifnot(xor(is.null(error_tables), is.null(error_tables_sum)))
  stopifnot(!is.null(rank_fun) || !is.function(rank_fun))

  obj <- list()

  if (!is.null(plot_id))
    obj$plot_id <- plot_id

  obj$hyperparams <- extract_best(
    error_tables = error_tables,
    error_tables_sum = error_tables_sum,
    rank_fun = rank_fun,
    filtering_fun = filtering_fun,
    tier1_rank = tier1_rank
  )

  if (F) {
    # str(obj$hyperparams$tier2)
    ## The tier1 settings are encoded as a string, retrieve the
    ## individual covariates from the settings table
    obj$tier1_settings <-
      subset_grid_dtf(
        grid_dtf = settings_grid,
        remove_obj_refs = F,
        name = obj$hyperparams$tier1
      )
  }

  obj$tier1_settings <- lookup_hyperparams_f(obj$hyperparams$tier1)

  ## TODO portability: make this 'lookup' section modular, i.e.
  ## external to this function
  gamma_tit_row <-
    obj$hyperparams$tier1 %>%
    stringr::str_replace(
      'Nhood_error',
      'gamma_titration_embedding'
    ) %>%
    tar_read_raw() %>%
    dplyr::right_join(obj$hyperparams$tier2) %>%
    { . }

  if ('PV_cosine_sim' %in% colnames(gamma_tit_row)) {
    obj$CS_mat <-
      gamma_tit_row %>%
      dplyr::pull(PV_cosine_sim) %>%
      purrr::pluck(1)
  }

  obj$dtf <- gamma_tit_row$dtf[[1]]

  ## TODO portability
  obj$sce <-
    obj$tier1_settings %>%
    dplyr::pull(NH_agg_neighbourhoods_obj) %>%
    purrr::pluck(1) %>%
    as.character() %>%
    tar_read_raw() %>%
    { . }

  obj$query <- obj$tier1_settings$query
  obj$ref_experiment <- obj$tier1_settings$reference

  obj$query_sa <-
    SummarizedExperiment::colData(obj$sce) %>%
    as.data.frame() %>%
    set_rownames(NULL) %>%
    extract_sa(meta_fields = c('condition_name', 'duration',
        'ifn_conc', 'tnf_conc', 'condition_i')) %>%
    dplyr::distinct() %>%
    dplyr::arrange(condition_i) %>%
    numerify_regressors() %>%
    norm_regressors() %>% {
      if (obj$query == '6369') {
        .$tnf_conc <- 0
      }
      .
    } %>%
    # add_binary_regressors(regressor_vars = reg_vars) %>%
    # {
    #   if (all(c('duration', 'tnf_conc', 'ifn_conc') %in%
    #       colnames(.))) {
    #     . <- dplyr::mutate(., duration =
    #       if_else(ifn_conc == 0 & tnf_conc == 0, NA_real_, duration))
    #   }
    #   .
    # } %>%
    # dplyr::select(-any_of(c('duration_bin'))) %>%
    # debug_pipe() %>%
    # {
    #   idxs <-
    #     match(
    #       levels(colData(sce)$condition_name),
    #       as.character(.$condition_name)
    #     )
    #   .[idxs, ] } %>%
    { . }

  obj$ref_so <-
    obj$tier1_settings %>%
    dplyr::pull(ref_so) %>%
    purrr::pluck(1) %>%
    as.character() %>%
    tar_read_raw() %>%
    { . }

  obj$NH_test_DA <-
    obj$tier1_settings %>%
    dplyr::pull(NH_test_DA_obj) %>%
    purrr::pluck(1) %>%
    as.character() %>%
    tar_read_raw() %>%
    { . }

  if ('NH_TRANSACT_reference_reconstruction_obj' %in%
    colnames(obj$tier1_settings)) {
    obj$ref_reconstruction <-
      obj$tier1_settings %>%
      dplyr::pull(NH_TRANSACT_reference_reconstruction_obj) %>%
      purrr::pluck(1) %>%
      as.character() %>%
      tar_read_raw() %>%
      { . }
  }

  ## Just a list, but with a 'DISTINCT flavour'
  class(obj) <- c('DISTINCT', class(obj))

  if (F) {
    obj$sce
    as_NhoodRefSim(obj)
    class(obj)
    saveRDS(obj, file.path(rds_dir, 'DISTINCT_test.rds'))
    # obj <- readRDS(file.path(rds_dir, 'DISTINCT_test.rds'))
  }
  # agg_error <- as_nrs.DISTINCT(obj)
  
  agg_error <- purrr::map(c('none', 'ref_reconstruction'), 
    function(elc) {
    agg_error <- compute_error(
      obj,
      ref_weights = ref_weights,
      error_level_computation = elc, 
      exclusion_affinity_bms = exclusion_affinity_bms,
      primary_ref_samples = primary_ref_samples
    )

    if (FALSE) {
      if ('CM' %in% names(agg_error)) {
        Nhood_importance_weights <- colSums(agg_error$CM)
      
        obj$dtf[[glue::glue('Nhood_importance_weights_{elc}')]] <-
          c(
            ## This assumes the reference to be the first rows in dtf
            # rep(NA_real_, table(obj$dtf$experiment)[1]),
            rep(NA_real_, rle(obj$dtf$experiment)$lengths[1]),
            Nhood_importance_weights
          )
      
        obj$dtf$Nhood_bandwidth_weights <-
          obj$dtf$Nhood_importance_weights / obj$dtf$N
      }
    }
  })

  obj$agg_error <- agg_error

  return(obj)
}


#' 
#'
#' @export
print.DISTINCT <- function(obj) {
  cat('DISTINCT object for experiments:',
    obj$ref_experiment, obj$query, '\n')
  cat(obj$hyperparams$tier1, '\n')
  cat('Tier 2 params:\n')
  for (cn in colnames(obj$hyperparams$tier2)) {
    cat('    ', cn, obj$hyperparams$tier2[[cn]], '\n')
  }
  cat('Prediction metrics:\n')
  for (cn in colnames(obj$hyperparams$agg_error)) {
    cat('    ', cn, obj$hyperparams$agg_error[[cn]], '\n')
  }
}


#' 
#'
#' @export
test_error <- function(...) UseMethod('test_error')


#' Make all DISTINCT plots
#'
#' @export
plot.DISTINCT <- function(
  obj,
  ltb = obj$hyperparams$tier2$ltb,
  out_dir = get_out_dir(obj),
  plot_id = get_plot_id(obj)) {

  dir.create(out_dir, showWarnings = FALSE)
  nrs <- as_NhoodRefSim(obj, ltb = ltb)

  if (F) {
    nrs$plot_gamma_titration_Nhood_retention(
      fn_app = plot_id,
      out_dir = out_dir
    )
  }

  if (F) {
    nrs$plot_ref_aff(fn_app = plot_id, out_dir = out_dir)
  }

  if (F) {
    # nrs$compute_error(scale_by_max_aff = TRUE)
    # nrs$plot_prediction_error_mat(
    #   out_dir = out_dir,
    #   fn_app = glue::glue('{plot_id}-scale_by_max_aff')
    # )
    nrs$compute_error(scale_by_max_aff = FALSE)
    nrs$plot_prediction_error_mat(
      fn_app = glue::glue('{plot_id}'),
      out_dir = out_dir
    )

    min_Nhood_size = 50L
    nrs$compute_error(scale_by_max_aff = FALSE,
      min_Nhood_size = min_Nhood_size)
    nrs$plot_prediction_error_mat(
      fn_app = glue::glue('{plot_id}{make_flag(min_Nhood_size)}'),
      out_dir = out_dir
    )
  }
}


DISTINCT_p_names <- c(
  'mean_error' = 'Mean error',
  'duration_IZ' = 'Duration CI',
  'tnf_conc_IZ' = 'TNFa CI',
  'ifn_conc_IZ' = 'IFNy CI',
  'tnf_conc' = '[TNFa] error',
  'ifn_conc' = '[IFNy] error',
  'duration' = 'Duration error',
  'frac_Nhoods_retained' = 'Included Nhoods',
  'frac_cells_retained' = 'Included cells',
  'frac_conditions_retained' = 'Included SC conditions',
  'ltb' = 'Transfer bandwidth'
)


hl_settings_heatmap <- function(
  error_tables_sum,
  max_rows = 100L,
  vis_vars = c('ifn_conc', 'tnf_conc', 'duration',
    'frac_Nhoods_retained', 'frac_cells_retained'),
  grid_dtf = NH_t_grid,
  # vis_vars = c('frac_Nhoods_retained', 'frac_cells_retained'),
  of = file.path(out_dir,
    glue::glue('optimal_settings.pdf'))) {

  if ('hl' %in% colnames(error_tables_sum) && 
      any(error_tables_sum$hl == 'Highlight setting')) {
    error_tables_sum <- error_tables_sum %>%
      dplyr::filter(hl == 'Highlight setting')
  }
  error_tables_sum <- 
    error_tables_sum[!is.na(error_tables_sum$name), ]
  if (!is.null(max_rows)) {
    error_tables_sum <- 
      error_tables_sum[1:min(max_rows, nrow(error_tables_sum)), ]
  }
  # error_tables_sum$frac_Nhoods_retained

  vis_vars <- intersect(vis_vars, colnames(error_tables_sum))

  ## Table of tier 1 parameters
  settings_comp <-
    error_tables_sum %>%
    # dplyr::pull(error_tables_sum, name) %>%
    # purrr::map_dfr(subset_grid_dtf, grid_dtf = grid_dtf) %>%
    { stopifnot(nrow(.) > 0); . } %>%
    # dplyr::select(-any_of(c('NH_gamma_titration_embedding_obj'))) %>%
    dplyr::mutate(across(where(is.integer), numeric2factor)) %>%
    dplyr::mutate(across(where(is.numeric), numeric2factor)) %>%
    dplyr::mutate(across(where(is.character), factor)) %>%
    dplyr::select(-query) %>%
    # dplyr::mutate(name =
    #   as.character(NH_gamma_titration_embedding_obj)) %>%
    # dplyr::mutate(name =
    #   as.character(NH_Nhood_error_obj)) %>%
    # dplyr::inner_join(error_tables_sum, by = 'name') %>%
    dplyr::select(-any_of(c('NH_gamma_titration_embedding_obj',
          'NH_Nhood_error_obj',
          'name', 'global_tnfa_conc_min', 'd', 'N_PV'))) %>%
    { . }

  ## Remove 'constant' columns
  settings_comp <-
    settings_comp[, which(
      imap_lgl(settings_comp,
        ~length(unique(.x)) > 1 | .y %in% vis_vars)
    )]

  l_param_types <- param_types
  l_param_types$extra <- c('ref_weights', 'integration_method')
  plot_vars <- list()
  for (i in 1:length(param_types)) {
    res <-
      types_to_names(names(param_types)[i]) %>%
      intersect(colnames(settings_comp)) %>%
      maartenutils::auto_name() %>%
      ## Remove variables that may have been mentioned already
      { setdiff(., unlist(plot_vars[1:(i-1)])) } %>%
      { . }
    if (!is.null(res)) plot_vars <- append(plot_vars, res)
    res <- NULL
  }

  set.seed(59)
  ## Create row annotation heatmaps
  HM_anns <-
    purrr::map(plot_vars, function(cns) {
      l_dat <- settings_comp[, cns, drop = F]
      colnames(l_dat) <- change_colnames(colnames(l_dat))
      ann <- rowAnnotation(
        # df = as.data.frame(t(l_dat))
        # df = set_rownames(l_dat, NULL)
        df = as.data.frame(l_dat)
        # , gp = map(cns, ~annotation_par),
        , col = purrr::map(auto_name(cns), function(x) {
            levs <- naturalsort::naturalsort(setdiff(unique(l_dat[[x]]), NA))
            suppressWarnings(RColorBrewer::brewer.pal(length(levs),
                'Paired')) %>%
              { .[1:length(levs)] } %>%
              setNames(levs)
          })
        , annotation_name_gp = annotation_par
      )
      return(ann)
    }) %>%
    purrr::discard(is.null)

  iz_anns <-
    purrr::map(c(vis_vars), function(x) {
      largs <- list(
        anno_empty(border = FALSE,
        which = 'row', width = unit(1, 'cm'))
      )
      names(largs) <- x
      do.call(rowAnnotation, largs)
    })

  HM_anns <- HM_anns %>%
    append(iz_anns) %>%
    { . }

  old_padding <- ht_opt$HEATMAP_LEGEND_PADDING
  ht_opt$HEATMAP_LEGEND_PADDING <- unit(2.8, 'cm')
  # ht_opt$ROW_ANNO_PADDING = unit(5.9, 'cm')
  # ht_opt$COLUMN_ANNO_PADDING = unit(5.9, 'cm')
  # ht_opt$ANNOTATION_LEGEND_PADDING = unit(2, 'cm')
  print_plot_eval(
    {
      draw(purrr::reduce(HM_anns, `+`), merge_legends = T,
        heatmap_legend_side = 'bottom',
        newpage = F, ht_gap = unit(5, 'mm'))
      eps <- 1/(2 *nrow(settings_comp))
      eps <- 1/(nrow(settings_comp))
      eps <- .5
      for (vn in c(vis_vars)) {
        # if (all(is.na(settings_comp[[vn]]))) {
        #   next()
        # }
        decorate_annotation(vn, {
          v <- factor_to_numeric(settings_comp[[vn]])
          x_range <- c(max(min(v), .0), max(max(v), .1))
          if (length(unique(v)) == 1 & v[1] > 0) {
            x_range[1] <- 0
          }
          x_range <- c(0, max(max(v), .1))
          pushViewport(
            viewport(
              # xscale = range(v),
              xscale = x_range,
              yscale = c(1-eps, length(v)+eps)
              # yscale = c(0, nrow(v))
            )
          )

          ## Draw background for each panel
          grid.rect(gp = gpar(fill = 'grey95'))

          ## Draw horizontal grid lines
          integral_x <- seq(0, 5) %>%
            { .[. >= x_range[1] & . <= x_range[2]] }
          x_points <- c(
            # round(x_range[1], 2),
            max(round(x_range[2], 2), 0)
          ) %>% c(integral_x) %>% unique() %>% sort
          if (F) {
            if (length(x_points) > 0) {
              for (j in x_points) {
                linecol <- ifelse(
                  j %in% integral_x, 'grey80', 'grey95')
                linecol <- ifelse(
                  j %in% integral_x && length(x_points) > 2,
                  'grey50', 'grey20')
                grid.lines(
                  x = c(j, j),
                  y = c(0.5, length(v)+.5),
                  # y = c(i-eps, i-eps),
                  gp = gpar(lty = 1, col = linecol),
                  default.units = 'native'
                )
              }
            }
          } else {
            if (length(integral_x) > 0) {
              for (j in integral_x) {
                grid.lines(
                  x = c(j, j),
                  y = c(0.5, length(v)+.5),
                  gp = gpar(lty = 1, col = 'grey80'),
                  default.units = 'native'
                )
              }
            }
          }

          # grid.lines(
          #   x = x_range,
          #   y = rep(length(v)+.5, 2),
          #   # y = c(i-eps, i-eps),
          #   gp = gpar(lty = 1, col = 'grey50'),
          #   default.units = 'native'
          # )

          CI_var <- glue::glue('{vn}_CI_l') %in% colnames(settings_comp)
          if (!CI_var) {
            for (i in 1:length(v)) {
              grid.points(
                x = v[i],
                y = length(v)-i+1,
                size = unit(1, 'mm'),
                pch = 19,
                gp = gpar(
                  col = ifelse(stringr::str_detect(vn, 'retained'),
                    'navyblue', 'indianred3')),
                default.units = 'native'
              )
            }
          } else {
            cn_l <- glue::glue('{vn}_CI_l')
            cn_h <- glue::glue('{vn}_CI_h')
            for (i in 1:length(v)) {
              grid.lines(
                x = c(
                  max(factor_to_numeric(settings_comp[[cn_l]][i]),
                    x_range[1]),
                  min(factor_to_numeric(settings_comp[[cn_h]][i]),
                    x_range[2])
                ),
                # y = c(i, i),
                y = c(length(v)-i+1, length(v)-i+1),
                # y = c(i-eps, i-eps),
                gp = gpar(lty = 1,
                  col = ifelse(stringr::str_detect(vn, 'retained'),
                    'navyblue', 'indianred3')),
                default.units = 'native'
              )
              grid.points(
                x = v[i],
                y = length(v)-i+1,
                size = unit(1, 'mm'),
                pch = 19,
                gp = gpar(
                  col = ifelse(stringr::str_detect(vn, 'retained'),
                    'navyblue', 'indianred3')),
                default.units = 'native'
              )
            }
          }
          # grid.polygon(
          #   x = c(v[, 1], rep(0, nrow(v))),
          #   y = c(nrow(v):1, 1:nrow(v)),
          #   gp = gpar(fill = 'indianred3', col = 'grey50'),
          #   default.units = 'native'
          # )
          grid.xaxis(
            at = x_points,
            gp = gpar(fontsize = 6, col = 'grey50')
          )
          grid.text(DISTINCT_p_names[vn],
            x = sum(x_range)/2, y = -0.04 * nrow(settings_comp),
            # vjust = -2.1,
            # hjust = 1.5,
            hjust = 1,
            gp = gpar(fontsize = 6, col = 'black'),
            rot = 90,
            default.units = 'native'
          )
          popViewport()
        })
      }
    },
    width = 17.4, height = 15, filename = of
  )
  ht_opt$HEATMAP_LEGEND_PADDING <- old_padding

  return(invisible(of))
}


test_DISTINCT_cache <- function(ds) {
  nrs <- as_NhoodRefSim(ds)
  nrs$compute_error(
    ltb = ds$hyperparams$tier2$ltb,
    min_Nhood_size = ds$hyperparams$tier2$min_Nhood_size
  )
  expected <- ds$hyperparams$agg_error$mean_error
  observed <- mean(as.matrix(nrs$agg_error$EM), na.rm = T)
  list('exp' = expected, 'obs' = observed,
    'delta' = abs(expected-observed))
}


#' Print an object of the error_tables_summary class
#' 
#' @param obj A error_tables_summary object
#' 
#' @export 
print.error_tables_summary <- function(obj) {
  cat('DISTINCT benchmark summary:\n')
  cat('  reference: ', obj$reference[1], '\n')
  cat('  query: ', obj$query[1], '\n')
  cat(' ', nrow(obj), 'rows \n')
}


#' Define s3 method
#' 
#' @export
settings_scatter <- function(...) UseMethod('settings_scatter')


#' Plot a scatter of all statistics for a summary of settings
#'
#' @param error_tables_sum A error_tables_summary object
#' @param of Output file, where to store the plot
#'
#' @export
settings_scatter.error_tables_summary <- function(
  error_tables_sum,
  of = NULL) {

  if (is.null(error_tables_sum$hl)) {
    error_tables_sum$hl <- 'none'
  }

  plots <-
    tibble(
      x_var = c(
        #'ltb',
        'duration', 'ifn_conc', 'tnf_conc', 'tnf_conc',
        'tnf_conc', 'tnf_conc'),
      y_var = c(
        #'frac_Nhoods_retained',
        'frac_Nhoods_retained', 'frac_Nhoods_retained',
        'frac_Nhoods_retained', 'mean_error', 'duration', 'ifn_conc'),
      type = c(
        # 'boxplot',
        rep('scatter', 6L))
    ) %>%
    purrr::pmap(function(x_var, y_var, type) {
      if (type == 'scatter') {
        p <-
          error_tables_sum %>%
          dplyr::arrange(desc(hl)) %>%
          ggplot(aes_string(x = x_var, y = y_var,
            label = 'hl_index', alpha = 'hl', size = 'hl',
            colour = 'hl')) +
          xlab(DISTINCT_p_names[x_var]) +
          ylab(DISTINCT_p_names[y_var]) +
          # geom_point(data = bm_res_flat, color = 'grey50', alpha = .5) +
          # geom_point(data = error_tables_sum, color = 'indianred3', size = 3, alpha = .5)
          scale_colour_manual(name = 'Highlighted setting',
            values = c(
              'Highlight setting' = 'indianred3',
              'none' = 'grey70')) +
          scale_size_manual(name = 'Highlighted setting',
            values = c('Highlight setting' = 2, 'none' = 1)) +
          scale_alpha_manual(name = 'Highlighted setting',
            values = c('Highlight setting' = .7, 'none' = .2)) +
          ggrastr::rasterise(geom_point(), dpi = 300) +
          ggrepel::geom_text_repel(
            data = dplyr::filter(error_tables_sum, !is.na(hl_index)),
            min.segment.length = 0,
            max.iter = 1e9,
            force = 10,
            max.overlaps = 20L,
            size = 2, colour = 'grey2', show.legend = FALSE) +
          theme()
      } else if (type == 'boxplot') {
        p <-
          error_tables_sum %>%
          dplyr::mutate(x_var = numeric2factor(.data[[x_var]])) %>%
          dplyr::filter(mean_error <= 1) %>%
          ggplot(aes_string(x = 'x_var', y = y_var)) +
          xlab(DISTINCT_p_names[x_var]) +
          ylab(DISTINCT_p_names[y_var]) +
          geom_violin(draw_quantiles = c(.25, .5, .75)) +
          rotate_x_labels(45) +
          labs(subtitle = 'Mean error <= 1') +
          theme()
      }
      return(p)
    })

  # pw_plot_panel_layout <- function(designs = NULL, ncol = 2, nrow = 3)
  design <-
    'AA
     BC
     DE'
  design <-
    'AB
     CD
     EF'
  p1 <- wrap_plots(plots[1:6], design = design, guides = 'collect')
  # design <-
  #   'AB
  #    CD
  #    EF'
  # p2 <- wrap_plots(plots[7:length(plots)],
  #   design = design, guides = 'collect')
  maartenutils::print_plot_eval(
    {
      print(p1);
      # print(p2)
    },
    width = 17.4, height = 25,
    filename = of)

  return(invisible(of))
}
