post_bm_res <- function(bm_res) {
  bm_res <-
    bm_res %>%
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
  return(bm_res)
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


rank_weights <- list(
  'TNFa_CI' =
    list(
      frac_Nhoods_retained_p = 2,
      tnf_conc_CI_h_p = 100,
      tnf_conc_p = 10,
      ifn_conc_p = 10,
      duration_p = 10
    ),
  'TNFa_CI_max_Nhoods' =
    list(
      frac_Nhoods_retained_p = 1,
      tnf_conc_CI_h_p = 100
    ),
  'duration_CI' =
    list(
      frac_Nhoods_retained_p = 5,
      duration_CI_h_p = 100,
      tnf_conc_p = 10,
      ifn_conc_p = 10,
      duration_p = 10
    )
  ) %>%
  purrr::map(function(obj) {
    all_parms <- c('frac_Nhoods_retained_p',
      'duration_CI_l_p', 'tnf_conc_CI_l_p', 'ifn_conc_CI_l_p',
      'duration_p', 'tnf_conc_p', 'ifn_conc_p',
      'duration_CI_h_p', 'tnf_conc_CI_h_p', 'ifn_conc_CI_h_p')
    missing_parms <- setdiff(all_parms, names(obj))
    c(obj, map(auto_name(missing_parms), ~0))
  })


rank_funs <- map(rank_weights, function(penalty_weights) {
  force(penalty_weights)
  out <- function(obj) {
    penalty <-
      setdiff(names(obj), 'frac_Nhoods_retained_p') %>%
      intersect(stringr::str_replace(names(penalty_weights),
          '_p', '')) %>%
      intersect(colnames(obj)) %>%
      auto_name() %>%
      purrr::map_dfc(function(cn) {
        w <- penalty_weights[[glue::glue('{cn}_p')]]
        if (w > 0) {
          return(obj[[cn]] * w)
        } else {
          return(NULL)
        }
      }) %>%
      # purrr::reduce(`*`) %>%
      rowMeans() %>%
      { . }
    exp(-penalty_weights$frac_Nhoods_retained_p * 
      obj$frac_Nhoods_retained) * penalty
  }
  class(out) <- c('rank_fun', class(out))
  return(out)
})


if (F) {
  #' 
  #'
  #'
  as.character.rank_fun <- function(obj, sep='-', sub_print_names = F) {
    penalty_weights <- environment(obj)$penalty_weights %>%
      purrr::keep(~!is.na(.x) && .x > 0)
    names(penalty_weights) <-
      stringr::str_replace(names(penalty_weights), '_p', '')
    # if (sub_print_names) {
    #   names(penalty_weights) <- DISTINCT_p_names
    # }
    imap_chr(penalty_weights, ~glue('{sep}{.y}={.x}')) %>%
      paste0(collapse = '')
  }
  # glue('{rank_funs[[1]]}')
  # as.character(rank_funs[[1]], ' ')
}


#' Print a human-friendly version of a ranking function 
#'
#' @param obj List of parameters
#'
#' @export
print.rank_fun <- function(obj) {
  penalty_weights <- environment(obj)$penalty_weights %>%
    purrr::keep(~!is.na(.x) && .x > 0)
  cat('DISTINCT benchmark ranking function:\n')
  for (cn in names(penalty_weights)) {
    cn_p <- stringr::str_replace(cn, '_p', '')
    cat('  ', cn_p, ':', penalty_weights[[cn]], '\n')
  }
}
# print(rank_funs[[1]])


sort_bm_res <- function(
  bm_res,
  filtering_fun = filtering_funs[[1]],
  rank_fun = rank_funs[['TNFa_duration']]) {

  if ('gamma_D' %in% colnames(bm_res) && 
      !'ltb' %in% colnames(bm_res)) {
    bm_res <- dplyr::rename(bm_res, ltb = gamma_D)
  }

  bm_res <- filtering_fun(bm_res)
  scores <- rank_fun(bm_res)
  entry_ord <- order(scores, decreasing = FALSE)

  if (FALSE) {
    scores[entry_ord]
    cor(entry_ord, bm_res$mean_error[entry_ord], 
      method = 'spearman')
    cor(entry_ord, bm_res$mean_error, method = 'spearman')
    cor(1:length(entry_ord), bm_res$mean_error[entry_ord], 
      method = 'spearman')
    scores[rev(entry_ord)]
    bm_res[rev(entry_ord), ]
    bm_res[entry_ord[1:5], ]
  }

  bm_res[entry_ord, ]
}


#' Reduce a large set of benchmarking results to 1 setting per
#' combination of tier 1 parameters
#' 
#' @param bm_res A list of benchmarking results
#' @param include_non_hl Also include non-highlighted tier 2 setings.
#' This setting will add he 'hl' column to the function's output.
#' @param filtering_fun Function to pre-filter benchmarking results
#' with
#' @param rank_fun Function to rank benchmarking results with
#'
#' @return An ordered tibble of benchmarking result of class
#' 'DISTINCT_bm_summary'
#'
#' @export
summarize_bm_res <- function(
  bm_res,
  include_non_hl = FALSE,
  mod_names = FALSE,
  filtering_fun = filtering_funs[[1]],
  rank_fun = rank_funs[[1]]) {

  ## Aggregate a bunch of individual benchmarks (each for one set of
  ## tier 1 params)
  bm_res_flat <-
    purrr::imap_dfr(bm_res, function(.x, .y) {
      # settings_meta <- map_dfr(names(bm_res), subset_grid_dtf)
      ## TODO portability
      if (stringr::str_detect(.y, 'kpca')) {
        settings_meta <- subset_grid_dtf(.y, grid_dtf = NH_kpca_grid)
      } else {
        settings_meta <- subset_grid_dtf(.y, grid_dtf = NH_t_grid)
      }
      if (null_dat(settings_meta)) return(NULL)
      if (nrow(settings_meta) > 1) browser()
      dplyr::select(.x, everything()) %>%
        dplyr::mutate(name = .y) %>%
        { dplyr::bind_cols(settings_meta, .) } %>%
        { . }
    })

  if (mod_names) {
    ## Turn ref_weights into a tier 1 param by including it in 'name'
    bm_res_flat$name <- paste0(bm_res_flat$name, '_',
      bm_res_flat$ref_weights)
  }

  ## Sort the aggregate and reduce to one row per tier1 setting
  out <-
    sort_bm_res(bm_res_flat, rank_fun = rank_fun) %>%
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

  class(out) <- c('DISTINCT_bm_summary', class(out))
  return(out)
}


#' Sort a set of tier1 objects by a ranking function and extract the
#' i-th item, along with the error metrics for this item
#'
#' @export
extract_best <- function(
  bm_res = NULL,
  bm_res_sum = NULL,
  filtering_fun = filtering_funs[[1]],
  rank_fun = rank_funs[[1]],
  tier1_rank = 1L) {

  stopifnot(xor(is.null(bm_res), is.null(bm_res_sum)))

  if (is.null(bm_res_sum) && !is.null(bm_res)) {
    bm_res_sum <- summarize_bm_res(
      bm_res,
      rank_fun = rank_fun,
      filtering_fun = filtering_fun
    )
  }

  tier1_rank <- min(tier1_rank, nrow(bm_res_sum))

  optimal_tier1 <- bm_res_sum$name[tier1_rank]

  ## What are the associated tier2 settings?
  optimal_tier2 <-
    bm_res_sum %>%
    dplyr::slice(tier1_rank) %>%
    # dplyr::right_join(tibble(name = optimal_tier1), by = 'name') %>%
    dplyr::select(any_of(c('kernel', 'log_gamma', 'ltb'))) %>%
    # as.list()
    { . }

  agg_error <-
    bm_res_sum %>%
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
  bm_res = NULL,
  bm_res_sum = NULL,
  rank_fun = rank_funs[['TNFa_duration']],
  filtering_fun = filtering_funs[[1]],
  lookup_hyperparams_f = lookup_tier1,
  # settings_grid = exp6369_grid,
  settings_grid = NH_t_grid,
  tier1_rank = 1L,
  plot_id = NULL) {

  stopifnot(xor(is.null(bm_res), is.null(bm_res_sum)))
  stopifnot(!is.null(rank_fun) || !is.function(rank_fun))

  obj <- list()

  if (!is.null(plot_id))
    obj$plot_id <- plot_id

  if (!is.null(bm_res_sum)) {
    obj$bm_res_sum <- bm_res_sum
  } else {
    if (!is.null(bm_res)) {
      obj$bm_res_sum <- summarize_bm_res(
        bm_res,
        rank_fun = rank_fun,
        mod_names = FALSE,
        filtering_fun = filtering_fun
      )
    }
  }

  obj$hyperparams <- extract_best(
    bm_res = bm_res,
    bm_res_sum = obj$bm_res_sum,
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
    as_nrs(obj)
    class(obj)
    saveRDS(obj, file.path(rds_dir, 'DISTINCT_test.rds'))
    # obj <- readRDS(file.path(rds_dir, 'DISTINCT_test.rds'))
  }
  agg_error <- as_nrs.DISTINCT(obj)
  agg_error <- as_nrs(obj)$agg_Nhood_error()
  Nhood_importance_weights <- colSums(agg_error$CM)

  obj$dtf$Nhood_importance_weights <-
    c(
      ## This assumes the reference to be the first rows in dtf
      # rep(NA_real_, table(obj$dtf$experiment)[1]),
      rep(NA_real_, rle(obj$dtf$experiment)$lengths[1]),
      Nhood_importance_weights
    )

  obj$dtf$Nhood_bandwidth_weights <-
    obj$dtf$Nhood_importance_weights / obj$dtf$N

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


#' 
#'
#' @export
as_nrs <- function(...) UseMethod('as_nrs')


#' Convert a DISTINCT object to a nrs (NhoodRefSim sample) object 
#'
#' @export
as_nrs.DISTINCT <- function(
  obj,
  ltb = obj$hyperparams$tier2$ltb,
  primary_ref_samples = 
    eval(formals(gamma_tit_error_estimates)$primary_ref_samples),
  exclusion_affinity_bms = 
    eval(formals(gamma_tit_error_estimates)$exclusion_affinity_bms)) {

  nrs <- NhoodRefSim$new(
    dtf = obj$dtf,
    sce = obj$sce,
    query = SummarizedExperiment::colData(obj$sce)$exp[1],
    primary_ref_samples = primary_ref_samples,
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


#' Make all DISTINCT plots
#'
#' @export
plot.DISTINCT <- function(
  obj,
  ltb = obj$hyperparams$tier2$ltb,
  out_dir = get_out_dir(obj),
  plot_id = get_plot_id(obj)) {

  dir.create(out_dir, showWarnings = FALSE)
  nrs <- as_nrs(obj, ltb = ltb)

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
    # nrs$agg_Nhood_error(scale_by_max_aff = TRUE)
    # nrs$plot_prediction_error_mat(
    #   out_dir = out_dir,
    #   fn_app = glue::glue('{plot_id}-scale_by_max_aff')
    # )
    nrs$agg_Nhood_error(scale_by_max_aff = FALSE)
    nrs$plot_prediction_error_mat(
      fn_app = glue::glue('{plot_id}'),
      out_dir = out_dir
    )

    min_Nhood_size = 50L
    nrs$agg_Nhood_error(scale_by_max_aff = FALSE,
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
  bm_res_sum,
  max_rows = 100L,
  vis_vars = c('ifn_conc', 'tnf_conc', 'duration',
    'frac_Nhoods_retained', 'frac_cells_retained'),
  grid_dtf = NH_t_grid,
  # vis_vars = c('frac_Nhoods_retained', 'frac_cells_retained'),
  of = file.path(out_dir,
    glue::glue('optimal_settings.pdf'))) {

  if ('hl' %in% colnames(bm_res_sum) && any(bm_res_sum$hl ==
      'Highlight setting')) {
    bm_res_sum <- bm_res_sum %>%
      dplyr::filter(hl == 'Highlight setting')
  }
  bm_res_sum <- bm_res_sum[!is.na(bm_res_sum$name), ]
  if (!is.null(max_rows)) {
    bm_res_sum <- bm_res_sum[1:min(max_rows, nrow(bm_res_sum)), ]
  }
  # bm_res_sum$frac_Nhoods_retained

  vis_vars <- intersect(vis_vars, colnames(bm_res_sum))

  ## Table of tier 1 parameters
  settings_comp <-
    bm_res_sum %>%
    # dplyr::pull(bm_res_sum, name) %>%
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
    # dplyr::inner_join(bm_res_sum, by = 'name') %>%
    dplyr::select(-any_of(c('NH_gamma_titration_embedding_obj',
          'NH_Nhood_error_obj',
          'name', 'global_tnfa_conc_min', 'd', 'N_PV'))) %>%
    # bind_cols(bm_res_sum)
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
  nrs <- as_nrs(ds)
  nrs$agg_Nhood_error(
    ltb = ds$hyperparams$tier2$ltb,
    min_Nhood_size = ds$hyperparams$tier2$min_Nhood_size
  )
  expected <- ds$hyperparams$agg_error$mean_error
  observed <- mean(as.matrix(nrs$agg_error$EM), na.rm = T)
  list('exp' = expected, 'obs' = observed,
    'delta' = abs(expected-observed))
}


#' Print an object of the DISTINCT_bm_summary class
#' 
#' @param obj A DISTINCT_bm_summary object
#' 
#' @export 
print.DISTINCT_bm_summary <- function(obj) {
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
#' @param bm_res_sum A DISTINCT_bm_summary object
#' @param of Output file, where to store the plot
#'
#' @export
settings_scatter.DISTINCT_bm_summary <- function(bm_res_sum,
  of = NULL) {

  if (is.null(bm_res_sum$hl)) {
    bm_res_sum$hl <- 'none'
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
          bm_res_sum %>%
          dplyr::arrange(desc(hl)) %>%
          ggplot(aes_string(x = x_var, y = y_var,
            label = 'hl_index', alpha = 'hl', size = 'hl',
            colour = 'hl')) +
          xlab(DISTINCT_p_names[x_var]) +
          ylab(DISTINCT_p_names[y_var]) +
          # geom_point(data = bm_res_flat, color = 'grey50', alpha = .5) +
          # geom_point(data = bm_res_sum, color = 'indianred3', size = 3, alpha = .5)
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
            data = dplyr::filter(bm_res_sum, !is.na(hl_index)),
            min.segment.length = 0,
            max.iter = 1e9,
            force = 10,
            max.overlaps = 20L,
            size = 2, colour = 'grey2', show.legend = FALSE) +
          theme()
      } else if (type == 'boxplot') {
        p <-
          bm_res_sum %>%
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


#' Deprecated
#'
#'
make_all_plots <- function(query) {
  out_dir <- file.path(Sys.getenv('img_dir'),
    glue::glue('exp{query}_TRANSACT'))
  dir.create(out_dir, showWarnings = F)

  bm_res <-
    tar_read_regex(
      regex = glue::glue('NH_Nhood_error_.*_{query}\\
        .*(mallow|linear|rbf|cosine).*'),
      min_date = '2022-07-08 12:00:00 CEST'
    )

  bm_res <-
    purrr::map(bm_res, function(x)
      dplyr::filter(x, !is.na(duration_IZ)) %>%
      dplyr::mutate(global_error_min =
        min(mean_error[!is.na(duration_IZ)]))
    )

  for (i in 1:length(rank_funs)) {
    l_out_dir <- file.path(out_dir, names(rank_funs)[i])
    dir.create(l_out_dir, showWarnings = F)

    bm_res_sum <-
      summarize_bm_res(
        bm_res,
        rank_fun = rank_funs[[i]],
        filtering_fun = filtering_funs[[1]],
        include_non_hl = TRUE
      )

    D_obj <- DISTINCT(bm_res_sum = bm_res_sum)

    hl_settings_heatmap(
      bm_res_sum = dplyr::filter(bm_res_sum,
        hl == 'Highlight setting'),
      max_rows = 100L,
      vis_vars = c(
        # 'duration_IZ',
        'duration',
        # 'tnf_conc_IZ',
        'tnf_conc',
        # 'ifn_conc_IZ',
        'ifn_conc',
        # 'frac_conditions_retained',
        'frac_Nhoods_retained',
        'frac_cells_retained',
        NULL),
      of = file.path(l_out_dir,
        glue::glue('optimal_settings_{names(rank_funs)[i]}.pdf'))
    )

    walk(unique(bm_res_sum$kernel), function(kernel) {
      t_dat <- dplyr::filter(bm_res_sum, kernel == .env[['kernel']])
      settings_scatter(
        bm_res_sum = t_dat,
        of = file.path(l_out_dir,
          glue::glue('error_marginal_{names(rank_funs)[i]}\\
            {make_flag(kernel)}.pdf'))
      )
    })

    plot_gamma_titration_Nhood_retention(
      D_obj,
      plot_id = names(rank_funs)[i],
      out_dir = l_out_dir
    )

    plot_prediction_error_mat(
      D_obj,
      plot_id = names(rank_funs)[i],
      out_dir = l_out_dir
    )

    if (interactive() && !test_rendering())
      source('~/MirjamHoekstra/R/init.R')
    plot_ref_mds(
      D_obj,
      include_query = TRUE,
      plot_id = names(rank_funs)[i],
      out_dir = l_out_dir
    )

  }
}


plot_reference_compression <- function(obj) {
  sa <- extract_sa(extract_ref_dtf(obj$dtf),
    meta_fields = c('condition_name', 'stim_group', 'duration')
  )
  p_dat <- ref_reconstruction %>%
    dplyr::left_join(sa, by = 'condition_name') %>%
    dplyr::mutate(log_gamma = round(log_gamma, 2)) %>%
    dplyr::filter(log_gamma %in% unique(log_gamma)[1:3])

  p <-
    p_dat %>%
    ggplot(aes(x = 1/NLPC, y = 1/CF, colour = stim_group, shape = duration)) +
    geom_abline(slope = 1, intercept = 0) +
    geom_point(alpha = 1) +
    scale_colour_stim_group(p_dat) +
    xlab('Mean distance to unstimulated\nsamples in NLPC space') +
    ylab('Mean distance to unstimulated\nsamples in CF space') +
    scale_shape_duration(p_dat) +
    guides(color = guide_legend(ncol = 2)) +
    facet_wrap(~log_gamma)
  print_plot_eval(print(p),
    width = 17.4, height = 10,
    filename = file.path(out_dir,
      glue::glue('ref_compression.pdf')))
}
