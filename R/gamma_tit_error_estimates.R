#' Run compute_error for all dtfs in a gamma titration object
#'
#' @param gamma_tit Gamma titration, 1 row per value of gamma
gamma_tit_error_estimates <- function(
  gamma_tit,
  ref_so,
  q_sc_so,
  sce,
  ref_reconstruction = NULL,
  primary_ref_samples = list(
    c('Unstimulated in vitro - 2h',
      'Unstimulated in vitro - 6h',
      'Unstimulated in vitro - 12h',
      'Unstimulated in vitro - 24h',
      '0.1 ng/ml TNFa - 2h',
      '0.1 ng/ml TNFa - 6h',
      '0.1 ng/ml TNFa - 12h',
      '0.1 ng/ml TNFa - 24h',
      '1 ng/ml TNFa - 2h',
      '1 ng/ml TNFa - 6h',
      '1 ng/ml TNFa - 12h',
      '1 ng/ml TNFa - 24h',
      '10 ng/ml TNFa - 2h',
      '10 ng/ml TNFa - 6h',
      '10 ng/ml TNFa - 12h',
      '10 ng/ml TNFa - 24h',
      '1 ng/ml IFNy - 2h',
      '1 ng/ml IFNy - 6h',
      '1 ng/ml IFNy - 12h',
      '1 ng/ml IFNy - 24h',
      '10 ng/ml IFNy - 2h',
      '10 ng/ml IFNy - 6h',
      '10 ng/ml IFNy - 12h',
      '10 ng/ml IFNy - 24h',
      '100 ng/ml IFNy - 2h',
      '100 ng/ml IFNy - 6h',
      '100 ng/ml IFNy - 12h',
      '100 ng/ml IFNy - 24h')
  ),
  # exclusion_affinity_bms =
  #   list(
  #     ## Low risk of confusion in 6369
  #     c('100 ng/ml IFNy - 6h',
  #       '100 ng/ml IFNy 0.1 ng/ml TNFa - 6h'),
  #     ## High risk of confusion in 6369
  #     c('100 ng/ml IFNy - 24h',
  #       '100 ng/ml IFNy 0.1 ng/ml TNFa - 24h'),
  #     c('100 ng/ml IFNy - 24h',
  #       '100 ng/ml IFNy 0.1 ng/ml TNFa - 12h'),
  #     c('10 ng/ml IFNy - 24h',
  #       '10 ng/ml IFNy 0.1 ng/ml TNFa - 24h'),
  #     c('10 ng/ml IFNy - 24h',
  #       '10 ng/ml IFNy 0.1 ng/ml TNFa - 12h')
  # ),
  exclusion_affinity_bms = c(
    compile_all_pairs(c(
        'Unstimulated in vitro - 2h',
        'Unstimulated in vitro - 6h' ,
        'Unstimulated in vitro - 12h',
        'Unstimulated in vitro - 24h')),
    # compile_all_pairs(c(
    #     '100 ng/ml IFNy - 24h',
    #     '10 ng/ml IFNy - 24h',
    #     '1 ng/ml IFNy - 24h',
    #     '100 ng/ml IFNy 0.1 ng/ml TNFa - 24h',
    #     '100 ng/ml IFNy 0.1 ng/ml TNFa - 12h',
    #     '10 ng/ml IFNy 0.1 ng/ml TNFa - 24h',
    #     '10 ng/ml IFNy 0.1 ng/ml TNFa - 12h',
    #     '1 ng/ml IFNy 0.1 ng/ml TNFa - 24h',
    #     '1 ng/ml IFNy 0.1 ng/ml TNFa - 12h')),
    NULL
  ),
  min_Nhood_sizes = as.integer(c(0, 5, 25, 50, 100)),
  do_debug = FALSE,
  ref_weights = 'none',
  verbose = TRUE) {

  gamma_tit <-
    gamma_tit %>%
    ## Assume these columns are constant
    dplyr::distinct(across(
        any_of(c('center_ref', 'center_query', 'log_gamma'))),
      .keep_all = T) %>%
    dplyr::select(-any_of(c('PV_cosine_sim', 'u_dtf', 'N_PV_umap'))) %>%
    { . }

  if (do_debug) {
    gamma_tit <-
      gamma_tit[sample(1:min(nrow(gamma_tit), 2L)), ] %>%
      dplyr::arrange(log_gamma)
  }

  ## Local helper function, fetch the data frame (coordinates +
  ## metadata) associated with a particular value for log_gamma
  fetch_dtf <- function(log_gamma) {
    gamma_tit %>%
      dplyr::filter(log_gamma == .env[['log_gamma']]) %>%
      dplyr::pull(dtf) %>%
      purrr::pluck(1)
  }

  gen_nrs <- function(ltb = 1e-1, log_gamma) {
    l_dtf <- fetch_dtf(log_gamma = log_gamma)
    if (maartenutils::null_dat(l_dtf)) return(NULL)

    if (!maartenutils::null_dat(ref_reconstruction)) {
      ref_reconstruction <- ref_reconstruction %>%
        dplyr::filter(log_gamma == .env[['log_gamma']])
    }

    nrs <- NhoodRefSim$new(
      sce = sce,
      query = SummarizedExperiment::colData(sce)$exp[1],
      ref_experiment = ref_so@meta.data$experiment[1],
      ref_sa =
        ref_so@meta.data %>%
        dplyr::select(any_of(c('sample_name')), stim_group, duration,
          condition_name, matches('conc|dilution')) %>%
        # { set_rownames(., tolower(.$sample_name)) } %>%
        order_duration() %>%
        dplyr::select(-stim_group) %>%
        { . },
      ref_reconstruction = ref_reconstruction,
      ref_CF_clustering = NULL,
      MEV_q = 0,
      min_aff_q = 0,
      cluster_h_q = 0.0,
      exclusion_affinity_bms = exclusion_affinity_bms,
      primary_ref_samples = primary_ref_samples,
      ltb = ltb,
      dtf = l_dtf
    )
    return(nrs)
  }

  error_table <-
    purrr::map_dfr(unique(gamma_tit$log_gamma), function(log_gamma) {
      nrs <- gen_nrs(log_gamma = log_gamma)
      if (null_dat(nrs)) return(NULL)
      out <-
        min_Nhood_sizes %>%
        { .[. <= max(nrs$Nhood_size)] } %>%
        purrr::map_dfr(function(min_Nhood_size) {
          nrs$compute_Nhood_stats(min_Nhood_size = min_Nhood_size)

          ## Determine which ltb are acceptably discriminatory
          ## (i.e. yield a median NN-probability of at least .1 and
          ## aren't too high such that only a handful of Nhoods will
          ## be close enough to any reference sample), such that we
          ## won't have to screen a lot of ltbs that are not going
          ## to be discriminatory anyway
          min_ltb <-
            which(nrs$NN_probs$med_NN_prob >= .10) %>%
            { .[1] } %>%
            { nrs$NN_probs$ltb[.] } %>%
            factor_to_numeric()

          max_ltb <-
            which(nrs$NN_probs$frac_gamma_Nhoods_retained >= .10) %>%
            { .[length(.)] } %>%
            { nrs$NN_probs$ltb[.] } %>%
            factor_to_numeric()

          settings <-
            as.vector(outer(c(1, 2.5, 5, 7.5), 10^(-1:6))) %>%
            { .[. >= min_ltb & . <= max_ltb] } %>%
            rev() %>%
            {
              tidyr::expand_grid(
                ltb = .,
                ref_weights = ref_weights
              )
            }

          purrr::pmap_dfr(settings, function(ltb, ref_weights) {
            agg_error <- compute_error(
              obj = nrs,
              ltb = ltb,
              min_Nhood_size = min_Nhood_size,
              ref_weights = ref_weights
            )
            if (is.null(agg_error))
              return(NULL)

            mean_error <-
              list(
                'mean_error' = mean(as.matrix(agg_error$EM), 
                  na.rm = T),
                'mean_error_CI_l' = CI_l(agg_error$EM),
                'mean_error_CI_h' = CI_h(agg_error$EM)
              )
            signal_spec_error <-
              list(
                apply(agg_error$EM, 2, mean, na.rm = T),
                append_name(apply(agg_error$EM, 2, CI_l), 
                  '_CI_l'),
                append_name(apply(agg_error$EM, 2, CI_h), 
                  '_CI_h')
              ) %>%
              unlist()
            mean_IZ_error <-
              colnames(agg_error$IZ) %>%
              auto_name() %>%
              purrr::map(~1-mean(as.integer(agg_error$IZ[, .x]),
                  na.rm = T)) %>%
              { setNames(., paste0(names(.), '_IZ')) } %>%
              { . }

            idx <- which.min(abs(ltb -
                factor_to_numeric(nrs$NN_probs$ltb)))

            N_non_zero_SC_conditions <-
              tryCatch({
                sum(rowSums(agg_error$WM) > 0)
              }, error = function(e) { 0 })

            out <-
              c(list(
                'log_gamma' = log_gamma,
                'ref_weights' = ref_weights,
                'min_Nhood_size' = min_Nhood_size,
                'ltb' = ltb,
                'N_non_zero_SC_conditions' =
                  N_non_zero_SC_conditions,
                'frac_gamma_Nhoods_retained' =
                  nrs$NN_probs[[idx, 'frac_gamma_Nhoods_retained']],
                'frac_gamma_cells_retained' =
                  nrs$NN_probs[[idx, 'frac_gamma_cells_retained']],
                'frac_size_Nhoods_retained' =
                  nrs$NN_probs[[idx, 'frac_size_Nhoods_retained']],
                'frac_size_cells_retained' =
                  nrs$NN_probs[[idx, 'frac_size_cells_retained']],
                'frac_Nhoods_retained' =
                  nrs$NN_probs[[idx, 'frac_Nhoods_retained']],
                'frac_cells_retained' =
                  nrs$NN_probs[[idx, 'frac_cells_retained']]
              )) %>%
              c(mean_IZ_error) %>%
              c(as.list(mean_error)) %>%
              c(as.list(signal_spec_error)) %>%
              { . }

            return(out)
          })
        })
      return(out)
    }) %>%
    { . }

  if (F && !maartenutils::null_dat(error_table)) {
    error_table <-
      error_table %>%
      dplyr::mutate(
        global_error_min = min(mean_error, na.rm = T),
        # global_tnfa_conc_min = min(mean_tnfa_conc, na.rm = T),
        optimal_error = abs(mean_error - global_error_min) < 1e-15
      ) %>%
      { . }

    if (verbose && 'tnf_conc_IZ' %in% colnames(error_table)) {
      print(error_table[which.min(error_table$tnf_conc_IZ), ])
    }
  }

  attr(error_table, 'exclusion_affinity_bms') <-
    exclusion_affinity_bms
  attr(error_table, 'primary_ref_samples') <-
    primary_ref_samples
  return(error_table)
}
