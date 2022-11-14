#' 
#'
#' @export
get_out_dir <- function(...) UseMethod('get_out_dir')


#' 
#'
#' @export
get_out_dir.DISTINCT <- function(obj) {
  file.path(Sys.getenv('img_dir'),
    with(c(obj, unlist(obj$hyperparams)),
      glue::glue('DISTINCT-{obj$hyperparams$tier1}')))
}


#' 
#'
#' @export
get_plot_id <- function(...) UseMethod('get_plot_id')


#' 
#'
#' @export
get_plot_id.DISTINCT <- function(obj) {
  if (is.null(obj$plot_id)) {
    out <- with(c(obj, obj$hyperparams$tier1, as.list(obj$hyperparams$tier2)),
      glue::glue('{make_flag(ref_experiment)}\\
        {make_flag(log_gamma)}{make_flag(ltb)}'))
  } else {
    out <- obj$plot_id
  }
  return(out)
}


#' 
#'
#' @export
plot_prediction_error_mat.DISTINCT <- function(obj, ...) {
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

  do.call(plot_prediction_error_mat, c(list('obj' = nrs), dots))
}



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
plot_aff_M <- function(...) UseMethod('plot_aff_M')


#' 
#'
#' @export
plot_aff_M.DISTINCT <- function(
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

  do.call(plot_aff_M, c(list('obj' = nrs), dots))
}


#' Resimulate from a probability mass function and return metadata
#' (sample annotation) corresponding to simulated classes
#'
#'
resim_pmf <- function(X, sa, N_resamples = 1e3, N_draws = 1e2) {
  stopifnot(rowSums(X) - 1 <= 1e-6)
  # stopifnot(nrow(sa) == nrow(XR))
  if (length(N_resamples) == 1) {
    N_resamples <- rep(N_resamples, nrow(X))
  }
  if (length(N_draws) == 1) {
    N_draws <- rep(N_draws, nrow(X))
  }
  purrr::map(1:nrow(X), function(i) {
    XR <- stats::rmultinom(N_resamples[i], N_draws[i], prob = X[i, ])
    XR <- XR / colSums(XR)
    # colSums(XR)
    pred_sa <- t(XR) %*% sa
  })
}


#' Run agg_Nhood_error for all dtfs in a gamma titration object
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
            nrs$agg_Nhood_error(
              ltb = ltb,
              min_Nhood_size = min_Nhood_size,
              ref_weights = ref_weights
            )
            if (is.null(nrs$agg_error))
              return(NULL)

            mean_error <-
              list(
                'mean_error' = mean(as.matrix(nrs$agg_error$EM), 
                  na.rm = T),
                'mean_error_CI_l' = CI_l(nrs$agg_error$EM),
                'mean_error_CI_h' = CI_h(nrs$agg_error$EM)
              )
            signal_spec_error <-
              list(
                apply(nrs$agg_error$EM, 2, mean, na.rm = T),
                append_name(apply(nrs$agg_error$EM, 2, CI_l), 
                  '_CI_l'),
                append_name(apply(nrs$agg_error$EM, 2, CI_h), 
                  '_CI_h')
              ) %>%
              unlist()
            # mean_tnfa_conc <-
            #   mean(as.matrix(nrs$agg_error$PM)[, 'tnf_conc'])
            mean_IZ_error <-
              colnames(nrs$agg_error$IZ) %>%
              auto_name() %>%
              purrr::map(~1-mean(as.integer(nrs$agg_error$IZ[, .x]),
                  na.rm = T)) %>%
              { setNames(., paste0(names(.), '_IZ')) } %>%
              { . }

            idx <- which.min(abs(ltb -
                factor_to_numeric(nrs$NN_probs$ltb)))

            N_non_zero_SC_conditions <-
              tryCatch({
                sum(rowSums(nrs$agg_error$WM) > 0)
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


#' 
#'
#' @export
plot_gamma_titration_Nhood_retention <- function(...)
  UseMethod('plot_gamma_titration_Nhood_retention')


#' 
#'
#' @export
plot_gamma_titration_Nhood_retention.NhoodRefSim <-
  function(obj, plot_id, out_dir, include_error = F) {

  obj$compute_Nhood_stats()

  median_probs <-
    obj$NN_probs %>%
    group_by(ltb) %>%
    dplyr::summarize(
      med = median(NN_prob, na.rm = T),
      q25 = quantile(NN_prob, .25, na.rm = T),
      q75 = quantile(NN_prob, .75, na.rm = T)
    ) %>%
    { . }

  plots <- list(
    p1 =
      ggplot(obj$NN_probs, aes(x = ltb,
          group = 1,
          y = 100 * frac_Nhoods_retained)) +
      # ggplot(obj$NN_probs, aes(x = ltb, y = N)) +
      geom_line(alpha = .5, colour = 'grey50', size = 1) +
      geom_point(alpha = .5, colour = 'grey50', size = 1) +
      ylab('% identifiable\nneighbourhoods') +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      ),

    p2 =
      ggplot(obj$NN_probs, aes(x = ltb, y = NN_prob,
          group = Nhood_idx)) +
      geom_line(alpha = .5, colour = 'grey50') +
      geom_line(data = median_probs,
        mapping = aes(x = ltb, y = med, group = 1),
        inherit.aes = FALSE, alpha = .7, size = 2, colour = 'indianred3') +
      geom_ribbon(data = median_probs,
        mapping = aes(x = ltb, ymin = q25, ymax = q75, group = 1),
        inherit.aes = FALSE, alpha = .2, fill = 'indianred3') +
      ylab('Assigned weight of\nnearest bulk reference sample') +
      xlab(expression(gamma)) +
      rotate_x_labels(45)
  )

  if (!include_error) {
    heights <- c(.2, .8)
  } else {
    plots <- c(list(p0), plots)
    heights <- c(.2, .2, .8)
  }
  p <- wrap_plots(plots, heights = heights, ncol = 1)

  if (T) {
    o_fn <- file.path(
      out_dir,
      glue::glue('gamma_titration_Nhood_retention\\
        {prepend_hyphen(plot_id)}.pdf')
    )
    print_plot_eval(print(p),
      # width = 17.4, height = 15, filename = o_fn)
      width = 12, height = 15, filename = o_fn)
  } else {
    return(p)
  }
}


#' 
#'
#' @export
plot_gamma_titration_Nhood_retention.DISTINCT <- function(obj, ...) {
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

  do.call(plot_gamma_titration_Nhood_retention, c(list('obj' = nrs), dots))
}


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
