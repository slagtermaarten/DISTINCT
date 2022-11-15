compute_error <- function(...) UseMethod('compute_error')


compute_error.DISTINCT <- function(
  obj, 
  ltb = obj$ltb %||% 10,
  exclusion_affinity_bms = obj$exclusion_affinity_bms %||%
    eval(formals(compute_error_estimates)$exclusion_affinity_bms),
  primary_ref_samples = obj$primary_ref_samples %||%
    eval(formals(compute_error_estimates)$primary_ref_samples),
  reg_vars = obj$reg_vars %||% c('duration', 'ifn_conc', 'tnf_conc'),
  by_var = obj$by_var %||% 'condition_name',
  min_Nhood_size = obj$min_Nhood_size %||% 0,
  ref_weights = obj$ref_weights %||% 'none',
  include_bs_pres = obj$include_bs_pres %||% FALSE,
  error_level_computation = 'none',
  verbose = FALSE) {

  nrs <- as_nrs(obj)

  agg_error <- compute_error(
    obj = nrs, 
    ltb = ltb,
    exclusion_affinity_bms = exclusion_affinity_bms,
    primary_ref_samples = primary_ref_samples,
    reg_vars = reg_vars,
    by_var = by_var,
    min_Nhood_size = min_Nhood_size,
    ref_weights = ref_weights,
    include_bs_pres = include_bs_pres,
    error_level_computation = error_level_computation,
    verbose = verbose
  )

  return(agg_error)
}


compute_error.NhoodRefSim = function(
  obj,
  # reg_vars = c('duration',
  #   'ifn_conc', 'ifn_conc_bin', 'tnf_conc', 'tnf_conc_bin'),
  ltb = obj$ltb,
  exclusion_affinity_bms = obj$exclusion_affinity_bms %||%
    eval(formals(compute_error_estimates)$exclusion_affinity_bms),
  primary_ref_samples = obj$primary_ref_samples %||%
    eval(formals(compute_error_estimates)$primary_ref_samples),
  ref_factors = obj$ref_reconstruction %||% NULL,
  reg_vars = c('duration', 'ifn_conc', 'tnf_conc'),
  by_var = 'condition_name',
  min_Nhood_size = NULL,
  scale_by_max_aff = FALSE,
  ref_weights = 'none',
  include_bs_pres = FALSE,
  error_level_computation = 'conditions',
  verbose = F) {

  ref_weights <- match.arg(ref_weights,
    choices = c('none', 'uniqueness', 'ref_reconstruction'))

  library(SingleCellExperiment)

  ref_dtf <- extract_ref_dtf(obj$dtf)
  ref_M <- extract_CF_M(ref_dtf)
  query_dtf <- extract_query_dtf(obj$dtf)
  # query_M <- extract_CF_M(query_dtf)

  # obj$compute_affM(ltb = ltb)
  query_sa <-
    SummarizedExperiment::colData(obj$sce) %>%
    as.data.frame() %>%
    set_rownames(NULL) %>%
    extract_sa(meta_fields = c(by_var, reg_vars, 'condition_i')) %>%
    dplyr::distinct() %>%
    # order_condition_name() %>%
    dplyr::arrange(condition_i) %>%
    numerify_regressors() %>%
    norm_regressors() %>% {
      if (obj$query == '6369') {
        .$tnf_conc <- 0
      }
      .
    } %>%
    add_binary_regressors(regressor_vars = reg_vars) %>%
    {
      if (all(c('duration', 'tnf_conc', 'ifn_conc') %in%
          colnames(.))) {
        . <- dplyr::mutate(., duration =
          if_else(ifn_conc == 0 & tnf_conc == 0, NA_real_, duration))
      }
      .
    } %>%
    dplyr::select(-any_of(c('duration_bin'))) %>%
    # debug_pipe() %>%
    # {
    #   idxs <-
    #     match(
    #       levels(SummarizedExperiment::colData(sce)$condition_name),
    #       as.character(.$condition_name)
    #     )
    #   .[idxs, ] } %>%
    { . }

  ## Condition names of the query experiment
  query_cn <-
    query_sa[, c('condition_i', 'condition_name')] %>%
    dplyr::distinct() %>%
    dplyr::arrange(condition_i) %>%
    pull(condition_name)

  ## (Unnormalized) Nhoods counts matrix, how do the query conditions
  ## distribute over the neighbourhoods?
  ## Dimensions: [query experimental conditions x Nhoods]
  if (F) {
    SummarizedExperiment::colData(sce)$condition_i <-
      factor(SummarizedExperiment::colData(sce)$condition_i,
        levels = 
          sort(unique(SummarizedExperiment::colData(sce)$condition_i)))
    sce <-
      countCells(
        sce,
        meta.data = data.frame(SummarizedExperiment::colData(sce)),
        samples = c('condition_i')
        # samples = c('stim_group')
      )
    UCM <- t(as.matrix(nhoodCounts(sce)))
    rownames(UCM) <- query_cn
    # nhoodCounts(sce)
    # UCM <- UCM[match(1:nrow(UCM), as.integer(rownames(UCM))), ]
    # dn <- dimnames(UCM)
  } else {
    UCM <-
      query_dtf %>%
      dplyr::select(matches('^CN\\d+$')) %>%
      as.matrix() %>%
      t()
  }

  if (!is.null(min_Nhood_size) && min_Nhood_size > 0) {
    valid_Nhoods <- which(colSums(UCM) >= min_Nhood_size)
    UCM <- UCM[, valid_Nhoods, drop = F]
  } else {
    valid_Nhoods <- 1:ncol(UCM)
  }
  obj$valid_Nhoods <- valid_Nhoods

  sel_Nhood_names <- colnames(nhoods(obj$sce))[valid_Nhoods]
  if (frac_included_cells(obj$sce, sel_Nhood_names) <= .025) {
    obj$agg_error <- NULL
    return(NULL)
  }

  ## Normalized counts matrix, scale the rows (each condition sums
  ## to 1)
  CM <- tryCatch(diag(1/rowSums(UCM)) %*% UCM,
    error = function(e) { print(e) })

  if (scale_by_max_aff) {
    max_aff_n <- obj$max_aff[valid_Nhoods] /
      sum(obj$max_aff[valid_Nhoods])
    CM <- CM %*% diag(max_aff_n)
    CM <- diag(1/rowSums(CM)) %*% CM
    stopifnot(rowSums(CM) - 1 < 1e-9)
  }

  ## Get distance matrix, which stores how close (in latent space)
  ## SC query Nhoods are to reference conditions
  ## Unmodified distances are maintained in obj$distM
  DM <- obj$distM
  dn <- dimnames(DM)

  ## Potentially modulate distance matrix
  if (ref_weights == 'uniqueness') {
    ref_factors <- compute_USF(ref_M)
    stopifnot(colnames(DM) == names(ref_factors))
    dn <- dimnames(DM)
    DM <- DM %*% diag(ref_factors)
    dimnames(DM) <- dn
  } else if (ref_weights == 'ref_reconstruction') {
    if (maartenutils::null_dat(ref_factors)) {
      rlang::warn('No ref_reconstruction found')
      return(NULL)
    }
    ref_factors <-
      tibble(condition_name = colnames(DM)) %>%
      dplyr::left_join(ref_factors, by = 'condition_name') %>%
      dplyr::mutate(recon_f = ifelse(is.na(recon_f), 1, recon_f))
    stopifnot(ref_factors$condition_name == colnames(DM))
    DM <- DM %*% diag(ref_factors$recon_f)
    dimnames(DM) <- dn
  }

  ## Compute Nhood to experimental ref affinity using (scaled)
  ## distances
  affM <- obj$compute_affM(distM = DM, ltb = ltb)

  ref_sa <-
    extract_sa(ref_dtf, meta_fields = reg_vars) %>%
    numerify_regressors() %>%
    add_binary_regressors(regressor_vars = reg_vars) %>%
    norm_regressors() %>%
    dplyr::select(-any_of(c('duration_bin'))) %>%
    as.matrix() %>%
    { . }
  if (F) {
    ## This should be redundant, as ref_sa and ref_M both derive from
    ## ref_dtf and their rows are not reordered
    match_col <- find_match_col(rownames(ref_M), dtf = ref_dtf)
    idxs <- match(colnames(WM), ref_dtf[[match_col]])
    ref_sa <- ref_sa[idxs, ]
  }

  if (!is.null(reg_vars)) {
    ## Penalty matrices, storing the penalty for each query
    ## experimental condition to each reference experimental
    ## condition
    penalty_Ms <-
      purrr::map(auto_name(reg_vars), function(rv) {
        if (is.null(query_sa[[rv]])) return(NULL)
        t(outer(ref_sa[, rv], query_sa[[rv]], '-'))
      }) %>%
      purrr::discard(is.null)
  }

  if (error_level_computation == 'conditions') {
    ## Weights matrix, how affine is each query condition to each
    ## ref condition? Marginalize the Nhoods out:
    ## [ref experimental conditions x Nhoods] times
    ## [Nhoods x query experimental conditions]
    WM <- CM %*% affM[valid_Nhoods, ]
    if (F) {
      apply(CM, 1, sum)
      apply(CM, 1, max)
      apply(CM, 2, sum)
      apply(affM[valid_Nhoods, ], 2, sum)
    }
    if (FALSE) {
      ltb_retained_Nhoods <-
        1-mean(apply(affM[valid_Nhoods, ], 1, sum) == 0)
    }
    rownames(WM) <- query_cn
    WM_1 <- (rowSums(WM) - 1) < 1e-2
    # if (any(!WM_1, na.rm = T) && interactive()) {
    # browser(expr = any(!WM_1, na.rm = T))
    if (!all(WM_1, na.rm = T) || any(!WM_1, na.rm = T)) { browser() }
    stopifnot(all(WM_1, na.rm = T))

    if (verbose) {
      setNames(1:nrow(WM), rownames(WM)) %>%
        map(~sort(WM[.x, ]) %>% { .[length(.)] }) %>%
        print()
    }

    ## Mean prediction matrix
    PM <-
      { WM %*% as.matrix(ref_sa) } %>%
      set_rownames(query_cn)

    ## Does the CI include zero?
    IZM <- matrix(nrow = nrow(WM), ncol = ncol(ref_sa))
    CI_l <- matrix(nrow = nrow(WM), ncol = ncol(ref_sa))
    CI_h <- matrix(nrow = nrow(WM), ncol = ncol(ref_sa))
    for (i in 1:nrow(IZM)) {
      for (j in 1:ncol(IZM)) {
        if (!any(WM[i, ] > 0, na.rm = T)) {
          CI_l[i, j] <- NA_real_
          CI_h[i, j] <- NA_real_
          IZM[i, j] <- NA
        } else {
          CI <- weighted_t_ci(ref_sa[, j], weights = WM[i, ])
          CI_l[i, j] <- CI[1]
          CI_h[i, j] <- CI[2]
          IZM[i, j] <- CI[1] <= 0 && CI[2] >= 0
        }
      }
    }
    dimnames(IZM) <- dimnames(PM)
    dimnames(CI_l) <- dimnames(PM)
    dimnames(CI_h) <- dimnames(PM)

    out <- list('CM' = CM, 'WM' = WM, 'PM' = PM, 'IZM' = IZM)

    if (!is.null(reg_vars)) {
      ## Observed error matrix, Hadamard multiplication of errors and
      ## weights
      EM <-
        purrr::map_dfc(penalty_Ms, function(pM) {
          { WM * abs(pM) } %>%
            rowSums(na.rm = F)
      }) %>%
      as.matrix() %>%
      set_rownames(query_cn)

      iz_vars <- intersect(reg_vars, colnames(IZM))
      IZ <- IZM[, reg_vars] == (!query_sa[, reg_vars])

      out <- c(out, list('EM' = EM, 'IZ' = IZ))
    }

    if (include_bs_pres) {
      ## Bootstrap resampled sample annotation for each row in AM
      ## (i.e. each neighbourhood)
      pred_sa <- resim_pmf(AM, sa = ref_sa, N_draws = colSums(UCM))

      ## Summarize the resampled sample annotation
      pred_sum <-
        pred_sa %>%
        # { .[1:3] } %>%
        imap_dfr(function(.x, i) {
          out <- purrr::map_dfr(1:ncol(.x), function(j) {
            tibble(
              i = rep(i, 3),
              rv = rep(colnames(.x)[j], 3),
              var = c('q1', 'q5', 'q9'),
              value = quantile(.x[, j], probs = c(.1, .5, .9), na.rm = T)
            )
          })
          return(out)
        })
        out[['pred_sum']] <- pred_sum
    }
  } else if (error_level_computation == 'neighbourhoods') {
    tUCM <- t(UCM[, valid_Nhoods]) %>%
      { diag(1/rowSums(.)) %*% .}
    stopifnot(abs(rowSums(tUCM) - 1) <= 1e-6)

    finisher <- function(x) {
      out <- 
        bind_cols(x) %>% 
        set_colnames(reg_vars) %>% 
        as.matrix()
      rownames(out) <- unlist(nhoodIndex(obj$sce))[valid_Nhoods]
      return(out)
    }

    PM <- purrr::map(reg_vars, function(rv) {
      # rv <- 'duration'
      ## Compute mean predicted score of each neighborhood
      ## (weighted predicted regression var)
      apply(affM[valid_Nhoods, ], 1, function(x) 
        weighted.mean(ref_sa[, rv], w = x))
    }) %>% finisher()

    MES <- purrr::map(reg_vars, function(rv) {
      ## Compute mean expected score MES for each neighborhood
      ## (weighed-mean of regression var)
      apply(tUCM, 1, function(x) 
        weighted.mean(query_sa[[rv]], w = x))
    }) %>% finisher()

    EM <- purrr::map(reg_vars, function(rv) {
      ## Compute the deviation between expected and predicted
      PM[, rv] - MES[, rv]
    }) %>% finisher()

    out <- list('PM' = PM, 'MES' = MES, 'EM' = EM)
  }

  # obj$agg_error <- out
  return(out)
}
