test_dim <- function(M) {
  return(!is.null(M) || any(dim(M) == 0))
}


gen_TRANSACT_wrapper <- function(reference_M, query_M) {
  tmp_dir <- stringr::str_replace(tempfile(), 'file', 'dir')
  dir.create(tmp_dir)

  force(reference_M)
  force(query_M)

  if (!test_dim(reference_M) || !test_dim(query_M))
    return(NULL)

  reference_M <- as.data.frame(reference_M)
  query_M <- as.data.frame(query_M)

  readr::write_csv(reference_M, file.path(tmp_dir, 'reference.csv'))
  readr::write_csv(query_M, file.path(tmp_dir, 'query.csv'))

  TRANSACT_caller <- function(
    TRANSACT_params = list(
      log_gamma = -5,
      N_source_components = 10L,
      N_target_components = 10L,
      N_PV = 5L, kernel = 'mallow'
    ),
    post_TRANSACT_params = list(
      center_ref = c(TRUE),
      center_query = c(FALSE)
    ), verbose = F,
    perform_validity_tests = T,
    ncores = 1L) {

    TRANSACT_params <- format_TRANSACT_params(TRANSACT_params)

    id_string <-
      c(TRANSACT_params, post_TRANSACT_params) %>%
      {
        if (.$kernel == 'linear') {
          . <- modifyList(., list(log_gamma = NA))
        }
        .
      } %>%
      modifyList(c(list(ps = '_'))) %>%
      with(paste(
        kernel,
        prepend_string(log_gamma, ps),
        prepend_string(N_source_components, ps),
        prepend_string(N_target_components, ps),
        prepend_string(N_PV, ps),
        prepend_string(center_ref, ps),
        prepend_string(center_query, ps),
        sep = '')
      )

    py_file <- file.path(Sys.getenv('python_dir'), 'run_transact.py')
    stopifnot(file.exists(py_file))
    command <- glue::glue('python {py_file} \\
      {tmp_dir} \\
      {ncores} \\
      {TRANSACT_params$kernel} \\
      {TRANSACT_params$log_gamma} \\
      {TRANSACT_params$N_source_components} \\
      {TRANSACT_params$N_target_components} \\
      {TRANSACT_params$N_PV} \\
      {post_TRANSACT_params$center_ref} \\
      {post_TRANSACT_params$center_query} \\
      {id_string}')

    system(command, wait = T)
    # print(list.files(tmp_dir))

    detected_types <-
      stringr::str_subset(list.files(tmp_dir), id_string) %>%
      stringr::str_replace(id_string, '') %>%
      stringr::str_replace('_\\.csv$', '')

    if (!all(c('reference_CF', 'query_CF') %in% detected_types)) {
      rlang::abort(glue::glue('reference_CF and/or query_CF not \\
          found in output, cache in {tmp_dir}'))
    }

    tryCatch({
      out <-
        detected_types %>%
        maartenutils::auto_name() %>%
        purrr::map(function(on) {
          readr::read_csv(
            file.path(tmp_dir, glue::glue('{on}_{id_string}.csv')),
            show_col_types = FALSE
          )
        })

      out$reference_CF <- as.matrix(out$reference_CF)
      out$query_CF <- as.matrix(out$query_CF)
      colnames(out$reference_CF) <- colnames(out$query_CF) <-
        glue::glue('CF{1:ncol(out$reference_CF)}')
      rownames(out$reference_CF) <- rownames(reference_M)
      rownames(out$query_CF) <- rownames(query_M)

      if (perform_validity_tests) {
        if ('CS_mat' %in% detected_types) {
          cs <- svd(out$CS_mat)$d
          out$valid_cs <- all(abs(cs) <= 1)
        }

        if ('EV_mat' %in% detected_types) {
          norms <- out$EV_mat %>%
            apply(2, function(x) norm(matrix(x, ncol = 1), type = 'F'))
          out$valid_norm <- all(norms <= 1.000001)
        }
      }

      if ('Kst_mat' %in% detected_types) {
        out$Kst_var <- var(as.vector(as.matrix(out$Kst_mat)))
      }
    }, error = function(e) { print(e); NULL })

    return(out)
  }

  return(TRANSACT_caller)
}


extract_transact_Ms <- function(
  reference_so,
  query_so,
  reference_M = NULL,
  query_M = NULL,
  normalize_sample_list_params = list(
    norm_method = 'logCPM',
    genes = NULL,
    genelist = NULL,
    GDR_thresh = NULL,
    Z_scale = FALSE
  )) {

  if (F) {
    ## Somehow doesn't work
    sample_list <-
      rlang::exec(
        normalize_sample_list,
        sample_list = list(
          'reference' = reference_so,
          'query' = query_so
        ),
        !!!normalize_sample_list_params
      )
  } else {
    sample_list <- normalize_sample_list(
      sample_list = list(
        'reference' = reference_so,
        'query' = query_so
      ),
      norm_method = normalize_sample_list_params$norm_method,
      genelist = normalize_sample_list_params$genelist,
      genes = normalize_sample_list_params$genes,
      GDR_thresh = normalize_sample_list_params$GDR_thresh,
      Z_scale = normalize_sample_list_params$Z_scale
    )
  }

  if (is.null(sample_list)) return(NULL)

  shared_genes <- find_shared_genes(sample_list)
  sample_list <- purrr::map(sample_list,
    ~as.data.frame(t(so2M(.x, datatype = 'data'))))
  out <- purrr::map(sample_list, ~.x[, shared_genes])

  stopifnot(names(out) == c('reference', 'query'))
  return(out)
}


#' Deprecated
#'
#'
read_Ms <- function(
  experiments,
  sc_mode = 'pseudobulk',
  genes = NULL,
  genelist = NULL,
  norm_method = 'none') {

  # experiments <-
  #   stringr::str_replace(experiments, '(\\d{4}).*', '\\1')

  Ms <-
    read_preproc_experiments(
      experiments = experiments,
      sc_mode = sc_mode
    ) %>%
    normalize_sample_list(
      norm_method = norm_method,
      merge_experiments = F,
      genes = genes,
      genelist = genelist,
    )
  # M <- purrr::exec(cbind, !!!map(Ms, so2M)) %>%
  #   t() %>%
  #   as.data.frame
  # map(Ms, dim)
  # map(Ms, colSums)
  # map(Ms, colSums)
  Ms <- purrr::discard(Ms, is.null)
  out <- Ms %>%
    purrr::map(~as.data.frame(t(so2M(.x, datatype = 'data'))))
  # map(out, max)
  return(out)
}


default_TRANSACT_params <-
  list(
    N_source_components = 10L,
    N_target_components = 10L,
    N_PV = 5L,
    log_gamma = -5L,
    kernel = 'mallow'
  )


format_TRANSACT_params <- function(TRANSACT_params) {
  TRANSACT_params <- modifyList(
    default_TRANSACT_params,
    TRANSACT_params
  )

  TRANSACT_params$N_PV <-
    as.integer(TRANSACT_params$N_PV)
  TRANSACT_params$N_source_components <-
    as.integer(TRANSACT_params$N_source_components)
  TRANSACT_params$N_target_components <-
    as.integer(TRANSACT_params$N_target_components)

  # if (TRANSACT_params$kernel == 'linear') {
  #   TRANSACT_params$log_gamma <- NULL
  # } else {
  #   TRANSACT_params$log_gamma <- round(TRANSACT_params$log_gamma, 3)
  # }

  return(TRANSACT_params)
}


default_norm_sample_list_args <- 
  list(
    genes = NULL,
    genelist = 'informativeV15',
    norm_method = 'CPM',
    Z_scale = FALSE
  )


compute_TRANSACT_embedding <- function(
  reference_so = NULL,
  query_so = NULL,
  reference_M = NULL,
  query_M = NULL,
  reference = extract_experiment(reference_so),
  query = extract_experiment(query_so),
  reference_sa = reference_so@meta.data,
  query_sa = query_so@meta.data,
  normalize_sample_list_params = default_norm_sample_list_args,
  TRANSACT_params = list(
    kernel = 'mallow',
    log_gamma = -5,
    N_PV = 5L,
    N_source_components = 15L,
    N_target_components = 10L
  ),
  post_TRANSACT_params = list(
    center_ref = c(TRUE),
    center_query = c(FALSE)
  ),
  add_cosine_info = T,
  ncores = 1L) {

  stopifnot(!is.null(reference))
  stopifnot(!is.null(query))

  stopifnot(!is.null(reference_sa))
  stopifnot(!is.null(query_sa))

  if (is.null(reference_M) || is.null(query_M)) {
    Ms <- extract_transact_Ms(
      reference_so = reference_so,
      query_so = query_so,
      normalize_sample_list_params = normalize_sample_list_params
    )
    if (is.null(reference_M)) {
      reference_M <- Ms[['reference']]
    }
    if (is.null(query_M)) {
      query_M <- Ms[['query']]
    }
    rm(Ms); gc()
  }
  scn <- intersect(colnames(reference_M), colnames(query_M))
  reference_M <- reference_M[, scn]
  query_M <- query_M[, scn]

  TRANSACT_caller <- gen_TRANSACT_wrapper(
    reference_M = reference_M,
    query_M = query_M
  )
  if (is.null(TRANSACT_caller)) return(NULL)

  TRANSACT_output <- TRANSACT_caller(
    TRANSACT_params = TRANSACT_params,
    post_TRANSACT_params = post_TRANSACT_params,
    ncores = ncores
  )
  if (is.null(TRANSACT_output)) return(NULL)

  ref_dtf <- merge_sa_on_rownames(
    M = TRANSACT_output$reference_CF,
    experiment = reference,
    sa = reference_sa
  )

  query_dtf <- merge_sa_on_rownames(
    M = TRANSACT_output$query_CF,
    experiment = query,
    sa = query_sa
  )

  dtf <- harmonize_bind_rows(ref_dtf, query_dtf)
  dtf <- order_duration(dtf)
  dtf <- order_stim_group(dtf)
  dtf <- order_condition_name(dtf)

  if (add_cosine_info) {
    cs <- svd(TRANSACT_output$CS_mat)$d
    cosine_dtf <- tibble(
      log_gamma = TRANSACT_params$log_gamma,
      i = seq_along(cs),
      cs = cs,
      Kst_var = TRANSACT_output$Kst_var,
      valid_norm = TRANSACT_output$valid_norm,
      valid_cs = TRANSACT_output$valid_cs
    )
    attr(dtf, 'cosine_dtf') <- cosine_dtf
  }

  return(dtf)
}


plot_gamma_tit_cosine_sim <- function(gamma_tit, return_plot = T) {
  if (maartenutils::null_dat(gamma_tit))
    return(NULL)

  if ('log_gamma_label' %in% colnames(gamma_tit)) {
    cv <- 'log_gamma_label'
  } else {
    cv <- 'log_gamma'
  }

  p1 <- gamma_tit %>%
    dplyr::mutate(i = ordered(i)) %>%
    dplyr::mutate(log_gamma = ordered(log_gamma)) %>%
    ggplot(aes_string(x = 'i', y = 'cs', colour = cv, 
        group = cv)) +
    geom_line() +
    # geom_point(data = dplyr::filter(gamma_tit, optimal_gamma == T)) +
    ylab('Cosine similarity') +
    scale_x_discrete(name = 'PV index', expand = c(0, 0)) +
    # scale_x_continuous(name = 'PV index', expand = c(0, 0)) +
    scale_colour_discrete(name = 'log10(gamma)') +
    guides(colour = guide_legend(ncol = 3))

  return(p1)
  # o_fn <- file.path(Sys.getenv('img_dir'),
  #   glue::glue('mallow_cs_{reference}_{query}_{genelist}\\
  #     {make_flag(N_source_components)}\\
  #     {make_flag(N_target_components)}.pdf'
  #   )
  # )
  # print_plot_eval({ print(p1) },
  #   width = 8.7, height = 10,
  #   filename = o_fn)
  # return(o_fn)
}


compute_TRANSACT_embedding_grid <- function(
  log_gamma,
  reference,
  query,
  reference_M,
  query_M,
  reference_sa,
  query_sa,
  TRANSACT_params = tibble::tibble(
    kernel = 'mallow',
    kernel_params = list('log_gamma' = -6),
    N_source_components = 10L,
    N_target_components = 10L,
    N_PV = pmin(5L, N_source_components, N_target_components)
  ) %>% as.list(),
  N_PV_umap_UL = 10L,
  N_PV_umap_LL = 2L,
  ncores = 1L,
  verbose = F) {

  log_gamma <- sort(force_numeric(log_gamma))
  
  # if (TRANSACT_params$kernel == 'cosine') {
  #   log_gamma <- log_gamma[1]
  # }

  param_grid <- tidyr::expand_grid(
    # center_ref = c(TRUE, FALSE),
    # center_query = c(TRUE, FALSE),
    tibble(
      center_ref = c(TRUE),
      center_query = c(FALSE)
    ),
    # tibble(
    #   center_ref = c(FALSE),
    #   center_query = c(FALSE)
    # ),
    # gamma_v = c(1e-5),
    log_gamma = log_gamma
  )
  ## N_mismatches
  ## N_combinations
  # exp(-gamma / N_combinations * N_mismatches)


  if (F) {
    param_grid <-
      param_grid %>%
      dplyr::mutate(dtf = purrr::pmap(
        list(center_ref, center_query, log_gamma),
        function(center_ref, center_query, log_gamma) {
          dtf <- tryCatch(compute_TRANSACT_embedding(
            reference = reference,
            query = query,
            reference_M = reference_M,
            query_M = query_M,
            reference_sa = reference_sa,
            query_sa = query_sa,
            ncores = ncores,
            TRANSACT_params = TRANSACT_params %>%
              modifyList(list('log_gamma' = log_gamma)),
            normalize_sample_list_params = list(
              genes = NULL,
              genelist = genelist
            ),
            post_TRANSACT_params = list(
              center_ref = center_ref,
              center_query = center_query
            )
          ), error = function(e) { print(e); NULL })
          return(dtf)
        })
      )
  } else {
    ## Early stopping; detect for which gamma the cosine similarity
    ## between PVs begins to be acceptable; stop iterating over the
    ## ascending-sorted gammas
    ## as soon as i) it has already been acceptable in a past
    ## iteration (a reasonable minimum value for gamma has
    ## already been surpassed) and ii) the cosine similarity is no
    ## longer acceptable.

    ## Init column list, don't know how else to do this at this point
    param_grid$dtf <- map(1:nrow(param_grid), ~.x)
    ## Hit the acceptable range, where cosine similarity of first PV
    ## is at least .1
    left_AR <- hit_AR <- F
    i <- 1
    while (i <= nrow(param_grid) && !left_AR) {
      if (verbose) {
        message(paste('i', i, 'left_AR', left_AR, 'hit_AR', hit_AR))
      }
      center_ref <- param_grid[[i, 'center_ref']]
      center_query <- param_grid[[i, 'center_query']]
      log_gamma <- param_grid[[i, 'log_gamma']]

      param_grid$dtf[[i]] <-
        tryCatch(compute_TRANSACT_embedding(
            reference = reference,
            query = query,
            reference_M = reference_M,
            query_M = query_M,
            reference_sa = reference_sa,
            query_sa = query_sa,
            ncores = ncores,
            TRANSACT_params = TRANSACT_params %>%
              modifyList(list('log_gamma' = log_gamma)),
            normalize_sample_list_params = list(
              genes = NULL,
              genelist = genelist
              ),
            post_TRANSACT_params = list(
              center_ref = center_ref,
              center_query = center_query
            )
        ), error = function(e) { print(e); NULL })

      CM <- attr(param_grid[i, ][['dtf']][[1]], 'cosine_dtf')
      ## Does this index hit the acceptable range?
      ci_hit_AR <- 
        CM %>%
        dplyr::filter(.data[['i']] == 1) %>%
        dplyr::pull(cs) %>%
        { any(. >= .1, na.rm = T) }
      if (ci_hit_AR && !hit_AR) {
        hit_AR <- T
      } else if (!ci_hit_AR && hit_AR) {
        left_AR <- T
      }
      i <- i + 1
    }

    if (!hit_AR) {
      ## Not a single gamma worked out
      return(NULL)
    } else if (left_AR && i < nrow(param_grid)) {
      ## At least one valid gamma was found and higher gammas result
      ## in invalid results
      param_grid <- param_grid[1:(max(1, i-1)), ]
    }
  }

  ## Add UMAP coordinates to each dtf

  ## Upper limit (UL) of the number of CFs to use
  N_PV_umap_UL <- with(TRANSACT_params,
    min(N_source_components, N_target_components, N_PV, N_PV_umap_UL))
  if (N_PV_umap_UL < 2)
    return(param_grid)

  ## Lower limit (LL) of the number of CFs to use
  if (F) {
    N_PV_umap_LL <- with(TRANSACT_params,
      min(N_source_components, N_target_components, 
        N_PV, N_PV_umap_LL))
  } else {
    ## Only compute UMAP for ALL CFs, I never use the ones with fewer
    ## CFs anyway
    N_PV_umap_LL <- N_PV_umap_UL
  }

  param_grid <- 
    param_grid %>%
    tidyr::expand_grid(N_PV_umap = N_PV_umap_LL:N_PV_umap_UL) %>%
    dplyr::mutate(u_dtf = purrr::pmap(list(dtf, N_PV_umap),
        function(dtf, N_PV_umap) {
      if (maartenutils::null_dat(dtf)) return(NULL)
      cns <- as.character(glue::glue('CF{1:N_PV_umap}')) %>%
        intersect(colnames(dtf))
      if (!length(cns)) {
        return(dtf)
      }
      all_NA <- apply(dtf[, cns], 2, function(x) all(is.na(x)))
      cns <- cns[!all_NA]
      if (length(cns) < 2) {
        return(dtf)
      }
      dtf <- dtf %>%
        add_umap(
          column_selector = any_of(cns),
          add_source_type = F,
          verbose = T
        )
      return(dtf)
    }))

  ## Pre-extract the cosine similarity matrices for easier viewing
  param_grid$PV_cosine_sim <- map(1:nrow(param_grid), function(i) {
      attr(param_grid[i, ][['dtf']][[1]], 'cosine_dtf')
  })

  return(param_grid)
}


#' DEPRECATED
#'
#'
assess_reference_reconstruction_dep <- function(
  log_gamma, reference_M, query_M,
  reference_sa, N_source_components, N_target_components) {

  stopifnot(is.data.frame(reference_sa))
  stopifnot(is.data.frame(reference_M))

  TRANSACT_caller <- gen_TRANSACT_wrapper(
    reference_M = reference_M,
    query_M = query_M
  )
  if (is.null(TRANSACT_caller)) return(NULL)
  max_N_PV <- min(N_source_components, N_target_components)
  # max_N_PV <- 5
  TRANSACT_output <- TRANSACT_caller(
    TRANSACT_params = list(
      kernel = 'mallow',
      log_gamma = log_gamma,
      N_PV = max_N_PV,
      N_source_components = N_source_components,
      N_target_components = N_target_components
    ),
    post_TRANSACT_params = list(
      'center_ref' = FALSE,
      'center_query' = FALSE
    )
  )

  if (T) {
    k_fn <-
      list.files(
        environment(TRANSACT_caller)$tmp_dir,
        full.names = T
      ) %>%
      stringr::str_subset('Ks_mat')
  } else {
    k_fn <- '/tmp/RtmpiacG4m/dir3c6fc678fc8b80/Ks_mat_mallow-5_10_10_5_TRUE_FALSE.csv'
  }
  # readr::read_csv(k_fn)

  py_file <- file.path(Sys.getenv('python_dir'), 'kernel_pca.py')
  reticulate::source_python(py_file)
  ref_kPCA <-
    run_kPCA(
      kernel_file = k_fn,
      n_jobs = 1L,
      # gamma = 10^doc_par$NH_optimal_log_gamma,
      gamma = 10^log_gamma,
      N_NLPC = max_N_PV
    )
  ref_NLPCs <- ref_kPCA[['eigenvectors_']] %*%
    sqrt(diag(ref_kPCA[['eigenvalues_']]))

  U_abs <- abs(svd(t(TRANSACT_output$reference_CF) %*% ref_NLPCs)$u)
  N_dim <- dim(U_abs)[1]
  PV_mapping <- rep(0, N_dim)
  for (i in 1:dim(U_abs)[1]) {
    remaining_samples <- setdiff(1:N_dim, PV_mapping)
    idx <- which.max(U_abs[i, remaining_samples])
    # idx <- which(ord == (N_dim-i+1))
    # while (!idx %in% PV_mapping) {
    # }
    PV_mapping[i] <- remaining_samples[idx]
    print(PV_mapping)
  }
  stopifnot(sort(unique(PV_mapping)) == 1:N_dim)

  ref_names <- rownames(reference_M)
  elongate_scores <- function(M, early_var_idx = floor(max_N_PV/2),
    permute_PV = NULL) {
    out <- M^2 %>%
      as.data.frame() %>%
      tidyr::pivot_longer(everything()) %>%
      dplyr::mutate(sample_name = rep(ref_names, each = max_N_PV)) %>%
      dplyr::mutate(sample_name = stringr::str_replace(sample_name,
          '^\\d{4}_', '')) %>%
      dplyr::mutate(sample_name = stringr::str_replace(sample_name,
          '_0_ng_ml_tnfa', '')) %>%
      dplyr::mutate(name = stringr::str_replace(name, '^[A-Za-z]+', '')) %>%
      dplyr::rename(PV = name) %>%
      dplyr::mutate(PV = as.integer(PV)) %>%
      ## Ensure PV is 1-based and not 0-based
      dplyr::mutate(PV = PV - min(PV) + 1)

    if (!is.null(permute_PV)) {
      out <- out %>%
        dplyr::mutate(PV = recode(PV, !!!permute_PV)) %>%
        dplyr::arrange(sample_name, PV)
    }

    out %>%
      dplyr::group_by(sample_name) %>%
      dplyr::mutate(cum_score = cumsum(value)) %>%
      dplyr::summarize(across(),
        early_var = cum_score[early_var_idx] / cum_score[max_N_PV])
  }

  p_dat <- elongate_scores(TRANSACT_output$PV_projected_source)
  p_dat_source <-
    elongate_scores(ref_NLPCs, 2, permute_PV = PV_mapping) %>%
    rename_with(~paste0(.x, '_ref'))
  p_dat <- cbind(p_dat, p_dat_source) %>%
    dplyr::mutate(early_var_norm = early_var / early_var_ref)

  # so <- p_dat %>%
  #   dplyr::filter(PV == 1) %>%
  #   dplyr::arrange(value) %>%
  #   dplyr::pull(sample_name)
  so <- p_dat %>%
    dplyr::arrange(early_var) %>%
    dplyr::pull(sample_name) %>%
    unique()
  p_dat$sample_name <- factor(p_dat$sample_name, levels = so)

  sa <- reference_sa[, 'condition_name', drop=F] %>%
    tibble::rownames_to_column('sample_name')
  p_dat <- dplyr::left_join(p_dat, sa, by = 'sample_name')
  so <- p_dat %>%
    dplyr::arrange(early_var) %>%
    dplyr::pull(condition_name) %>%
    unique()
  p_dat$condition_name <- factor(p_dat$condition_name, levels = so)

  return(p_dat)
}


assess_reference_reconstruction <- function(
  log_gamma, reference_M, query_M,
  reference_sa, 
  N_source_components, N_target_components, N_PV,
  kernel = 'mallow', ncores = 1L) {

  stopifnot(is.data.frame(reference_sa))
  stopifnot(is.data.frame(reference_M))

  TRANSACT_caller <- gen_TRANSACT_wrapper(
    reference_M = reference_M,
    query_M = query_M
  )
  if (is.null(TRANSACT_caller)) return(NULL)
  TRANSACT_output <- tryCatch({
    TRANSACT_caller(
      TRANSACT_params = list(
        kernel = kernel,
        log_gamma = log_gamma,
        N_PV = N_PV,
        N_source_components = N_source_components,
        N_target_components = N_target_components
      ),
      ncores = ncores,
      post_TRANSACT_params = list(
        'center_ref' = TRUE,
        'center_query' = FALSE
      )
    )
  }, error = function(e) { print(e); NULL }) 
  if (is.null(TRANSACT_output) || is.null(TRANSACT_output$reference_CF)) return(NULL)

  ## Scaling the vectors by the eigenvalues makes the total sample
  ## composition relevant; factors describing a lot of variance (i.e.
  ## samples) will become more heavily weighted in downstream analyses
  ref_CF <- TRANSACT_output$reference_CF
  # ref_CF <- 
  #   ref_CF %*% diag(1/apply(ref_CF, 2, function(x) norm(as.matrix(x), type = 'F')))
  # apply(ref_CF, 2, function(x) norm(as.matrix(x), type = 'F'))

  k_fn <-
    list.files(
      environment(TRANSACT_caller)$tmp_dir,
      full.names = T
    ) %>%
    stringr::str_subset('Ks_mat')
  
  py_file <- file.path(Sys.getenv('python_dir'), 'kernel_pca.py')
  reticulate::source_python(py_file)
  ref_kPCA <-
    run_kPCA(
      kernel_file = k_fn,
      n_jobs = 1L,
      # gamma = 10^doc_par$NH_optimal_log_gamma,
      gamma = 10^log_gamma,
      N_NLPC = N_source_components
    )

  ref_NLPCs <- 
    ref_kPCA[['eigenvectors_']] %*% sqrt(diag(ref_kPCA[['eigenvalues_']]))
  # ref_NLPCs <- ref_kPCA[['eigenvectors_']] 
  # apply(ref_NLPCs, 2, function(x) 
    # norm(as.matrix(x), type = 'F'))

  stopifnot(rownames(reference_M) == rownames(reference_sa))
  rownames(reference_M) <- reference_sa$condition_name
  rownames(ref_NLPCs) <- reference_sa$condition_name
  rownames(ref_CF) <- reference_sa$condition_name
  ref_names <- reference_sa$condition_name

  unstim_samples <- stringr::str_subset(ref_names, '(U|u)nstimulated')
  too_low_ifny <- stringr::str_subset(ref_names, '0.01 ng/ml IFNy')

  primary_samples <- 
    reference_sa %>%
    dplyr::filter(!(ifn_conc == 0.01)) %>%
    dplyr::filter(tnf_rank == 1 | ifn_rank == 1) %>%
    pull(condition_name) %>%
    as.character() %>%
    naturalsort::naturalsort()

  secondary_samples <- 
    ref_names %>%
    setdiff(primary_samples) %>%
    setdiff(too_low_ifny) %>%
    naturalsort::naturalsort()

  compute_scores <- function(X) {
    ## X has unit 'proximity' or 'affinity'
    ## TODO revert this back to distance, which is easier to
    ## conceptually
    X <- 1/as.matrix(dist(X))
    dimnames(X) <- list(ref_names, ref_names)
    X <- X[!rownames(X) %in% too_low_ifny, !colnames(X) %in% too_low_ifny]
    ## The highest concentration of IFNy is assumed to be a good
    ## landmark; it should always substantially differ from the
    ## unstimulated sample
    ref_value <- mean(X['100 ng/ml IFNy - 24h', unstim_samples])
    X <- X[
      which(!rownames(X) %in% unstim_samples), 
      which(colnames(X) %in% unstim_samples)]
    X <- X / ref_value
    return(X)
  }
  D_NLPC <- compute_scores(ref_NLPCs)
  D_CF <- compute_scores(ref_CF)

  ## The higher this ratio, the more 'compressed' the sample
  ## is in CF space as compared to direct/reference only NLPC space. 
  D_R <- D_CF / D_NLPC
  score_trans <- identity
  D_R <- score_trans(D_R)
  if (is.matrix(D_R)) {
    scores <- apply(D_R, 1, mean, na.rm = T)
  } else {
    scores <- D_R
  }

  if (F) {
    sort(scores[primary_samples])
    sort(scores[secondary_samples])
    summary(scores[primary_samples])
    summary(scores[secondary_samples])
  }

  tibble::enframe(scores, 'condition_name', 'recon_f') %>%
    dplyr::right_join(tibble::enframe(rowMeans(D_NLPC),
        'condition_name', 'NLPC')) %>%
    dplyr::right_join(tibble::enframe(rowMeans(D_CF), 
        'condition_name', 'CF'))
}


#' Deprecated as long as targets and reticulate play nice together
#'
#'
kPCA_caller <- function(
  kernel_file = k_fn,
  ncores = 1L,
  gamma = 10^-6,
  N_NLPC = 10L) {

  py_file <- file.path(Sys.getenv('python_dir'), 'kernel_pca.py')
  stopifnot(file.exists(py_file))
  command <- glue::glue('python {py_file} \\
    {tmp_dir} \\
    {ncores} \\
    {N_NLPC}')
  system(command, wait = T)
}


expected_winners_by_exp <- list(
  '6600' = c('10 - 0 - 6', '10 - 100 - 6', '0 - 0 - 6'),
  '6369' = c(
    '0 - 1 - 2', '0 - 10 - 2', '0 - 100 - 2',
    '0 - 1 - 6', '0 - 10 - 6', '0 - 100 - 6',
    '0 - 1 - 12', '0 - 10 - 12', '0 - 100 - 12',
    '0 - 1 - 24', '0 - 10 - 24', '0 - 100 - 24'
  )
)

annotate_expected_winner <- function(dtf, experiment) {
  if (experiment %in% names(expected_winners_by_exp)) {
    dtf <- dtf %>% dplyr::mutate(
      expected_winner = sample_name %in%
        expected_winners_by_exp[experiment])
  } else {
    dtf <- dtf %>% dplyr::mutate(expected_winner = FALSE)
  }
  return(dtf)
}


annotate_ref_recon <- function(p_dat) {
  levels(p_dat$condition_name) <-
    levels(p_dat$condition_name) %>%
    stringr::str_replace('(?<=IFNy)( )(\\d)', '\n\\2') %>%
    stringr::str_replace(' - ', '\n') %>%
    stringr::str_replace('Unstimulated in vitro\n', '') %>%
    { . }

  ranking <-
    p_dat %>%
    dplyr::group_by(condition_name) %>%
    dplyr::summarize(
      C = cor(cum_score, cum_score_ref),
      C_norm = cor(
        cum_score/max(cum_score),
        cum_score_ref/max(cum_score_ref)
      ),
      S_norm = sum(abs(cum_score - cum_score_ref)[1:2]),
      # T_norm =
      #   mean(abs(cum_score/max(cum_score) -
      #     cum_score_ref/max(cum_score_ref))),
      # AUC = pracma::trapz(1:n(), cum_score),
      # AUC_ref = pracma::trapz(1:n(), cum_score_ref),
      # AUC_diff = AUC - AUC_ref,
      across(matches('expected_winner'), ~.x[1])
      ) %>%
    dplyr::arrange(S_norm)

  so <- ranking %>%
    dplyr::pull(condition_name)

  p_dat <- p_dat %>%
    dplyr::select(condition_name, cum_score, early_var,
      PV = PV_ref, cum_score_ref, expected_winner) %>%
    tidyr::pivot_longer(
      names_prefix = 'cum_score',
      cols = c(cum_score, cum_score_ref)
      # names_from = c(condition_name, PV),
    ) %>%
    dplyr::mutate(condition_name = factor(condition_name, levels = so)) %>%
    dplyr::rename(type = name) %>%
    dplyr::mutate(type = ifelse(
        type == '', '(Shared) PV', 'Ref-only kPCA'))

  return(p_dat)
}


gen_NH_report_filename <- function(reference, query, genelist, GDR,
  Z_scale, k, d, min_neighbourhood_size, N_source_components,
  N_target_components, N_PV, ...) {
  paste0('NH_TRANSACT',
          '-reference=', reference,
          '-query=', query,
          '-genelist=', genelist,
          '-GDR=',
          stringr::str_replace(
            as.character(round(GDR, 2)), '0\\.', ''),
          '-Z_scale=', Z_scale,
          '-k=', k,
          '-d=', d,
          '-min_neighbourhood_size=', min_neighbourhood_size,
          '-N_source_components=', N_source_components,
          '-N_target_components=', N_target_components,
          '-N_PV=', N_PV,
          # '-min_early_var=',
          # stringr::str_replace(
          #   as.character(round(min_early_var, 2)), '0\\.', ''),
          # '-D_gamma=', D_gamma,
          # '.html'
          # '.pdf'
          '')
}


extract_ref_dtf <- function(dtf) {
  ref_dtf <-
    dtf %>%
    dplyr::filter(sample_type == 'bulk') %>%
    order_condition_name() %>%
    dplyr::arrange(condition_name) %>%
    { . }
  return(ref_dtf)
}


extract_query_dtf <- function(dtf) {
  query_dtf <- dtf %>%
    dplyr::filter(is.na(sample_type) | sample_type != 'bulk') %>%
    { . }
  return(query_dtf)
}


extract_CF_M <- function(dtf) {
  if (!is.null(dtf$condition_name) &&
      !any(is.na(dtf$condition_name))) {
    rn <- dtf$condition_name
  } else if (!is.null(dtf$rn) && !any(is.na(dtf$rn))) {
    rn <- dtf$rn
  }
  ref_M <- dtf %>%
    dplyr::select(matches('CF|^NLPC_')) %>%
    as.matrix() %>%
    set_rownames(rn)
  return(ref_M)
}


transact_NH_object_loading <- rlang::expr({
  # stopifnot(!is.null(get_obj('ref_so')))

  ## LOAD OBJECTS
  ref_experiment <- unique(get_obj('ref_so')@meta.data$experiment)
  # get_obj('Nhood_stats')
  # get_obj('NH_TRANSACT_reference_reconstruction')
  # tar_read(TRANSACT_reference_reconstruction_comb_6369_informativeV15_15_10)

  stopifnot(!all(sapply(get_obj('NH_gamma_titration_embedding')$u_dtf,
        is.null)))

  dtf <- tryCatch({
    get_obj('NH_gamma_titration_embedding') %>%
      dplyr::mutate(log_gamma = factor_to_numeric(log_gamma)) %>%
      dplyr::filter(log_gamma == get_obj('NH_optimal_log_gamma')) %>%
      dplyr::filter(N_PV_umap == get_obj('N_PV')) %>%
      purrr::pluck('u_dtf', 1) %>%
      dplyr::mutate(tnf_conc = factor(tnf_conc)) %>%
      {
        ucn <- setdiff(colnames(get_obj('Nhood_stats')), colnames(.))
        dplyr::left_join(., get_obj('Nhood_stats')[, ucn],
          by = c('rn' = 'nhoodIndex'))
      } %>%
      add_umap(any_of(paste0('CF', 1:get_obj('N_PV')))) %>%
      dplyr::mutate(across(matches('^CN\\d+$'),
          function(x) ifelse(is.na(x), 0, x))) %>%
      dplyr::filter(!(experiment == get_obj('query') & N == 0)) %>%
      { . }
    }, error = function(e) { NULL })
  stopifnot(!is.null(dtf))

  coord_ranges <- gen_coord_ranges_fun(dtf)
  transact_id <- get_obj('transact_id')

  library(miloR)
  sce <- get_obj('agg_neighbourhoods')
  Nhood_names <- unlist(nhoodIndex(sce))
  idxs <- suppressWarnings(match(Nhood_names, as.numeric(dtf$rn)))
  query_UMAP_coords <-
    dtf[idxs, ] %>%
    dplyr::select(matches('UMAP')) %>%
    as.matrix() %>%
    set_rownames(Nhood_names)

  ref_dtf <- extract_ref_dtf(dtf)
  ref_M <- extract_CF_M(ref_dtf)
  query_dtf <- extract_query_dtf(dtf)
  query_M <- extract_CF_M(query_dtf)

  ref_CF_clustering <- gen_clust_object(t(ref_M), dist_f = 'euclidean')
  ref_CF_cut <- cutree(
    ref_CF_clustering,
    h = quantile(as.hclust(ref_CF_clustering)$height, .75)
  ) %>% sort()

  tryCatch({
    if (!test_rendering() || T) {
      # stopifnot(is.list(doc_par))
      stopifnot(is.integer(get_obj('N_PV')))
      # stopifnot(is.numeric(doc_par$min_early_var))
      # stopifnot(is.numeric(doc_par$D_gamma))
      stopifnot(is.numeric(get_obj('NH_optimal_log_gamma')))
      stopifnot(get_obj('NH_optimal_log_gamma') < 0)
      stopifnot(!maartenutils::null_dat(head(get_obj('NH_so')@meta.data)))
      stopifnot(get_obj('NH_so')@meta.data$experiment == get_obj('query'))
      emb_tit <- get_obj('NH_gamma_titration_embedding')
      if (maartenutils::null_dat(emb_tit)) {
        print(emb_tit)
      }
      if (all(sapply(emb_tit$dtf, is.null)))
        rlang::abort('All embeddings are NULL')
      stopifnot(any(dtf$experiment == get_obj('query')))
      dtf %>%
        dplyr::filter(experiment != ref_experiment) %>%
        dplyr::select(matches('CN\\d+')) %>%
        is.na() %>%
        { all(!.) } %>%
        stopifnot()
    }
  }, error = function(e) { print(e) })
})


#' Compute 'uniqueness' scaling factor
#'
#' This function aims to compute scaling factors to downweight objects
#' (rows) that are very similar to all other objects in the set
#' (matrix). Approach: compute Euclidean distance (dissimilarity)
#' between a bunch of objects (samples) and then compute the fold
#' difference between each object's mean distance to all other objects
#' and the global mean distance.
compute_USF <- function(M) {
  if (F) {
    corM <- cor(t(M))
    diag(corM) <- 0
    ## Factor vector
    FV <- apply(corM, 1, mean) / mean(corM[lower.tri(corM)])
  } else {
    distM <- as.matrix(dist(M))
    diag(distM) <- NA_real_
    FV <- apply(distM, 1, mean, na.rm = T) /
      mean(distM[lower.tri(distM)], na.rm = T)
  }
  return(FV)
}


extract_dtf_from_gamma_tit <- function(
  gamma_tit,
  lg = gamma_tit$log_gamma[1], 
  ref_experiment = 'comb',
  N_PV_umap = max(gamma_tit$N_PV_umap)) {

  dtf <- tryCatch({
    gamma_tit %>%
      dplyr::mutate(log_gamma = factor_to_numeric(log_gamma)) %>%
      dplyr::filter(log_gamma == lg) %>%
      dplyr::filter(N_PV_umap == .env[['N_PV_umap']]) %>%
      purrr::pluck('u_dtf', 1) %>%
      # add_umap(any_of(paste0('CF', 1:N_PV))) %>%
      dplyr::mutate(across(
          matches('^CN\\d+$'),
          function(x) ifelse(is.na(x), 0, x)
          )) %>%
      # dplyr::filter(!(experiment == query & N == 0)) %>%
      { . }
  }, error = function(e) { print(e); NULL })
  dtf$sample_type[dtf$experiment == ref_experiment] <- 'bulk'
  return(dtf)
}


extract_cosine_from_gamma_tit <- function(gamma_tit) {
  out <-
    gamma_tit %>%
    distinct(center_ref, center_query, log_gamma, .keep_all = T) %>%
    dplyr::select(
      -any_of(c('dtf', 'N_PV_umap', 'u_dtf', 'log_gamma'))) %>%
    tidyr::unnest(PV_cosine_sim)
  return(out)
}


get_i_feats <- function(genelist, q_sc_so, 
  include_variable_feats = 'union') {
  feats <- read_geneset(genelist)
  if (include_variable_feats == 'union') {
    feats <- union(feats, VariableFeatures(q_sc_so))
  } else if (include_variable_feats == 'intersect') {
    feats <- intersect(feats, VariableFeatures(q_sc_so))
  } else {
  }
  return(feats)
}
