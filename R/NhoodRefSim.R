NhoodRefSim <- R6::R6Class('NhoodRefSim',
  lock_objects = FALSE,
  lock_class = FALSE,
  public = list(
  print = function() {
    pf <- function(x) cat(x, '\n')
    # pf(self$ref_experiment)
  },

  #' Init an Nhood to reference similarity object
  #'
  #' @param dtf Dataframe
  #' @param sce
  #' @param query Name of the query experiment
  #' @param ref_sa Meta data of the reference experiment
  #' @param ref_experiment Name of the reference experiment
  #' samples by their uniqueness. Will be deduced from @ref_sa if
  #' absent.
  #' @param exclusion_affinity_bms
  #' @param primary_ref_samples
  #' @param ref_factors Reference-specific factors to weight the query
  #' to reference distances by
  #' @param N_PV
  #' @param ltb
  #' @param ref_CF_clustering
  #'
  #' @export
  #'
  #' @returns an NhoodRefSim object
  #'
  #' @examples
  #' NhoodRefSim$initialize()
  #'
  initialize = function(
    dtf,
    sce,
    query,
    ref_sa,
    ref_experiment = ref_sa$experiment[1],
    primary_ref_samples = list(),
    exclusion_affinity_bms = list(),
    ref_factors = NULL,
    N_PV = NULL,
    ## 'Label transfer bandwidth', parameter of the exponential
    ## label transfer similarity function
    ltb = NULL,
    ref_CF_clustering = NULL,
    ...) {

    library(miloR)
    library(SingleCellExperiment)

    stopifnot(!null_dat(dtf))
    stopifnot(!maartenutils::null_dat(sce))
    stopifnot(!maartenutils::null_dat(ref_sa))
    stopifnot(is.null(primary_ref_samples) ||
      is.list(primary_ref_samples))

    self$ref_factors <- ref_factors
    # stopifnot(!null_dat(self$ref_factors))
    self$dtf <- dtf
    self$sce <- sce
    self$Nhood_size <- 
      rowSums(as.matrix(miloR::nhoodCounts(self$sce)))
    names(self$Nhood_size) <- unlist(miloR::nhoodIndex(self$sce))
    self$ref_experiment <- ref_experiment
    self$query <- query
    self$ref_sa <- ref_sa %>%
      order_duration() %>%
      order_concentration() %>%
      arrange(across(
          any_of(c('condition_name', 'stim_group', 'duration'))))
    self$ltb <- ltb
    self$ref_CF_clustering <- ref_CF_clustering
    self$dtf <- order_duration(self$dtf)
    self$query_M <- self$dtf %>%
      extract_query_dtf() %>%
      extract_CF_M()
    self$ref_M <-
      dtf %>%
      extract_ref_dtf() %>%
      extract_CF_M()
    self$N_PV <- N_PV %||% ncol(self$query_M)

    if (!is.null(primary_ref_samples)) {
      if (!is.null(exclusion_affinity_bms)) {
        ## Learn the max tolerated similarity between reference
        ## samples from a set of pre-specified combinations that are
        ## expected to be indistinguishable
        ## A sane (and my...) default would be a bunch of unstimulated
        ## samples, observed at different durations
        aff_M <- 1/as.matrix(dist(self$ref_M))
        diag(aff_M) <- NA
        # summary(as.vector(aff_M))
        self$max_tol_aff <-
          min(purrr::map_dbl(exclusion_affinity_bms, function(x) {
            aff_M[x[1], x[2]]
          }))
      } else {
        ## Learn the max tolerated similarity between reference
        ## samples from the max similarities from query SC Nhoods to
        ## reference samples; take the Xth percentile of this
        ## distribution
        self$max_tol_aff <-
          { 1/matrix_euclidean(t(self$ref_M), t(self$query_M)) } %>%
          { apply(., 2, quantile, probs = 1) } %>%
          # { quantile(., probs = seq(0, 1, by = .1)) }
          quantile(probs = .9)
      }

      ## primary_ref_samples might be a list of vectors, if so,
      ## iterate over this list
      if (is.list(primary_ref_samples)) {
        for (s in primary_ref_samples) {
          output <- exclude_ref_sample_cor(
            M = self$ref_M,
            sa = self$ref_sa,
            sample_names = s,
            max_tol_aff = self$max_tol_aff
          )
          self$ref_M <- output$M
          self$ref_sa <- output$sa
        }
      } else {
        stop('Unexpected type for primary_ref_samples, ', '
          expecting a list')
      }
    }

    self$allowed_ref_samples <- rownames(self$ref_M)
    match_col <- find_match_col(
      self$allowed_ref_samples, dtf = self$dtf)
    self$dtf <-
      self$dtf %>%
      dplyr::filter(
        experiment %in% self$query |
        .data[[match_col]] %in% self$allowed_ref_samples
      )

    self$compute_distM()
  },

  #' Subset the reference sample annotation with a specific set of
  #' samples
  #'
  #'
  subset_ref_sa = function(sample_names = colnames(self$distM)) {
    match_col <- find_match_col(sample_names, dtf = self$ref_sa)
    stopifnot(!is.null(self$ref_sa[[match_col]]))
    idxs <- match(sample_names, self$ref_sa[[match_col]],
      nomatch = NULL) %>% setdiff(NA)

    out <-
      self$ref_sa %>%
      { .[idxs, ] } %>%
      # dplyr::select(-any_of(match_col)) %>%
      # tibble::rownames_to_column('sample_name') %>%
      # dplyr::select(sample_name, stim_group, duration) %>%
      dplyr::select(
        any_of(c(match_col, 'sample_name', 'stim_group')),
        matches('conc|duration')
      ) %>%
      # {
      #   tmp <- tibble(temp = sample_names) %>%
      #     set_colnames(match_col)
      #   dplyr::left_join(tmp, ., by = match_col)
      # } %>%
      # dplyr::select(-sample_name) %>%
      # order_stim_group() %>%
      order_concentration() %>%
      order_duration() %>%
      dplyr::select(-any_of(c(match_col, 'sample_name'))) %>%
      { . }

    return(out)
  },

  compute_distM = function() {
    self$distM <- matrix_euclidean(t(self$query_M), t(self$ref_M))
  },

  compute_affM = function(
    distM = self$distM,
    ltb = self$ltb,
    min_sum_prob = .0) {

    if (!test_dim(distM)) return(NULL)

    affM <- exp(-ltb * distM)
    self$max_aff <- apply(affM, 1, max)
    self$min_aff <- apply(affM, 1, min)
    self$sum_aff <- apply(affM, 1, sum)

    affM <-
      apply(affM, 1, function(x) x / sum(x, na.rm = T)) %>%
      t()
    affM[!is.finite(affM)] <- 0

    if (min_sum_prob > 0) {
      ref_idxs <- which(apply(self$affM, 2, sum) >= min_sum_prob)
      affM <- affM[, ref_idxs]
    }

    ## Make sure ordering is harmonious at this point, this assumes
    ## ref_sa is logically ordered
    match_col <- find_match_col(
      colnames(affM),
      dtf = self$ref_sa
    )
    stopifnot(!is.null(self$ref_sa[[match_col]]))
    idxs <- match(self$ref_sa[[match_col]], colnames(affM),
      nomatch = NULL) %>% setdiff(NA)
    affM <- affM[, idxs]

    if (maartenutils::null_dat(affM) || !test_dim(affM)) {
      rlang::warn('NULL affM')
    }
    self$affM <- affM
    return(affM)
  },

    compute_Nhood_stats = function(min_Nhood_size = 0) {
    ## Column, rows
    nearest_ref_idxs <-
      cbind(apply(self$distM, 1, which.min), 1:nrow(self$distM))

    NN_probs <-
      factor(as.vector(outer(c(1, 2.5, 5, 7.5), 10^(-1:6)))) %>%
      purrr::map_dfr(function(ltb) {
        ltb_N <- factor_to_numeric(ltb)
        ## Compute affinities
        affM <-
          exp(-ltb_N * self$distM) %>%
          apply(1, function(x) x / sum(x, na.rm = T))
        tibble(
          ltb = ltb,
          NN_prob = affM[nearest_ref_idxs],
          Nhood_idx = 1:nrow(self$distM),
          Nhood_name = rownames(self$distM)[Nhood_idx],
          Nhood_size = self$Nhood_size[Nhood_idx]
        )
      })

    gamma_stats <-
      NN_probs %>%
      group_by(ltb) %>%
      dplyr::summarize(
        med_NN_prob = median(NN_prob, na.rm = T),
        gamma_Nhoods_retained = sum(NN_prob > 0, na.rm = T),
        frac_gamma_Nhoods_retained = gamma_Nhoods_retained / n(),
        frac_gamma_cells_retained = frac_included_cells(self$sce,
          allowed_Nhoods = which(NN_prob > 0)),
        size_Nhoods_retained = sum(Nhood_size > min_Nhood_size, na.rm = T),
        frac_size_Nhoods_retained = size_Nhoods_retained / n(),
        frac_size_cells_retained = frac_included_cells(self$sce,
          allowed_Nhoods = which(Nhood_size > min_Nhood_size)),
        Nhoods_retained = sum(
          NN_prob > 0 & Nhood_size > min_Nhood_size, na.rm = T),
        frac_Nhoods_retained = Nhoods_retained / n(),
        frac_cells_retained = frac_included_cells(self$sce,
          allowed_Nhoods =
            which(NN_prob > 0 & Nhood_size > min_Nhood_size))
      ) %>%
      dplyr::filter(Nhoods_retained > 0 & med_NN_prob >= .01) %>%
      { . }

    self$NN_probs <-
      NN_probs %>% dplyr::right_join(gamma_stats, by = 'ltb')
  },

  overlay_cluster_cuts = function(heatmap_name) {
    decorate_column_dend(
      heatmap = heatmap_name,
      code = {
        N_clusters <- max(as.integer(self$tree_cut))
        stopifnot(is.finite(N_clusters))
        xl = sapply(unique(self$ind), function(i)
          data.table::first(which(self$ind == i)) - 1)
        xr = sapply(unique(self$ind), function(i)
          data.table::last(which(self$ind == i)))
        ##
        cols <- rep(okabe,
          ceiling(N_clusters/length(okabe)))[1:N_clusters]
        grid.rect(
          x = xl/length(self$ind),
          width = (xr - xl)/length(self$ind),
          just = 'left',
          default.units = 'npc',
          gp = gpar(fill = add_alpha(cols, 0.5), col = NA)
        )
        grid.text(
          label = unique(self$ind),
          x = 1/2*(xl+xr)/length(self$ind),
          y = .5,
          just = 'center',
          default.units = 'npc',
          gp = gpar(fontsize = 6)
        )
      }
    )
  },

  overlay_cluster_labels = function(x_loc = .9, y_loc = .1) {
    if (T || !is.null(self$E_tab_p)) {
      self$compute_cluster_enrichment()
    }
    y_loc = max(y_loc, nrow(self$E_tab_p)*.025)
    pushViewport(viewport(
        x = unit(x_loc, 'npc'),
        y = unit(y_loc, 'npc'),
        width = unit(1-x_loc, 'npc'),
        height = unit(y_loc, 'npc'),
        just = c(1, 1)))
    # grid.rect(gp=gpar(fill="blue"))
    draw_table(self$E_tab_p[unique(self$ind), ])
  },

  compute_cluster_enrichment = function() {
    self$ref_sa <- self$ref_sa %>%
      dplyr::mutate(stim_no_conc = factor(case_when(
            tnf_conc > 0 & ifn_conc > 0 ~ 'combo',
            tnf_conc > 0 & ifn_conc == 0 ~ 'TNFa',
            tnf_conc == 0 & ifn_conc > 0 ~ 'IFNy',
            tnf_conc == 0 & ifn_conc == 0 ~ 'Unstim'
        ))
      )

    f_ref_sa <- self$ref_sa %>%
      { .[match(colnames(self$affM), rownames(.)), ] } %>%
      { . }

    ## large enough clusters
    LEC <- names(which(map_int(auto_name(unique(self$tree_cut)),
          ~sum(self$tree_cut == .x)) > 1))
    ## Remove NAs
    rn <- function(x) {
      if (is.na(x)) return(0)
      else return(x)
    }
    ## log2 fraction success
    l2fs <- function(ta)
      log2(rn(ta['TRUE'])+1) - log2(sum(ta, na.rm = T)+1)
    ## Cluster enrichment of factors of interest
    self$E_tab <-
      tidyr::expand_grid(
        cluster = LEC,
        var = c('duration', 'stim_no_conc')
      ) %>%
      purrr::pmap_dfr(function(cluster, var) {
        tibble::tibble(
          cluster = cluster,
          var = var,
          lev = levels(self$ref_sa[, var]),
          purrr::map_dfr(lev, function(l) {
            cl_tab <- table(f_ref_sa[which(self$tree_cut == cluster), var] == l)
            all_tab <- table(self$ref_sa[, var] == l)
            tibble(
              OE = l2fs(cl_tab) - l2fs(all_tab),
              p = phyper(
                q = rn(cl_tab['TRUE'])-1,
                m = rn(all_tab['TRUE']),
                n = rn(all_tab['FALSE']),
                k = sum(cl_tab, na.rm = T),
                lower.tail=FALSE
              )
            )
          })
        )
      })

    if ('p' %in% colnames(self$E_tab)) {
      self$E_tab <- dplyr::mutate(self$E_tab,
        p_adj = p.adjust(p))
    }

    if ('cluster' %in% colnames(self$E_tab)) {
      self$E_tab <- self$E_tab %>%
      dplyr::mutate(cluster = as.integer(cluster)) %>%
      dplyr::arrange(cluster)
    }

    self$E_tab_p <- tryCatch({
      self$E_tab %>%
      dplyr::filter(p_adj <= .25 & OE > 0) %>%
      # dplyr::filter(OE >= 0.2) %>%
      dplyr::mutate(lev = ifelse(
          var == 'duration', paste0(lev, 'h'), lev)) %>%
      dplyr::group_by(cluster) %>%
      dplyr::summarize(label = paste(lev, collapse = ' - ')) %>%
      dplyr::ungroup() %>%
      { .[naturalsort::naturalorder(.$cluster), ] } %>%
      { . }

    }, error = function(e) { print(e); NULL })
  },

  compute_rep_ref_distM = function(ltb = 20) {
    stopifnot(!is.null(self$tree_cut))

    self$ref_cluster_centroids <-
      self$dtf[match(names(self$tree_cut), self$dtf$rn), ] %>%
      dplyr::select(matches('^CF')) %>%
      dplyr::mutate(sample_name = glue::glue('rc_{self$tree_cut}')) %>%
      dplyr::group_by(sample_name) %>%
      dplyr::summarize(across(matches('CF'), median)) %>%
      dplyr::ungroup()

    query_Nhood_coords <-
      self$dtf %>%
      dplyr::filter(experiment == self$query) %>%
      dplyr::select(matches('^CF'), rn) %>%
      dplyr::rename(sample_name = rn) %>%
      { . }

    ## Find the ref sample is most closely located to the centroids
    ref_and_centroids <-
      self$dtf[match(names(self$tree_cut), self$dtf$rn), ] %>%
      dplyr::select(matches('^CF'), rn) %>%
      dplyr::rename(sample_name = rn) %>%
      rbind(self$ref_cluster_centroids)

    rep_ref_sample_idxs <-
      ref_and_centroids %>%
      dplyr::select(matches('^CF')) %>%
      dist() %>%
      as.matrix() %>%
      {
        .[!ref_and_centroids$sample_name %in% names(self$tree_cut),
          ref_and_centroids$sample_name %in% names(self$tree_cut)]
      } %>%
      apply(1, which.min) %>%
      { . }

    self$rep_ref_samples <-
      self$dtf[match(names(self$tree_cut), self$dtf$rn), ] %>%
      { .[rep_ref_sample_idxs, ] }

    ## Combine the query Nhoods and rep ref samples in order to
    ## compute distances between them
    dtf_c <-
      self$rep_ref_samples %>%
      dplyr::select(matches('^CF'), rn) %>%
      dplyr::rename(sample_name = rn) %>%
      rbind(query_Nhood_coords)

    ## ith column corresponds to ith ref sample cluster
    self$rep_ref_affM <-
      dtf_c %>%
      { set_rownames(as.matrix(dplyr::select(., matches('^CF'))),
        .$sample_name) } %>%
      dist() %>%
      as.matrix() %>%
      {
        .[!dtf_c$sample_name %in% names(self$tree_cut),
          dtf_c$sample_name %in% names(self$tree_cut)]
      } %>%
      { exp(-ltb * .) } %>%
      { set_colnames(., glue::glue('rc_{1:ncol(.)}')) } %>%
      { . }

    self$dtf_u <-
      self$dtf %>%
      dplyr::filter(experiment == self$query) %>%
      bind_cols(self$rep_ref_affM)
  },

  plot_ref_affM_UMAP = function() {
    if (is.null(self$dtf_u)) {
      self$compute_rep_ref_distM()
    }

    plots <-
      unique(self$ind) %>%
      purrr::map(function(i) {
        ref_sample <- self$rep_ref_samples[i, ]
        plot_dim_reduc(dtf = self$dtf_u,
          colour_var = colnames(self$rep_ref_distM)[i]) +
          annotate_npc(
            npcx = .9, npcy = .9,
            hjust = 1, vjust = 1,
            label = unlist(self$E_tab_p[i, 'label'])
          ) +
          geom_point(data = ref_sample,
            colour = 'indianred3', size = 3, alpha = .5) +
          geom_point(data = ref_sample,
            colour = 'indianred3', size = 2, alpha = .7) +
          remove_x +
          remove_y +
          theme(aspect.ratio = 1) +
          guides(colour = guide_colourbar())
      })
    pe <- rlang::expr({print(patchwork::wrap_plots(plots))})
    N_plots <- length(plots)
    self$distM_UMAP_plot_height <- 12 * (N_plots %% 2 + 1)
    print_plot_eval(!!pe,
      width = 17.4,
      height = self$distM_UMAP_plot_height,
      filename = file.path(Sys.getenv('img_dir'),
        glue::glue('rep_ref_sims_UMAP.pdf')))
  },

  compute_UMAP_coords = function(cap = 1e3) {
    self$full_distM <-
      self$dtf %>%
      dplyr::select(any_of(paste0('CF', 1:self$N_PV))) %>%
      as.matrix() %>%
      dist() %>%
      as.matrix() %>%
      # { exp(-self$ltb * .) } %>%
      # { exp(-self$ltb * .) } %>%
      # { 1 / . } %>%
      set_rownames(self$dtf$rn) %>%
      set_colnames(self$dtf$rn) %>%
      { . }
    # self$full_distM <- 1 / self$full_distM
    # self$full_distM[self$full_distM >= cap] <- cap
    self$full_distM <- self$full_distM /
      quantile(self$full_distM[is.finite(self$full_distM)], .9)
    self$full_distM[!is.finite(self$full_distM)] <- 1e9

    library(reticulate)
    umap <- import('umap')
    umap_pc <- umap$UMAP(metric = 'precomputed')
    umap_pc <- umap_pc$fit_transform(self$full_distM)
    colnames(umap_pc) <- c('UMAP1', 'UMAP2')
    umap_pc
  },

  simplify_affM = function(max_allowed_size = 1e-2) {
    affM_S <- self$affM
    affM_S[affM_S < max_allowed_size] <- 0
    dn <- dimnames(affM_S)
    affM_S <- affM_S %>%
      apply(1, function(x) x / sum(x, na.rm = T)) %>%
      t()
    affM_S[!is.finite(self$affM_S)] <- 0
    dimnames(affM_S) <- dn
    # apply(nrs$affM, 1, max)
    affM_S <- affM_S[
      apply(affM_S, 1, function(x) all(is.finite(x))), ]
    self$affM_S <- affM_S
  },

  plot_Nhood_ref_affinity_graph = function() {
    library('ggraph')
    library('tidygraph')

    query_sa <-
      SummarizedExperiment::colData(self$sce) %>%
      as.data.frame() %>%
      set_rownames(NULL) %>%
      extract_sa(meta_fields = c('condition_i', 'condition_name')) %>%
      dplyr::distinct() %>%
      # order_condition_name() %>%
      dplyr::arrange(condition_i) %>%
      { . }

    ## Condition names of the query experiment
    query_cn <-
      query_sa[, c('condition_i', 'condition_name')] %>%
      dplyr::distinct() %>%
      dplyr::arrange(condition_i) %>%
      pull(condition_name)

    ref_idxs <- which(apply(self$affM_S, 2, sum) > 0)

    SC_to_Nhoods <- as.matrix(nhoodCounts(self$sce))
    SC_to_Nhoods <-
      SC_to_Nhoods[, naturalsort::naturalsort(colnames(SC_to_Nhoods))]
    rownames(SC_to_Nhoods) <- unlist(nhoodIndex(self$sce))

    SC_to_Nhoods_graph <-
      purrr::map(1:ncol(SC_to_Nhoods), function(i) {
        idxs <- which(SC_to_Nhoods[, i] > 0)
        tbl_graph(
          edges =
            tibble(
              from = paste0('CN', i),
              to = rownames(SC_to_Nhoods)[idxs],
              affinity = SC_to_Nhoods[idxs, i] /
                sum(SC_to_Nhoods[idxs, i])
            ),
          nodes =
            tibble(
              name = c(paste0('CN', i), rownames(SC_to_Nhoods)[idxs]),
              lvl = c(3, rep(2, length(idxs)))
            )
        )
      }) %>%
      purrr::discard(is.null) %>%
      purrr::reduce(graph_join, by = c('name', 'lvl')) %>%
      { . }

    graph <-
      ref_idxs %>%
      purrr::map(function(i) {
        cur <- self$affM_S[, i] %>% { .[. > 0.05] }
        if (!length(cur)) return(NULL)
        tbl_graph(
          edges =
            tibble(
              from = names(cur),
              to = colnames(self$affM_S)[i],
              affinity = unname(cur)
            ),
          nodes =
            tibble(
              name = c(names(cur), colnames(self$affM_S)[i]),
              type = c(rep('SC Nhood', length(cur)), 'reference'),
              lvl = c(rep(2, length(cur)), 1),
              label = c(rep('', length(cur)), colnames(self$affM_S)[i]),
              stim_group = c(rep(NA_character_, length(cur)),
                as.character(ref_sa[match(colnames(self$affM_S)[i],
                    ref_sa$condition_name), 'stim_group'])),
              duration = c(rep('Unknown', length(cur)),
                as.character(ref_sa[match(colnames(self$affM_S)[i],
                    ref_sa$condition_name), 'duration']))
            )
        )
      }) %>%
      purrr::discard(is.null) %>%
      purrr::reduce(graph_join,
        by = c('name', 'type', 'lvl', 'label',
          'stim_group', 'duration')) %>%
      graph_join(SC_to_Nhoods_graph, by = c('name', 'lvl')) %>%
      { . }

    # as_tibble(activate(graph, nodes))
    # ref_sa$duration
    # ref_sa$stim_group
    # table(as_tibble(activate(graph, nodes))$duration)
    # table(as_tibble(activate(graph, nodes))$stim_group)

    if (F) {
      graph <- graph_join(graph,
        tbl_graph(
          edges =
            tidyr::expand_grid(
              tibble(
                i = ref_idxs,
                from = colnames(affM_S)[i],
              ),
              tibble(
                j = ref_idxs,
                to = colnames(affM_S)[j],
              ),
              affinity = 0.01
            ) %>%
            dplyr::filter(from < to) %>%
            dplyr::select(-i, -j),
          nodes =
            tibble(
              name = colnames(self$affM_S)[ref_idxs],
              type = c('reference'),
              lvl = 1
            )
        )
      )
    }

    cl <- igraph::components(graph)
    graph <- activate(graph, nodes) %>%
      mutate(cluster = cl$membership)
    graph <- activate(graph, nodes) %>%
      filter(cluster == 1) %>%
      # as_tibble()
      { . }

    pacman::p_load('ggnetwork')

    # p <- ggraph(graph, 'centrality',
    #   centrality = igraph::graph.strength(graph)) +
    xy <- layout_as_multilevel(graph,
      FUN1 = igraph::layout_on_grid,
      FUN2 = igraph::layout_on_grid,
      type = 'all'
      # type = 'separate',
      # , alpha = 25, beta = 45
    )
    p <- ggraph(graph, 'manual', x = xy[, 1], y = xy[, 2]) +
      geom_edge_link(aes(edge_width = affinity), colour = 'grey50',
        alpha = .2) +
      geom_node_point(
        mapping = aes(
          size = igraph::graph.strength(graph),
          shape = duration,
          colour = stim_group
        ), alpha = .5
      ) +
      # ggnetwork::geom_nodetext_repel(
      #   aes(label = duration, x = x, y = y),
      #   size = 2) +
      # geom_node_point(data = activate(ref_nodes, nodes)) +
      scale_size_continuous(name = '# proximal SC Nhoods') +
      scale_colour_stim_group(as_tibble(activate(graph, nodes))) +
      scale_shape_duration(as_tibble(activate(graph, nodes))) +
      scale_edge_width_continuous(range = c(0.2, 1)) +
      coord_fixed() +
      theme_graph(base_family = 'Helvetica') +
      set_graph_style(family = 'Helvetica') +
      theme_cyto_inf(
        legend.direction = 'vertical',
        legend.position = 'right'
      ) +
      gg_tabula_rasa
    print_plot_eval(print(p),
      width = 17.4, height = 15,
      filename = file.path(out_dir,
        glue::glue('graph.pdf')))

    return(p)
  },

  compute_ref_tree = function() {
    self$ref_M_tree <- gen_clust_object(
      M = t(self$ref_M), dist_f = 'euclidean',
      clust_method = 'complete'
    )
    # dend <- as.dendrogram(self$ref_M_tree)
    # self$ind = self$tree_cut[order.dendrogram(dend)]
  }

  )
)



#' Simplify a collection of reference samples using Euclidean distance
#'
#' Reduce the number of reference samples by specifying a set of
#' 'anchor' samples and assessing their Euclidean distance to the
#' remaining 'query' samples. Any of the query samples sufficiently
#' far away from the anchor samples are kept.
#'
#' @param M Matrix of sample coordinates in latent or raw space. Row
#' names must be sample identifiers
#' @param sa Sample annotation. At least one the columns must contain
#' the rownames of M. 
#' @param sample_name_regex Regex to use to identify 'anchor' samples
#' @param sample_names Vector of anchor sample names. Either this
#' argument of \code{sample_name_regex} mus be NULL
#' @param max_tol_aff Maximum tolerated affinity for a query sample
#' and any of the anchor samples in order for the query sample not to
#' be thrown out.
#' @param verbose Whether to print intermediate information
#' @return A list with the cleaned up matrix M and sample annotation
#' sa
exclude_ref_sample_cor <- function(
  M, sa,
  # sample_name_regex = '(U|u)nstim',
  sample_name_regex = NULL,
  sample_names = NULL,
  max_tol_aff = .05,
  verbose = FALSE) {

  stopifnot(!is.null(M))
  stopifnot(!is.null(sa))
  stopifnot(xor(is.null(sample_name_regex), is.null(sample_names)))

  ## Make sure 'M' and 'sa' are in identical order
  match_col <- find_match_col(rownames(M), dtf = sa)
  idxs <- match(tolower(sa[[match_col]]), tolower(rownames(M)))
  M <- M[idxs, ]

  aff_M <- 1/as.matrix(dist(M))

  if (!is.null(sample_name_regex)) {
    ## Indices of samples that are to be compared against
    ref_idx <-
      which(stringr::str_detect(rownames(aff_M), sample_name_regex))
    aff_M <- aff_M[ref_idx, -ref_idx, drop=F]
  } else if (!is.null(sample_names)){
    sample_names <- intersect(sample_names, colnames(aff_M))
    aff_M <-
      aff_M[sample_names, setdiff(colnames(aff_M), sample_names),
      drop=F]
  } else {
    rlang::abort('Define either sample_name_regex or sample_names')
  }

  ## Compute the max affinity for each tested sample to each of the
  ## 'reference' samples
  if (nrow(aff_M) > 1) {
    max_aff <- apply(aff_M, 2, max)
  } else {
    max_aff <- aff_M[1, ]
  }

  if (is.character(max_tol_aff)) {
    if (!max_tol_aff %in% names(max_aff)) {
      browser()
    }
    max_tol_aff <- unname(max_aff[max_tol_aff]) - 1e-8
    message(max_tol_aff)
  }

  # exclude_sn <- names(which(max_aff >= .35))
  exclude_sn <- names(which(max_aff >= max_tol_aff))
  match_col <- find_match_col(exclude_sn, dtf = sa)
  exclude_idxs <- match(tolower(exclude_sn), tolower(sa[[match_col]]))

  if (length(exclude_idxs) > 0) {
    if (verbose) {
      cat('Including:\n')
      print(max_aff[setdiff(names(max_aff), exclude_sn)])
      if (!is.null(sample_names))
        print(sample_names)
      cat('Excluding:\n')
      # id_vars <- c('stim_group', 'duration',
      #         'condition_name') %>%
      id_vars <- c('stim_group', 'duration', 'condition_name') %>%
        intersect(colnames(sa))
      # if (all(c('stim_group', 'condition_name') %in% colnames(sa))) {
      #   id_vars <- setdiff(id_vars, 'condition_name')
      # }
      if ('condition_name' %in% colnames(sa)) {
        id_vars <- 'condition_name'
      }
      set_rownames(sa[exclude_idxs, ], NULL) %>%
        dplyr::select(any_of(id_vars)) %>%
        arrange_at(id_vars) %>%
        print()
      # exclude_idxs <- match(exclude_sn, sa$condition_name)
    }

    if (!is.null(sa)) {
      sa <- sa[-exclude_idxs, ]
      if ('stim_group' %in% colnames(sa))
        sa$stim_group <- droplevels(sa$stim_group)
    } else {
      sa <- NULL
    }
    M <- M[-exclude_idxs, ]
  }

  return(list(M = M, sa = sa))
}
