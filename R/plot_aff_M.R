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


#' 
#'
#' @export
plot_aff_M.NhoodRefSim <- function(
    obj,
    compo_HM_data = NULL,
    out_dir = Sys.getenv('img_dir'),
    plot_id = '',
    ## What fraction of Nhoods should be left out, based on their
    ## affinity for any reference sample?
    min_aff_q = NULL,
    order_ref_columns = 'Nhood_diagonal',
    ...
  ) {

  of <- file.path(out_dir, glue::glue('affM_HM{prepend_hyphen(plot_id)}.pdf'))
  obj$compute_affM()
  # obj$simplify_affM()

  if (!is.null(min_aff_q) && min_aff_q > 0) {
    min_aff <- quantile(obj$max_aff, min_aff_q)
    idxs <- which(obj$max_aff >= min_aff)
    l_affM <- obj$affM[idxs, ]
    allowed_Nhoods <-
      intersect(rownames(obj$affM), rownames(l_affM))
  } else {
    l_affM <- obj$affM
    allowed_Nhoods <- rownames(obj$affM)
  }
  if (F) {
    rowSums(l_affM)
    colSums(l_affM)
  }

  if (is.null(compo_HM_data)) {
    NH_composition_HM <- gen_NH_CN_HM(
      sce = obj$sce,
      dtf = obj$dtf,
      side = 'left',
      allowed_Nhoods = allowed_Nhoods,
      show_row_dend = FALSE,
      show_column_names = FALSE,
      # ref_experiment = obj$ref_experiment
      forego_sa = FALSE,
    )
    if (interactive() && F) {
      print_plot_eval(draw(NH_composition_HM),
        width = 17.4, height = 10,
        filename = file.path(out_dir, glue::glue('NH_comp.pdf')))
    }
    main_heatmap <- 'Neighbourhood abundance\nof condition'
  } else {
    compo_HM_data <- compo_HM_data[
      match(rownames(l_affM), rownames(compo_HM_data)),
      , drop = F]
    # row_clustering <- gen_clust_object(t(compo_HM_data[,
    #   stringr::str_detect(colnames(compo_HM_data), 'order')]))
    # row_clustering <- gen_clust_object(t(compo_HM_data[,
    #   !stringr::str_detect(colnames(compo_HM_data), 'order')]))
    NH_composition_HM <- gen_HM(
      compo_HM_data[,
      !stringr::str_detect(colnames(compo_HM_data), 'order')],
      cluster_rows = T,
      show_row_dend = F,
      show_column_dend = F,
      show_column_names = T,
      name = 'DA logFC'
    )
    if (T) {
      print_plot_eval(draw(NH_composition_HM),
        width = 17.4, height = 10,
        filename = file.path(out_dir,
          glue::glue('milo_coef_test.pdf')))
    }

    main_heatmap <- 'DA logFC'
  }

  if (order_ref_columns == 'cluster') {
    cluster_columns <- gen_clust_object(obj$affM,
      dist_f = 'euclidean')
    ref_sample_order <- colnames(l_affM)
  } else if (order_ref_columns == 'Nhood_diagonal') {
    cluster_columns <- FALSE
    ## Once the rows are ordered, how do we order the columns such
    ## that the sparse matrix most strongly resembles a diagonal
    ## matrix?
    print_plot_eval({ r_order <<- row_order(NH_composition_HM) },
      width = 17.4, height = 10,
      filename = file.path(out_dir,
        glue::glue('temp.pdf')))
    orig <- l_affM[r_order, ]
    tr <- diag(nrow(l_affM):1-floor(nrow(l_affM)/2)) %*% orig
    ref_sample_order <- names(sort(apply(tr, 2, sum)))
  } else if (order_ref_columns == 'natural') {
    cluster_columns <- FALSE
    ref_sample_order <- colnames(l_affM)
  } else if (order_ref_columns == 'ifny_duration_tnfa') {
    match_col <- find_match_col(colnames(l_affM), dtf = obj$ref_sa)
    ref_sample_order <-
      dplyr::arrange(obj$ref_sa, ifn_conc, duration, tnf_conc) %>%
      pull(match_col) %>%
      intersect(colnames(l_affM))
    cluster_columns <- FALSE
  }

  if (F) {
    gst2gs <-
      list(
        'default' = tar_read(default_gene_sets)
      ) %>%
      # unlist(recursive = F) %>%
      purrr::flatten() %>%
      { .[!stringr::str_detect(names(.), 'HK')] } %>%
      { .[!stringr::str_detect(names(.), 'CXCL10')] } %>%
      { .[!stringr::str_detect(names(.), '6h.max')] } %>%
      { .[!stringr::str_detect(names(.), 'TNFa.plateau')] } %>%
      { . }
    MRs <- unname(unlist(gst2gs)) %>%
      intersect(rownames(nhoodExpression(obj$sce)))

    geneset_labels = tibble::enframe(gst2gs, 'gs', 'gene') %>%
      dplyr::mutate(gs = factor(gs,
          levels = c('TNFa', 'IFNy', 'synergy'))) %>%
      tidyr::unnest(gene) %>%
      dplyr::right_join(tibble(gene = MRs))


    M <- t(subset_feats(nhoodExpression(obj$sce), MRs))
    M <- scale(M)

    set.seed(3)
    gene_tree <- gen_clust_object(M[r_order, ])
    HA <- HeatmapAnnotation(gs = geneset_labels$gs, which = 'col')
    HM <- gen_HM(M[r_order, ], cluster_rows = F,
      cluster_columns = gene_tree)
    print_plot_eval({
      draw(HA %v% HM)
    }, width = 17.4, height = 10,
      filename = file.path(out_dir,
        glue::glue('MR_HM.pdf')))
  }

  cell_size = 18/max(ncol(l_affM), nrow(l_affM))
  aff_name <- 'Query sample to reference sample\naffinity'
  HM <- gen_HM(
    l_affM[, ref_sample_order],
    # ca = obj$subset_ref_sa(colnames(l_affM)),
    ca = obj$subset_ref_sa(ref_sample_order),
    column_annotation_name_side = 'left',
    row_dend_side = 'left',
    # width = unit(ncol(l_affM)*cell_size*.5, 'cm'),
    # height = unit(nrow(l_affM)*cell_size*.5, 'cm'),
    # width = unit(10, 'cm'),
    # height = unit(16, 'cm'),
    show_row_names = F,
    show_column_names = F,
    cluster_columns = F,
    # cluster_columns = T,
    # cluster_columns = obj$ref_M_tree,
    column_dend_reorder = FALSE,
    # cluster_rows = T,
    show_column_dend = T,
    show_row_dend = T,
    name = aff_name,
    column_gap = unit(1, 'mm'),
    row_gap = unit(1, 'mm'),
    ...
  )

  RA <- rowAnnotation(
    `Max affinity (log10)` = log10(obj$max_aff[rownames(l_affM)]),
    `Nhood size (log10)` =
      log10(obj$Nhood_size[rownames(l_affM)]),
    # annotation_name_gp = legend_params,
    annotation_legend_param = legend_params
  )
  print_plot_eval(
    {
      draw(
        NH_composition_HM + HM,
        merge_legend = F,
        # main_heatmap = aff_name,
        main_heatmap = main_heatmap,
        heatmap_legend_side = 'bottom'
      )
      # obj$overlay_cluster_cuts(heatmap_name = aff_name)
      # obj$overlay_cluster_labels()
    },
    width = 17.4,
    # height = min(25, nrow(l_affM)/1.2),
    height = min(15, nrow(l_affM)/1.2),
    filename = of
  )
}


