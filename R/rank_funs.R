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
