% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare_settings.R
\name{extract_best}
\alias{extract_best}
\title{Sort a set of tier1 objects by a ranking function and extract the
i-th item, along with the error metrics for this item}
\usage{
extract_best(
  bm_res = NULL,
  bm_res_sum = NULL,
  filtering_fun = filtering_funs[[1]],
  rank_fun = rank_funs[[1]],
  tier1_rank = 1L
)
}
\description{
Sort a set of tier1 objects by a ranking function and extract the
i-th item, along with the error metrics for this item
}
