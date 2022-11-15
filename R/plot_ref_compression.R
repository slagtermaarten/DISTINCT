## TODO fix me
plot_reference_compression <- function(obj) {
  sa <- extract_sa(extract_ref_dtf(obj$dtf),
    meta_fields = c('condition_name', 'stim_group', 'duration')
  )

  p_dat <- obj$ref_reconstruction %>%
    dplyr::left_join(sa, by = 'condition_name') %>%
    dplyr::mutate(log_gamma = round(log_gamma, 2)) %>%
    dplyr::filter(log_gamma %in% unique(log_gamma)[1:3])

  p <-
    p_dat %>%
    ggplot(aes(x = 1/NLPC, y = 1/CF, 
        colour = stim_group, shape = duration)) +
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
