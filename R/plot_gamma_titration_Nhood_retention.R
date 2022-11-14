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
