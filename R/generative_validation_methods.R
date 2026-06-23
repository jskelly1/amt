#' @export
print.amt_generative_validation <- function(x, ...) {
  cat("amt generative validation\n")
  cat("Observed locations:", x$observed_summary$n_locations, "\n")
  cat("Simulated tracks:", x$simulated_summary$n_simulations, "\n")
  cat("Validation pillars:", paste(generative_metric_pillar(x$metrics), collapse = ", "), "\n\n")
  print(generative_print_table(x), n = Inf, width = Inf)
  invisible(x)
}

#' @export
summary.amt_generative_validation <- function(object, ...) {
  rows <- lapply(names(object$metric_results), function(metric_name) {
    result <- object$metric_results[[metric_name]]
    observed_statistic <- switch(
      metric_name,
      ud = result$statistic,
      msd = result$statistic,
      straightness = result$observed_value,
      barrier = result$observed_count,
      NA_real_
    )

    discrepancy_statistic <- switch(
      metric_name,
      ud = result$statistic,
      msd = result$statistic,
      straightness = result$statistic,
      barrier = result$observed_count,
      NA_real_
    )

    statistic_name <- switch(
      metric_name,
      ud = "mean observed-simulated grid W1",
      msd = "MSD integrated squared error",
      straightness = "absolute straightness deviation",
      barrier = "segment-intersection count",
      NA_character_
    )

    tibble::tibble(
      metric = metric_name,
      statistic_name = statistic_name,
      observed_statistic = observed_statistic,
      discrepancy_statistic = discrepancy_statistic,
      p_value = generative_metric_p_value(result),
      alternative = result$rank_test$alternative %||% NA_character_
    )
  })

  tibble::as_tibble(do.call(rbind, lapply(rows, as.data.frame)))
}

generative_print_table <- function(x) {
  summary_tbl <- summary(x)
  tibble::tibble(
    pillar = generative_metric_pillar(summary_tbl$metric),
    metric = summary_tbl$metric,
    diagnostic = summary_tbl$statistic_name,
    observed = summary_tbl$observed_statistic,
    discrepancy = summary_tbl$discrepancy_statistic,
    p_value = summary_tbl$p_value,
    alternative = summary_tbl$alternative
  )
}

#' Plot generative validation diagnostics
#'
#' Produces `ggplot2` summaries for available generative validation metrics.
#' `ggplot2` is optional and must be installed to use this method.
#'
#' @param x An object returned by [validate_generative()].
#' @param metric Optional metric to plot. If `NULL`, all available plots are
#'   printed and returned as a named list.
#' @param ... Currently unused.
#'
#' @return A ggplot object when `metric` is supplied, otherwise an invisible
#'   named list of ggplot objects.
#' @export
plot.amt_generative_validation <- function(x, metric = NULL, ...) {
  require_generative_ggplot2()
  plots <- lapply(x$metric_results, plot_generative_metric)

  if (!is.null(metric)) {
    metric <- match.arg(metric, names(plots))
    print(plots[[metric]])
    return(invisible(plots[[metric]]))
  }

  for (plot_i in plots) {
    print(plot_i)
  }
  invisible(plots)
}

plot_generative_metric <- function(x) {
  if (inherits(x, "amt_generative_metric_ud")) {
    return(plot_generative_metric_ud(x))
  }
  if (inherits(x, "amt_generative_metric_msd")) {
    return(plot_generative_metric_msd(x))
  }
  if (inherits(x, "amt_generative_metric_straightness")) {
    return(plot_generative_metric_straightness(x))
  }
  if (inherits(x, "amt_generative_metric_barrier")) {
    return(plot_generative_metric_barrier(x))
  }
  stop("No plot method is available for this metric result.", call. = FALSE)
}

plot_generative_metric_ud <- function(x) {
  cols <- generative_plot_colours()
  ggplot2::ggplot(x$simulated_statistics, ggplot2::aes(x = .data$statistic)) +
    ggplot2::geom_histogram(bins = 15, fill = cols$simulated,
                            color = "white", alpha = 0.75) +
    ggplot2::geom_vline(xintercept = x$statistic, color = cols$observed,
                        linewidth = 1) +
    theme_generative_diagnostic() +
    ggplot2::labs(
      x = "Mean 1-Wasserstein distance",
      y = "Number of simulations",
      title = "Emergent utilization distribution",
      subtitle = paste0("Monte Carlo p = ", format_generative_p(generative_metric_p_value(x)))
    )
}

plot_generative_metric_msd <- function(x) {
  cols <- generative_plot_colours()
  ggplot2::ggplot(x$envelope, ggplot2::aes(x = .data$lag)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$lo, ymax = .data$hi),
      fill = cols$simulated,
      alpha = 0.22
    ) +
    ggplot2::geom_line(ggplot2::aes(y = .data$mean), color = cols$simulated,
                       linetype = "dashed") +
    ggplot2::geom_line(
      data = x$observed_curve,
      ggplot2::aes(y = .data$msd),
      color = cols$observed,
      linewidth = 0.9
    ) +
    theme_generative_diagnostic() +
    ggplot2::labs(
      x = "Lag",
      y = "Mean squared displacement",
      title = "Mean squared displacement",
      subtitle = paste0("Monte Carlo p = ", format_generative_p(generative_metric_p_value(x)))
    )
}

plot_generative_metric_straightness <- function(x) {
  cols <- generative_plot_colours()
  ggplot2::ggplot(x$simulated_values, ggplot2::aes(x = .data$value)) +
    ggplot2::geom_histogram(bins = 15, fill = cols$simulated,
                            color = "white", alpha = 0.75) +
    ggplot2::geom_vline(xintercept = x$observed_value, color = cols$observed,
                        linewidth = 1) +
    theme_generative_diagnostic() +
    ggplot2::labs(
      x = "Straightness index",
      y = "Number of simulations",
      title = "Path straightness",
      subtitle = paste0("Monte Carlo p = ", format_generative_p(generative_metric_p_value(x)))
    )
}

plot_generative_metric_barrier <- function(x) {
  cols <- generative_plot_colours()
  ggplot2::ggplot(x$simulated_counts, ggplot2::aes(x = .data$count)) +
    ggplot2::geom_histogram(
      binwidth = 1,
      fill = cols$simulated,
      color = "white",
      boundary = -0.5,
      alpha = 0.75
    ) +
    ggplot2::geom_vline(xintercept = x$observed_count, color = cols$observed,
                        linewidth = 1) +
    theme_generative_diagnostic() +
    ggplot2::labs(
      x = "Movement segments intersecting barrier",
      y = "Number of simulations",
      title = "Barrier interactions",
      subtitle = paste0("Monte Carlo p = ", format_generative_p(generative_metric_p_value(x)))
    )
}

require_generative_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The `ggplot2` package is required for this plot. Install it to continue.",
         call. = FALSE)
  }
}

generative_plot_colours <- function() {
  list(
    simulated = "#4C78A8",
    observed = "#D55E00",
    grid = "grey90"
  )
}

theme_generative_diagnostic <- function(base_size = 11) {
  cols <- generative_plot_colours()
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = cols$grid, linewidth = 0.3),
      plot.title = ggplot2::element_text(face = "plain"),
      legend.position = "bottom"
    )
}
