library(amt)

observed <- track(c(0, 1, 1, 2), c(0, 0, 1, 1))
simulated <- list(
  sim_1 = track(c(0, 1, 2, 3), c(0, 0, 0, 0)),
  sim_2 = track(c(0, 0, 1, 1), c(0, 1, 1, 2)),
  sim_3 = track(c(0, 1, 1, 1), c(0, 0, 1, 2))
)

test <- mc_rank_test(5, c(1, 2, 3, 4), alternative = "greater")
expect_equal(test$p_value, 0.2)
expect_equal(test$rank, 5)
expect_equal(test$alternative, "greater")

res <- validate_generative(
  observed,
  simulated,
  metrics = c("ud", "msd", "straightness"),
  ud_args = list(grid_size = 4),
  msd_args = list(max_lag = 2)
)

expect_true(inherits(res, "amt_generative_validation"))
expect_equal(res$metrics, c("ud", "msd", "straightness"))
expect_true(is(summary(res), "tbl_df"))
expect_true(is.finite(res$metric_results$ud$statistic))
expect_equal(res$metric_results$ud$grid_size, c(4L, 4L))
expect_equal(res$metric_results$msd$max_lag, 2L)
expect_true(is.finite(res$metric_results$straightness$observed_value))

res_alias <- validate_generative(
  observed,
  simulated,
  metrics = "sinuosity"
)
expect_equal(res_alias$metrics, "straightness")

sim_df <- data.frame(
  sim_id = rep(c("a", "b"), each = 4),
  x_ = c(0, 1, 2, 3, 0, 0, 1, 1),
  y_ = c(0, 0, 0, 0, 0, 1, 1, 2)
)
res_df <- validate_generative(
  observed,
  sim_df,
  metrics = "straightness"
)
expect_equal(res_df$simulated_summary$n_simulations, 2)

steps_observed <- steps(observed)
res_steps <- validate_generative(
  steps_observed,
  simulated,
  metrics = "straightness"
)
expect_equal(res_steps$observed_summary$n_locations, nrow(observed))

barrier <- sf::st_sfc(
  sf::st_linestring(matrix(c(0.5, -1, 0.5, 2), ncol = 2, byrow = TRUE)),
  crs = 3857
)
barrier_observed <- track(c(0, 0.25, 0.4), c(0, 0, 0))
barrier_sims <- list(
  track(c(0, 1), c(0, 0)),
  track(c(0, 0.75), c(0.2, 0.2)),
  track(c(0, 0.25, 0.4), c(1, 1, 1))
)
barrier_res <- validate_generative(
  barrier_observed,
  barrier_sims,
  metrics = "barrier",
  barrier = barrier
)
expect_equal(barrier_res$metric_results$barrier$observed_count, 0L)
expect_equal(barrier_res$metric_results$barrier$simulated_counts$count, c(1L, 1L, 0L))

touching_track <- track(c(0, 0.5, 0.5), c(0, 0, 1))
touching_res <- validate_generative(
  touching_track,
  barrier_sims,
  metrics = "barrier",
  barrier = barrier
)
expect_equal(touching_res$metric_results$barrier$observed_count, 2L)

expect_error(validate_generative(observed, simulated[1], metrics = "msd"))
expect_error(validate_generative(observed, simulated, metrics = "unknown"))
expect_error(validate_generative(track(c(0, 0), c(0, 0)), simulated, metrics = "straightness"))

lonlat_track <- sf::st_sf(
  geometry = sf::st_sfc(
    sf::st_point(c(-71, 46)),
    sf::st_point(c(-71.1, 46.1)),
    crs = 4326
  )
)
expect_error(validate_generative(lonlat_track, simulated, metrics = "straightness"))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  p <- plot(res, metric = "msd")
  expect_true(inherits(p, "ggplot"))

  data(deer)
  sl_plot <- deer[1:20, ] |>
    steps_by_burst() |>
    random_steps() |>
    plot_sl(engine = "ggplot2", plot = FALSE)
  expect_true(inherits(sl_plot, "ggplot"))
}

