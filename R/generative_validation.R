#' Monte Carlo rank test
#'
#' Computes a rank-based Monte Carlo p-value with the standard +1 correction.
#' The helper compares one observed statistic with corresponding statistics
#' from simulated trajectories.
#'
#' @param observed `[numeric(1)]` Observed statistic.
#' @param simulated `[numeric]` Simulated statistics.
#' @param alternative `[character(1)]` One of `"two.sided"`, `"less"`, or
#'   `"greater"`.
#'
#' @return A list with the observed statistic, simulated statistics, rank,
#'   p-value, and alternative.
#' @export
#'
#' @references
#' North, B. V., Curtis, D., and Sham, P. C. (2002). A note on the calculation
#' of empirical P values from Monte Carlo procedures. *American Journal of
#' Human Genetics*, 71(2), 439-441. \doi{10.1086/341527}
#'
#' @examples
#' mc_rank_test(5, c(1, 2, 3, 4), alternative = "greater")
mc_rank_test <- function(observed, simulated,
                         alternative = c("two.sided", "less", "greater")) {
  alternative <- match.arg(alternative)

  if (!is.numeric(observed) || length(observed) != 1 || !is.finite(observed)) {
    stop("`observed` must be one finite numeric value.", call. = FALSE)
  }
  if (!is.numeric(simulated)) {
    stop("`simulated` must be numeric.", call. = FALSE)
  }

  simulated <- simulated[is.finite(simulated)]
  if (length(simulated) < 1) {
    stop("`simulated` must contain at least one finite value.", call. = FALSE)
  }

  k <- length(simulated)
  p_less <- (sum(simulated <= observed) + 1) / (k + 1)
  p_greater <- (sum(simulated >= observed) + 1) / (k + 1)

  p_value <- switch(
    alternative,
    less = p_less,
    greater = p_greater,
    two.sided = min(1, 2 * min(p_less, p_greater))
  )

  out <- list(
    observed = observed,
    simulated = simulated,
    rank = sum(simulated <= observed) + 1,
    p_value = p_value,
    alternative = alternative
  )
  class(out) <- c("amt_mc_rank_test", "list")
  out
}

#' Validate generative behavior of movement simulations
#'
#' Runs trajectory-level generative validation diagnostics by comparing an
#' observed track with simulated tracks. The fitted model is not required by
#' this function; users fit and simulate with `amt` or another workflow, then
#' pass the observed and simulated trajectories here.
#'
#' Coordinates are treated as planar Euclidean coordinates. Longitude/latitude
#' CRS are rejected when the CRS is known. If no CRS is available, coordinates
#' are assumed to be planar and already comparable across observed, simulated,
#' and barrier data.
#'
#' @param observed Observed track-like object. Supported inputs include
#'   `track_xy`, `track_xyt`, `steps_xyt`, `steps_xy`, and data frames with
#'   `x_` and `y_` columns.
#' @param simulated A list of simulated track-like objects, or a data frame with
#'   one row per simulated location and a simulation identifier column.
#' @param metrics `[character]` Metrics to compute. Supported values are
#'   `"ud"`, `"msd"`, `"straightness"`, `"sinuosity"`, and `"barrier"`.
#'   `"sinuosity"` is accepted as an alias for `"straightness"`.
#' @param barrier Optional `sf` LINESTRING or MULTILINESTRING object used for
#'   barrier crossing validation.
#' @param ud_args,msd_args,straightness_args,barrier_args Lists of additional
#'   arguments passed to the corresponding metric functions.
#' @param ... Passed to simulation coercion helpers.
#'
#' @return An object of class `amt_generative_validation`.
#' @export
#'
#' @references
#' Nicosia, A. (2026). Beyond the next step: A multi-criteria generative
#' validation framework for step selection functions. *Methods in Ecology and
#' Evolution*. \doi{10.1111/2041-210x.70313}
#'
#' @examples
#' observed <- track(c(0, 1, 1, 2), c(0, 0, 1, 1))
#' simulated <- list(
#'   track(c(0, 1, 2, 3), c(0, 0, 0, 0)),
#'   track(c(0, 0, 1, 1), c(0, 1, 1, 2))
#' )
#' res <- validate_generative(
#'   observed,
#'   simulated,
#'   metrics = c("msd", "straightness"),
#'   msd_args = list(max_lag = 2)
#' )
#' summary(res)
validate_generative <- function(
  observed,
  simulated,
  metrics = c("ud", "msd", "straightness", "barrier"),
  barrier = NULL,
  ud_args = list(),
  msd_args = list(),
  straightness_args = list(),
  barrier_args = list(),
  ...
) {
  call <- match.call()
  metrics <- validate_generative_metric_names(metrics)

  if ("barrier" %in% metrics && is.null(barrier)) {
    warning("Skipping `barrier` because `barrier = NULL`.", call. = FALSE)
    metrics <- setdiff(metrics, "barrier")
  }
  if (length(metrics) < 1L) {
    stop("No validation metrics remain after checking inputs.", call. = FALSE)
  }

  obs <- generative_track(observed)
  sims <- generative_simulations(simulated, ...)

  metric_results <- list()
  if ("ud" %in% metrics) {
    metric_results$ud <- do.call(
      validate_generative_ud,
      c(list(observed = obs, simulated = sims), ud_args)
    )
  }
  if ("msd" %in% metrics) {
    metric_results$msd <- do.call(
      validate_generative_msd,
      c(list(observed = obs, simulated = sims), msd_args)
    )
  }
  if ("straightness" %in% metrics) {
    metric_results$straightness <- do.call(
      validate_generative_straightness,
      c(list(observed = obs, simulated = sims), straightness_args)
    )
  }
  if ("barrier" %in% metrics) {
    metric_results$barrier <- do.call(
      validate_generative_barrier,
      c(list(observed = obs, simulated = sims, barrier = barrier), barrier_args)
    )
  }

  out <- list(
    call = call,
    observed_summary = generative_track_summary(obs),
    simulated_summary = generative_simulation_summary(sims),
    metrics = names(metric_results),
    metric_results = metric_results,
    settings = list(
      ud_args = ud_args,
      msd_args = msd_args,
      straightness_args = straightness_args,
      barrier_args = barrier_args
    )
  )
  class(out) <- c("amt_generative_validation", "list")
  out
}

validate_generative_metric_names <- function(metrics) {
  valid <- c("ud", "msd", "straightness", "sinuosity", "barrier")
  metrics <- unique(tolower(metrics))
  unknown <- setdiff(metrics, valid)

  if (length(unknown) > 0) {
    stop("Unknown metric(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }

  metrics[metrics == "sinuosity"] <- "straightness"
  unique(metrics)
}

validate_generative_ud <- function(observed, simulated, grid_size = 50,
                                   bounds = NULL, ...) {
  obs <- generative_track(observed)
  sims <- generative_simulations(simulated, ...)
  check_generative_simulation_count(sims, minimum = 2)

  grid_size <- normalize_generative_grid_size(grid_size)
  bounds <- validate_generative_bounds(bounds %||% combined_generative_bounds(obs, sims))
  check_generative_tracks_within_bounds(obs, sims, bounds)

  ud_obs <- empirical_grid_ud(obs, grid_size = grid_size, bounds = bounds)
  ud_sims <- lapply(sims, empirical_grid_ud, grid_size = grid_size, bounds = bounds)

  obs_distances <- vapply(
    ud_sims,
    function(ud_sim) wasserstein_ud(ud_obs, ud_sim),
    numeric(1)
  )
  check_finite_generative_vector(obs_distances, "Observed-simulated Wasserstein distances")

  n_sims <- length(ud_sims)
  distance_matrix <- matrix(0, nrow = n_sims, ncol = n_sims)
  for (i in seq_len(n_sims - 1L)) {
    for (j in seq.int(i + 1L, n_sims)) {
      distance_matrix[i, j] <- wasserstein_ud(ud_sims[[i]], ud_sims[[j]])
      distance_matrix[j, i] <- distance_matrix[i, j]
    }
  }

  sim_statistics <- vapply(seq_len(n_sims), function(i) {
    mean(distance_matrix[i, -i], na.rm = TRUE)
  }, numeric(1))
  check_finite_generative_vector(sim_statistics, "Simulation-to-simulation Wasserstein statistics")

  observed_statistic <- mean(obs_distances, na.rm = TRUE)
  check_finite_generative_vector(observed_statistic, "Observed UD discrepancy statistic")

  out <- list(
    metric = "ud",
    method = "empirical_grid_wasserstein",
    grid_size = grid_size,
    bounds = bounds,
    observed_ud = ud_obs,
    simulated_uds = ud_sims,
    observed_distances = tibble::tibble(
      sim_id = names(sims),
      distance = unname(obs_distances)
    ),
    simulated_statistics = tibble::tibble(
      sim_id = names(sims),
      statistic = unname(sim_statistics)
    ),
    statistic = observed_statistic,
    rank_test = mc_rank_test(observed_statistic, sim_statistics, alternative = "greater")
  )
  class(out) <- c("amt_generative_metric_ud", "list")
  out
}

validate_generative_msd <- function(observed, simulated, max_lag = NULL,
                                    envelope_probs = c(0.025, 0.975), ...) {
  obs <- generative_track(observed)
  sims <- generative_simulations(simulated, ...)
  check_generative_simulation_count(sims, minimum = 2)

  common_max_lag <- min(c(nrow(obs), vapply(sims, nrow, integer(1)))) - 1L
  max_lag <- validate_generative_max_lag(max_lag %||% common_max_lag, common_max_lag)
  envelope_probs <- validate_generative_envelope_probs(envelope_probs)

  obs_curve <- generative_msd_curve(obs, max_lag = max_lag)
  sim_curves <- lapply(sims, generative_msd_curve, max_lag = max_lag)
  sim_matrix <- do.call(cbind, lapply(sim_curves, `[[`, "msd"))
  colnames(sim_matrix) <- names(sims)

  sim_mean <- rowMeans(sim_matrix)
  observed_statistic <- sum((obs_curve$msd - sim_mean)^2)
  sim_statistics <- vapply(seq_len(ncol(sim_matrix)), function(i) {
    reference_mean <- rowMeans(sim_matrix[, -i, drop = FALSE])
    sum((sim_matrix[, i] - reference_mean)^2)
  }, numeric(1))

  envelope <- tibble::tibble(
    lag = seq_len(max_lag),
    mean = sim_mean,
    lo = apply(sim_matrix, 1, stats::quantile, probs = envelope_probs[1]),
    hi = apply(sim_matrix, 1, stats::quantile, probs = envelope_probs[2])
  )

  simulated_curves <- tibble::tibble(
    sim_id = rep(names(sims), each = max_lag),
    lag = rep(seq_len(max_lag), times = length(sims)),
    msd = as.vector(sim_matrix)
  )

  out <- list(
    metric = "msd",
    max_lag = max_lag,
    observed_curve = obs_curve,
    simulated_curves = simulated_curves,
    envelope = envelope,
    simulated_statistics = tibble::tibble(
      sim_id = names(sims),
      statistic = unname(sim_statistics)
    ),
    statistic = observed_statistic,
    rank_test = mc_rank_test(observed_statistic, sim_statistics, alternative = "greater")
  )
  class(out) <- c("amt_generative_metric_msd", "list")
  out
}

validate_generative_straightness <- function(observed, simulated, ...) {
  obs <- generative_track(observed)
  sims <- generative_simulations(simulated, ...)
  check_generative_simulation_count(sims, minimum = 2)

  observed_value <- generative_straightness_index(obs)
  simulated_values <- vapply(sims, generative_straightness_index, numeric(1))

  if (!is.finite(observed_value) || any(!is.finite(simulated_values))) {
    stop("Straightness index is undefined for at least one zero-length path.", call. = FALSE)
  }

  observed_statistic <- abs(observed_value - mean(simulated_values, na.rm = TRUE))
  sim_statistics <- vapply(seq_along(simulated_values), function(i) {
    abs(simulated_values[i] - mean(simulated_values[-i], na.rm = TRUE))
  }, numeric(1))

  out <- list(
    metric = "straightness",
    measure = "straightness_index",
    observed_value = observed_value,
    simulated_values = tibble::tibble(
      sim_id = names(sims),
      value = unname(simulated_values)
    ),
    statistic = observed_statistic,
    simulated_statistics = tibble::tibble(
      sim_id = names(sims),
      statistic = unname(sim_statistics)
    ),
    rank_test = mc_rank_test(observed_statistic, sim_statistics, alternative = "greater")
  )
  class(out) <- c("amt_generative_metric_straightness", "list")
  out
}

validate_generative_barrier <- function(observed, simulated, barrier,
                                        alternative = c("less", "greater", "two.sided"),
                                        ...) {
  alternative <- match.arg(alternative)

  if (missing(barrier) || is.null(barrier)) {
    stop("`barrier` is required for barrier crossing validation.", call. = FALSE)
  }

  barrier_geom <- validate_generative_barrier_geometry(barrier)
  obs <- generative_track(observed, barrier_crs = sf::st_crs(barrier_geom))
  sims <- generative_simulations(simulated, barrier_crs = sf::st_crs(barrier_geom), ...)
  check_generative_simulation_count(sims, minimum = 2)

  observed_count <- count_generative_barrier_crossings(obs, barrier_geom)
  simulated_counts <- vapply(sims, count_generative_barrier_crossings, integer(1),
                             barrier_geom = barrier_geom)

  out <- list(
    metric = "barrier",
    method = "segment_intersection_count",
    observed_count = observed_count,
    simulated_counts = tibble::tibble(
      sim_id = names(sims),
      count = unname(simulated_counts)
    ),
    rank_test = mc_rank_test(observed_count, simulated_counts, alternative = alternative),
    barrier = barrier
  )
  class(out) <- c("amt_generative_metric_barrier", "list")
  out
}

