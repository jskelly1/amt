`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

generative_track <- function(x, barrier_crs = NULL) {
  if (inherits(x, "sf")) {
    return(generative_track_sf(x, barrier_crs = barrier_crs))
  }

  if (inherits(x, "steps_xyt")) {
    return(generative_track(as_track(x), barrier_crs = barrier_crs))
  }

  if (!is.data.frame(x)) {
    stop("Track inputs must be data frames, sf point objects, or amt track/step objects.",
         call. = FALSE)
  }

  out <- tibble::as_tibble(x)
  out <- standardize_generative_track_columns(out)
  validate_generative_track(out)
}

generative_track_sf <- function(x, barrier_crs = NULL) {
  geom <- sf::st_geometry(x)
  geom_type <- unique(as.character(sf::st_geometry_type(geom)))

  if (!all(geom_type %in% "POINT")) {
    stop("`sf` tracks must contain only POINT geometries.", call. = FALSE)
  }
  if (isTRUE(sf::st_is_longlat(geom))) {
    stop("`sf` tracks must use projected planar coordinates; longitude/latitude CRS are not supported.",
         call. = FALSE)
  }
  if (!is.null(barrier_crs) && !is.na(sf::st_crs(geom)) &&
      !is.na(barrier_crs) && sf::st_crs(geom) != barrier_crs) {
    stop("Track and `barrier` must use the same CRS.", call. = FALSE)
  }

  coords <- sf::st_coordinates(geom)
  out <- tibble::as_tibble(sf::st_drop_geometry(x))
  out$x_ <- coords[, 1]
  out$y_ <- coords[, 2]
  validate_generative_track(out)
}

standardize_generative_track_columns <- function(out) {
  if (all(c("x_", "y_") %in% names(out))) {
    return(out)
  }
  if (all(c("x1_", "y1_", "x2_", "y2_") %in% names(out))) {
    return(tibble::tibble(
      x_ = c(out$x1_[1], out$x2_),
      y_ = c(out$y1_[1], out$y2_)
    ))
  }

  stop("Could not find coordinate columns. Use `x_`/`y_` or amt step columns.",
       call. = FALSE)
}

validate_generative_track <- function(out) {
  if (!is.numeric(out$x_) || !is.numeric(out$y_)) {
    stop("Coordinate columns `x_` and `y_` must be numeric.", call. = FALSE)
  }

  keep <- is.finite(out$x_) & is.finite(out$y_)
  if (!all(keep)) {
    warning("Removed rows with missing or non-finite coordinates.", call. = FALSE)
    out <- out[keep, , drop = FALSE]
  }
  if (nrow(out) < 2) {
    stop("A track must contain at least two finite locations.", call. = FALSE)
  }

  tibble::as_tibble(out)
}

generative_simulations <- function(simulated, id_col = NULL, barrier_crs = NULL, ...) {
  if (inherits(simulated, "amt_generative_simulations")) {
    return(simulated)
  }

  if (is.list(simulated) && !is.data.frame(simulated) && !inherits(simulated, "sf")) {
    sims <- lapply(simulated, generative_track, barrier_crs = barrier_crs)
  } else if (is.data.frame(simulated)) {
    sim_tbl <- tibble::as_tibble(simulated)
    id_col <- id_col %||% detect_generative_simulation_id(sim_tbl)

    if (is.null(id_col)) {
      stop("`simulated` is a data frame, so it must contain a simulation identifier column.",
           call. = FALSE)
    }
    if (!is.character(id_col) || length(id_col) != 1L || is.na(id_col) ||
        !id_col %in% names(sim_tbl)) {
      stop("`id_col` must name one simulation identifier column in `simulated`.",
           call. = FALSE)
    }

    sims <- lapply(split(sim_tbl, sim_tbl[[id_col]]), generative_track,
                   barrier_crs = barrier_crs)
  } else {
    stop("`simulated` must be a list of tracks or a data frame of simulated tracks.",
         call. = FALSE)
  }

  if (length(sims) < 1) {
    stop("At least one simulated track is required.", call. = FALSE)
  }
  if (is.null(names(sims)) || any(names(sims) == "")) {
    names(sims) <- paste0("sim_", seq_along(sims))
  }

  class(sims) <- c("amt_generative_simulations", "list")
  sims
}

detect_generative_simulation_id <- function(x) {
  candidates <- c("sim_id", "sim_id_", ".simulation", "simulation",
                  ".replicate", "replicate")
  matched <- candidates[candidates %in% names(x)]
  if (length(matched) < 1L) {
    return(NULL)
  }
  matched[1]
}

generative_straightness_index <- function(track) {
  trk <- generative_track(track)
  dx <- diff(trk$x_)
  dy <- diff(trk$y_)
  path_length <- sum(sqrt(dx^2 + dy^2))

  if (!is.finite(path_length) || path_length <= 0) {
    return(NA_real_)
  }

  net_displacement <- sqrt(
    (trk$x_[nrow(trk)] - trk$x_[1])^2 +
      (trk$y_[nrow(trk)] - trk$y_[1])^2
  )

  net_displacement / path_length
}

generative_msd_curve <- function(track, max_lag = NULL) {
  trk <- generative_track(track)
  n <- nrow(trk)

  if (is.null(max_lag)) {
    max_lag <- n - 1L
  }
  max_lag <- min(as.integer(max_lag), n - 1L)

  if (!is.finite(max_lag) || max_lag < 1L) {
    stop("`max_lag` must be at least 1 and smaller than the track length.",
         call. = FALSE)
  }

  values <- vapply(seq_len(max_lag), function(lag_i) {
    dx <- trk$x_[(lag_i + 1L):n] - trk$x_[seq_len(n - lag_i)]
    dy <- trk$y_[(lag_i + 1L):n] - trk$y_[seq_len(n - lag_i)]
    mean(dx^2 + dy^2)
  }, numeric(1))

  tibble::tibble(lag = seq_len(max_lag), msd = values)
}

generative_track_summary <- function(x) {
  trk <- generative_track(x)
  tibble::tibble(
    n_locations = nrow(trk),
    x_min = min(trk$x_),
    x_max = max(trk$x_),
    y_min = min(trk$y_),
    y_max = max(trk$y_),
    straightness = generative_straightness_index(trk)
  )
}

generative_simulation_summary <- function(sims) {
  n_locations <- vapply(sims, nrow, integer(1))
  tibble::tibble(
    n_simulations = length(sims),
    min_locations = min(n_locations),
    median_locations = stats::median(n_locations),
    max_locations = max(n_locations)
  )
}

check_generative_simulation_count <- function(sims, minimum = 2) {
  if (length(sims) < minimum) {
    stop("At least ", minimum, " simulated tracks are required for this validation.",
         call. = FALSE)
  }
}

normalize_generative_grid_size <- function(grid_size) {
  if (!is.numeric(grid_size) || !length(grid_size) %in% c(1L, 2L)) {
    stop("`grid_size` must be an integer scalar or length-two integer vector.",
         call. = FALSE)
  }
  if (any(!is.finite(grid_size)) || any(grid_size != floor(grid_size))) {
    stop("`grid_size` values must be whole numbers.", call. = FALSE)
  }

  grid_size <- as.integer(grid_size)
  if (length(grid_size) == 1L) {
    grid_size <- rep(grid_size, 2)
  }
  if (any(grid_size < 2L)) {
    stop("`grid_size` values must be at least 2.", call. = FALSE)
  }

  grid_size
}

validate_generative_bounds <- function(bounds) {
  required <- c("xmin", "xmax", "ymin", "ymax")
  if (!is.numeric(bounds) || !all(required %in% names(bounds))) {
    stop("`bounds` must be a named numeric vector with xmin, xmax, ymin, and ymax.",
         call. = FALSE)
  }

  bounds <- bounds[required]
  if (any(!is.finite(bounds))) {
    stop("All `bounds` values must be finite.", call. = FALSE)
  }
  if (bounds["xmin"] >= bounds["xmax"] || bounds["ymin"] >= bounds["ymax"]) {
    stop("`bounds` must satisfy xmin < xmax and ymin < ymax.", call. = FALSE)
  }

  bounds
}

combined_generative_bounds <- function(obs, sims) {
  all_x <- c(obs$x_, unlist(lapply(sims, `[[`, "x_"), use.names = FALSE))
  all_y <- c(obs$y_, unlist(lapply(sims, `[[`, "y_"), use.names = FALSE))
  bounds <- c(
    xmin = min(all_x),
    xmax = max(all_x),
    ymin = min(all_y),
    ymax = max(all_y)
  )

  if (bounds["xmin"] == bounds["xmax"]) {
    bounds[c("xmin", "xmax")] <- bounds[c("xmin", "xmax")] + c(-0.5, 0.5)
  }
  if (bounds["ymin"] == bounds["ymax"]) {
    bounds[c("ymin", "ymax")] <- bounds[c("ymin", "ymax")] + c(-0.5, 0.5)
  }

  bounds
}

check_generative_tracks_within_bounds <- function(obs, sims, bounds) {
  tracks <- c(list(observed = obs), sims)
  outside <- vapply(tracks, function(track) {
    trk <- generative_track(track)
    any(
      trk$x_ < bounds["xmin"] |
        trk$x_ > bounds["xmax"] |
        trk$y_ < bounds["ymin"] |
        trk$y_ > bounds["ymax"]
    )
  }, logical(1))

  if (any(outside)) {
    stop("`bounds` must contain all observed and simulated track coordinates.",
         call. = FALSE)
  }

  invisible(TRUE)
}

empirical_grid_ud <- function(track, grid_size, bounds) {
  trk <- generative_track(track)
  x_breaks <- seq(bounds["xmin"], bounds["xmax"], length.out = grid_size[1] + 1L)
  y_breaks <- seq(bounds["ymin"], bounds["ymax"], length.out = grid_size[2] + 1L)

  x_bin <- pmin(pmax(findInterval(trk$x_, x_breaks, all.inside = TRUE), 1L),
                grid_size[1])
  y_bin <- pmin(pmax(findInterval(trk$y_, y_breaks, all.inside = TRUE), 1L),
                grid_size[2])

  cells <- tibble::tibble(x_bin = x_bin, y_bin = y_bin)
  counts <- as.data.frame(table(cells$x_bin, cells$y_bin), stringsAsFactors = FALSE)
  names(counts) <- c("x_bin", "y_bin", "n")
  counts$x_bin <- as.integer(as.character(counts$x_bin))
  counts$y_bin <- as.integer(as.character(counts$y_bin))
  counts <- counts[counts$n > 0, , drop = FALSE]

  if (nrow(counts) < 1L || sum(counts$n) <= 0) {
    stop("Empirical grid UD has no occupied cells.", call. = FALSE)
  }

  tibble::tibble(
    x = x_breaks[counts$x_bin] + diff(x_breaks)[counts$x_bin] / 2,
    y = y_breaks[counts$y_bin] + diff(y_breaks)[counts$y_bin] / 2,
    mass = counts$n / sum(counts$n)
  )
}

wasserstein_ud <- function(ud_a, ud_b) {
  check_empirical_ud(ud_a, "ud_a")
  check_empirical_ud(ud_b, "ud_b")

  a <- transport::wpp(as.matrix(ud_a[, c("x", "y")]), ud_a$mass)
  b <- transport::wpp(as.matrix(ud_b[, c("x", "y")]), ud_b$mass)
  distance <- transport::wasserstein(a, b, p = 1, method = "networkflow")

  if (!is.numeric(distance) || length(distance) != 1L || !is.finite(distance)) {
    stop("Wasserstein distance calculation returned a non-finite value.",
         call. = FALSE)
  }

  unname(distance)
}

check_empirical_ud <- function(ud, name) {
  required <- c("x", "y", "mass")
  if (!is.data.frame(ud) || !all(required %in% names(ud))) {
    stop("`", name, "` must be an empirical grid UD with x, y, and mass columns.",
         call. = FALSE)
  }

  values <- ud[, required]
  if (any(!is.finite(as.matrix(values)))) {
    stop("`", name, "` contains non-finite values.", call. = FALSE)
  }
  if (any(ud$mass < 0) || sum(ud$mass) <= 0) {
    stop("`", name, "` must have non-negative masses with positive total mass.",
         call. = FALSE)
  }
}

validate_generative_max_lag <- function(max_lag, common_max_lag) {
  if (!is.numeric(max_lag) || length(max_lag) != 1L || !is.finite(max_lag)) {
    stop("`max_lag` must be one finite whole number.", call. = FALSE)
  }
  if (max_lag != floor(max_lag)) {
    stop("`max_lag` must be a whole number.", call. = FALSE)
  }

  max_lag <- as.integer(max_lag)
  if (max_lag < 1L || max_lag > common_max_lag) {
    stop("`max_lag` must be at least 1 and no larger than the largest common lag.",
         call. = FALSE)
  }

  max_lag
}

validate_generative_envelope_probs <- function(envelope_probs) {
  if (!is.numeric(envelope_probs) || length(envelope_probs) != 2L) {
    stop("`envelope_probs` must be a numeric vector of length two.", call. = FALSE)
  }
  if (any(!is.finite(envelope_probs)) || any(envelope_probs < 0) ||
      any(envelope_probs > 1)) {
    stop("`envelope_probs` values must be finite probabilities.", call. = FALSE)
  }
  if (envelope_probs[1] >= envelope_probs[2]) {
    stop("`envelope_probs` must be ordered from lower to upper probability.",
         call. = FALSE)
  }

  envelope_probs
}

validate_generative_barrier_geometry <- function(barrier) {
  if (!inherits(barrier, "sf") && !inherits(barrier, "sfc")) {
    stop("`barrier` must be an sf or sfc LINESTRING or MULTILINESTRING object.",
         call. = FALSE)
  }

  barrier_geom <- if (inherits(barrier, "sf")) sf::st_geometry(barrier) else barrier
  geom_type <- unique(as.character(sf::st_geometry_type(barrier_geom)))

  if (!all(geom_type %in% c("LINESTRING", "MULTILINESTRING"))) {
    stop("`barrier` must contain only LINESTRING or MULTILINESTRING geometries.",
         call. = FALSE)
  }
  if (is.na(sf::st_crs(barrier_geom))) {
    warning("`barrier` has no CRS. Assuming it uses the same coordinates as the tracks.",
            call. = FALSE)
  } else if (isTRUE(sf::st_is_longlat(barrier_geom))) {
    stop("`barrier` must use projected planar coordinates; longitude/latitude CRS are not supported.",
         call. = FALSE)
  }

  barrier_geom
}

count_generative_barrier_crossings <- function(track, barrier_geom) {
  trk <- generative_track(track)
  if (nrow(trk) < 2L) {
    return(0L)
  }

  barrier_crs <- sf::st_crs(barrier_geom)
  segment_list <- lapply(seq_len(nrow(trk) - 1L), function(i) {
    sf::st_linestring(matrix(
      c(trk$x_[i], trk$y_[i], trk$x_[i + 1L], trk$y_[i + 1L]),
      ncol = 2,
      byrow = TRUE
    ))
  })

  segments <- sf::st_sfc(segment_list, crs = barrier_crs)
  intersections <- sf::st_intersects(segments, barrier_geom, sparse = FALSE)
  as.integer(sum(rowSums(intersections) > 0))
}

check_finite_generative_vector <- function(x, label) {
  if (!is.numeric(x) || length(x) < 1L || any(!is.finite(x))) {
    stop(label, " must contain finite numeric values.", call. = FALSE)
  }
}

generative_metric_p_value <- function(x) {
  if (is.null(x$rank_test)) {
    return(NA_real_)
  }
  x$rank_test$p_value
}

generative_metric_pillar <- function(metric) {
  unname(c(
    ud = "Emergent space use",
    msd = "Diffusion behavior",
    straightness = "Path structure",
    barrier = "Barrier interactions"
  )[metric])
}

format_generative_p <- function(x) {
  if (is.na(x)) {
    return("NA")
  }
  formatC(x, digits = 3, format = "f")
}

format_generative_number <- function(x) {
  formatC(x, digits = 4, format = "fg")
}
