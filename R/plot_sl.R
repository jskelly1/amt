#' Plot step-length distribution
#'
#' @param x `[fit_clogit|random_steps]` \cr A fitted step selection or random steps.
#' @param upper_quantile `[nummeric(1)=0.99]{0-1}` \cr The quantile until where the distribution should be plotted. Typically this will be `0.95` or `0.99`.
#' @param n `[numeric(1)=1000]{>0}` \cr The number of breaks between `0` and `upper_quantile`.
#' @param plot `[logical(1)=TRUE]` \cr Indicates if a plot should be drawn or not.
#' @param engine `[character(1)="base"]` \cr Plotting engine. `"base"` keeps
#' the existing base graphics behavior. `"ggplot2"` returns a ggplot object and
#' requires the optional `ggplot2` package.
#' @template dots_none
#' @return A plot of the step-length distribution.
#' @export
#' @name plot_sl
#' @examples
#' data(deer)
#'
#' # with random steps
#' deer[1:100, ] |> steps_by_burst() |> random_steps() |> plot_sl()
#' deer[1:100, ] |> steps_by_burst() |> random_steps() |> plot_sl(upper_quantile = 0.5)
#'
plot_sl <- function(x, ...) {
  UseMethod("plot_sl", x)
}

#' @export
#' @rdname plot_sl
plot_sl.fit_clogit <- function(x, n = 1000, upper_quantile = 0.99,
                               plot = TRUE, engine = c("base", "ggplot2"), ...) {
  plot_sl_base(
    x = x,
    n = n,
    upper_quantile = upper_quantile,
    plot = plot,
    engine = engine,
    ...
  )
}

#' @export
#' @rdname plot_sl
plot_sl.random_steps <- function(x, n = 1000, upper_quantile = 0.99,
                                 plot = TRUE, engine = c("base", "ggplot2"), ...) {
  plot_sl_base(
    x = x,
    n = n,
    upper_quantile = upper_quantile,
    plot = plot,
    engine = engine,
    ...
  )
}

plot_sl_base <- function(x, n, upper_quantile, plot, engine = c("base", "ggplot2"), ...) {
  engine <- match.arg(engine)
  xx <- sl_distr_params(x)
  if (sl_distr_name(x) == "gamma") {
    to <- qgamma(upper_quantile, shape = xx$shape, scale = xx$scale)
    xs <- seq(0, to, length.out = n)
    ys <- dgamma(xs, shape = xx$shape, scale = xx$scale)
    dat <- data.frame(sl = xs, d = ys)

    if (engine == "ggplot2") {
      if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("The `ggplot2` package is required when `engine = \"ggplot2\"`.",
             call. = FALSE)
      }
      p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$sl, y = .data$d)) +
        ggplot2::geom_line() +
        ggplot2::labs(x = "Distance", y = "Probability")
      if (plot) {
        print(p)
      }
      return(invisible(p))
    }

    if (plot) {
      plot(xs, ys, type = "l",
           ylab = "Probability",
           xlab = "Distance")
    }
    invisible(dat)
  } else {
    stop ("distr not implemented")
  }
}
