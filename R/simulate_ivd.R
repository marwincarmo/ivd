##' Simulate MELSM data with known deviant clusters
##'
##' Generates data from the intercept-only mixed-effects location-scale
##' model that [ivd()] fits, with a known subset of clusters whose residual
##' (within-cluster) SD deviates from the common baseline:
##' \deqn{y_{ij} = \beta_0 + u_j + e_{ij}, \quad
##'       e_{ij} \sim N(0, \sigma_j^2),}
##' where \eqn{u_j \sim N(0, \tau_{loc}^2)} are location random intercepts
##' and \eqn{\sigma_j = \sigma_{within} \cdot f_j} with \eqn{f_j = 1} for
##' regular clusters and \eqn{f_j =} `deviant_factor` for the deviant ones.
##'
##' Because the ground truth is known, the output is suited for power
##' analysis (how many clusters/observations to detect a given SD
##' inflation?), teaching, and validating recovery: fit with
##' \code{ivd(y ~ 1 + (1 | id), ~ 1 + (1 | id), data = sim$data, ...)} and
##' compare \code{pip(fit)} against \code{sim$truth}.
##' @title Simulate data with deviant within-cluster variances
##' @param J Number of clusters. Defaults to 100.
##' @param n_j Observations per cluster: a single number or a length-`J`
##'   vector. Defaults to 20.
##' @param n_deviant Number of deviant clusters, drawn at random. Defaults
##'   to 5% of `J` (rounded up). Ignored when `deviant_clusters` is given.
##' @param deviant_clusters Optional integer vector of cluster indices to
##'   make deviant (overrides `n_deviant`).
##' @param deviant_factor Multiplicative factor on the within-cluster SD of
##'   the deviant clusters: 2 doubles it, 0.5 halves it. A single number or
##'   one value per deviant cluster. Defaults to 2.
##' @param beta Location grand mean (fixed intercept). Defaults to 0.
##' @param sigma_within Baseline within-cluster SD of the regular clusters.
##'   Defaults to 1. With `family = "student"` this is the *scale* of the t
##'   (its SD is `sigma_within * sqrt(df/(df-2))`).
##' @param tau_loc SD of the location random intercepts. Defaults to 0.5.
##' @param family Residual distribution: `"gaussian"` (default) or
##'   `"student"` for heavy-tailed t residuals with `df` degrees of freedom
##'   -- useful for checking how heavy tails masquerade as variance
##'   heterogeneity under a gaussian [ivd()] fit.
##' @param df Degrees of freedom of the student-t residuals (must be > 2 so
##'   the variance exists). Only used with `family = "student"`. Defaults
##'   to 3.
##' @param seed Optional seed for reproducibility.
##' @return A list of class `ivd_sim`:
##' \itemize{
##'   \item \code{data}: data frame with `y` and cluster `id`, ready for
##'         [ivd()].
##'   \item \code{truth}: data frame with one row per cluster: `id`,
##'         `u_loc`, `sd_within`, `factor`, and `deviant`.
##'   \item \code{params}: the simulation parameters.
##' }
##' @examples
##' sim <- simulate_ivd(J = 20, n_j = 10, n_deviant = 2, seed = 1)
##' sim
##' subset(sim$truth, deviant)
##' \dontrun{
##' fit <- ivd(y ~ 1 + (1 | id), ~ 1 + (1 | id), data = sim$data,
##'            niter = 2000, nburnin = 2000)
##' merge(pip(fit), sim$truth, by.x = "cluster_id", by.y = "id")
##' }
##' @author Philippe Rast
##' @export
simulate_ivd <- function(J = 100, n_j = 20,
                         n_deviant = ceiling(0.05 * J),
                         deviant_clusters = NULL,
                         deviant_factor = 2,
                         beta = 0, sigma_within = 1, tau_loc = 0.5,
                         family = c("gaussian", "student"), df = 3,
                         seed = NULL) {
  family <- match.arg(family)
  if (!is.numeric(df) || length(df) != 1 || df <= 2) {
    stop("`df` must be a single number > 2 (so the residual variance exists).",
         call. = FALSE)
  }
  if (!is.numeric(J) || length(J) != 1 || J < 2 || J != round(J)) {
    stop("`J` must be a single integer >= 2.", call. = FALSE)
  }
  if (!is.numeric(n_j) || !length(n_j) %in% c(1, J) ||
      any(n_j < 1) || any(n_j != round(n_j))) {
    stop("`n_j` must be a positive integer, or a length-J vector of them.",
         call. = FALSE)
  }
  if (!is.numeric(sigma_within) || length(sigma_within) != 1 || sigma_within <= 0) {
    stop("`sigma_within` must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(tau_loc) || length(tau_loc) != 1 || tau_loc < 0) {
    stop("`tau_loc` must be a single non-negative number.", call. = FALSE)
  }
  if (!is.numeric(beta) || length(beta) != 1) {
    stop("`beta` must be a single number (the location grand mean).",
         call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)

  if (is.null(deviant_clusters)) {
    if (!is.numeric(n_deviant) || length(n_deviant) != 1 ||
        n_deviant < 0 || n_deviant > J || n_deviant != round(n_deviant)) {
      stop("`n_deviant` must be a single integer in 0..J.", call. = FALSE)
    }
    deviant_clusters <- sort(sample.int(J, n_deviant))
  } else {
    if (!is.numeric(deviant_clusters) ||
        any(deviant_clusters < 1) || any(deviant_clusters > J) ||
        any(deviant_clusters != round(deviant_clusters)) ||
        anyDuplicated(deviant_clusters)) {
      stop("`deviant_clusters` must be distinct cluster indices in 1..J.",
           call. = FALSE)
    }
    deviant_clusters <- sort(as.integer(deviant_clusters))
  }
  n_dev <- length(deviant_clusters)
  if (!is.numeric(deviant_factor) || any(deviant_factor <= 0) ||
      !length(deviant_factor) %in% c(1, max(n_dev, 1))) {
    stop("`deviant_factor` must be positive: a single value or one per ",
         "deviant cluster.", call. = FALSE)
  }

  n_j <- rep_len(as.integer(n_j), J)
  factor_j <- rep(1, J)
  factor_j[deviant_clusters] <- rep_len(deviant_factor, n_dev)
  sd_within <- sigma_within * factor_j
  u_loc <- rnorm(J, 0, tau_loc)

  id <- rep(seq_len(J), times = n_j)
  resid <- if (family == "student") {
    sd_within[id] * stats::rt(sum(n_j), df = df)
  } else {
    rnorm(sum(n_j), 0, sd_within[id])
  }
  y <- beta + u_loc[id] + resid

  truth <- data.frame(
    id = seq_len(J),
    n = n_j,
    u_loc = u_loc,
    sd_within = sd_within,
    factor = factor_j,
    deviant = factor_j != 1
  )

  structure(
    list(
      data = data.frame(y = y, id = id),
      truth = truth,
      params = list(J = J, n_j = n_j, deviant_clusters = deviant_clusters,
                    deviant_factor = deviant_factor, beta = beta,
                    sigma_within = sigma_within, tau_loc = tau_loc,
                    family = family, df = if (family == "student") df else NULL,
                    seed = seed)
    ),
    class = c("ivd_sim", "list")
  )
}

##' Print method for simulated ivd data
##' @title Print an ivd_sim object
##' @param x An object of class `ivd_sim`.
##' @param ... Not used.
##' @return `x`, invisibly.
##' @author Philippe Rast
##' @export
print.ivd_sim <- function(x, ...) {
  p <- x$params
  dev <- x$truth[x$truth$deviant, , drop = FALSE]
  cat(sprintf("Simulated MELSM data: %d clusters, %d observations%s\n",
              p$J, nrow(x$data),
              if (identical(p$family, "student")) {
                sprintf(" (student-t residuals, df = %g)", p$df)
              } else ""))
  cat(sprintf("  Baseline within-cluster %s: %g (location grand mean %g, tau_loc %g)\n",
              if (identical(p$family, "student")) "scale" else "SD",
              p$sigma_within, p$beta, p$tau_loc))
  if (nrow(dev) == 0) {
    cat("  No deviant clusters.\n")
  } else {
    cat(sprintf("  %d deviant cluster(s): %s (SD factor %s)\n",
                nrow(dev), paste(dev$id, collapse = ", "),
                paste(unique(dev$factor), collapse = ", ")))
  }
  cat("\nFit with, e.g.:\n  ivd(y ~ 1 + (1 | id), ~ 1 + (1 | id), data = <sim>$data, ...)\n")
  invisible(x)
}
