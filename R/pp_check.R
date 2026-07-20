##' Observed and replicated within-cluster SDs
##'
##' Internal core of [pp_check.ivd()]: reconstructs \code{mu} and \code{tau}
##' per posterior draw from the stored design matrices, simulates replicated
##' data \code{Y_rep ~ N(mu, tau)}, and computes the within-cluster SD of
##' every replicate. Pure R (no NIMBLE), so it is unit-testable under
##' coverage.
##' @param object An object of class `ivd`.
##' @param ndraws Number of posterior draws to simulate from (sampled
##'   without replacement from the pooled chains; capped at the number of
##'   available draws).
##' @param seed Optional seed for reproducibility of the draw sampling and
##'   the replicate noise.
##' @return List with `observed` (named length-J vector of within-cluster
##'   SDs), `rep` (an `ndraws` x J matrix of replicated SDs), and `n_obs`
##'   (cluster sizes). Clusters with fewer than 2 observations are dropped
##'   (their SD is undefined).
##' @keywords internal
##' @importFrom stats quantile rnorm
.pp_cluster_sds <- function(object, ndraws = 100, seed = NULL) {
  if (is.null(object$X) || is.null(object$Z) ||
      is.null(object$X_scale) || is.null(object$Z_scale) || is.null(object$Y)) {
    stop("Cannot run pp_check: design matrices are missing from the ivd ",
         "object. Refit with the current version of ivd().", call. = FALSE)
  }
  Kr <- object$nimble_constants$Kr
  Sr <- object$nimble_constants$Sr
  P <- Kr + Sr
  J <- object$nimble_constants$J
  group_id <- object$Y$group_id
  Y <- object$Y$Y
  N <- length(Y)

  draws <- .pooled_draws(object)
  cn <- colnames(draws)

  beta_cols <- grep("^beta\\[", cn, value = TRUE)
  beta_cols <- beta_cols[order(as.integer(gsub("\\D", "", beta_cols)))]
  zeta_cols <- grep("^zeta\\[", cn, value = TRUE)
  zeta_cols <- zeta_cols[order(as.integer(gsub("\\D", "", zeta_cols)))]
  ## u[j, p] as a J x P matrix of column names for direct reshaping
  u_names <- outer(seq_len(J), seq_len(P),
                   function(j, p) paste0("u[", j, ", ", p, "]"))
  if (!all(c(beta_cols, zeta_cols, u_names) %in% cn) ||
      length(beta_cols) != ncol(object$X) ||
      length(zeta_cols) != ncol(object$X_scale)) {
    stop("Cannot run pp_check: expected beta/zeta/u columns are missing ",
         "from the stored samples.", call. = FALSE)
  }

  ## student-t fits replicate from a t with the draw's df; tau is its scale
  student <- identical(object$family, "student")
  if (student && !("nu" %in% cn)) {
    stop("Cannot run pp_check: this student-t fit has no stored `nu` draws.",
         call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)
  ndraws <- min(ndraws, nrow(draws))
  draw_ids <- sample.int(nrow(draws), ndraws)

  ## observed statistic; clusters of size 1 have no SD
  observed <- c(tapply(Y, group_id, sd)) # c() drops tapply's array dim
  n_obs <- as.integer(table(factor(group_id, levels = seq_len(J))))
  keep <- n_obs >= 2

  rep_mat <- matrix(NA_real_, nrow = ndraws, ncol = J)
  Xb_cols <- seq_len(Kr)
  Zs_cols <- Kr + seq_len(Sr)
  for (s in seq_len(ndraws)) {
    d <- draws[draw_ids[s], ]
    u_mat <- matrix(d[u_names], nrow = J, ncol = P)
    mu <- as.numeric(object$X %*% d[beta_cols]) +
      rowSums(object$Z * u_mat[group_id, Xb_cols, drop = FALSE])
    tau <- exp(as.numeric(object$X_scale %*% d[zeta_cols]) +
                 rowSums(object$Z_scale * u_mat[group_id, Zs_cols, drop = FALSE]))
    y_rep <- if (student) {
      mu + tau * stats::rt(N, df = d[["nu"]])
    } else {
      rnorm(N, mu, tau)
    }
    rep_mat[s, ] <- tapply(y_rep, group_id, sd)
  }

  list(observed = observed[keep],
       rep = rep_mat[, keep, drop = FALSE],
       n_obs = n_obs[keep],
       cluster_index = which(keep),
       dropped = which(!keep))
}

##' Posterior predictive check
##'
##' Generic; see [pp_check.ivd()]. Defined in-package so that `ivd` does not
##' depend on \pkg{bayesplot}; if bayesplot is attached its same-named
##' generic dispatches to the ivd method just the same.
##' @title Posterior predictive check
##' @param object A fitted model object.
##' @param ... Passed to methods.
##' @return See the method documentation.
##' @author Philippe Rast
##' @export
pp_check <- function(object, ...) UseMethod("pp_check")

##' Posterior predictive check of the within-cluster SDs
##'
##' The central claim of an `ivd` fit concerns the residual (within-cluster)
##' variability, so the check replicates data from the posterior and
##' compares each cluster's *observed* within-cluster SD with its posterior
##' predictive distribution:
##' \itemize{
##'   \item \code{type = "intervals"} (default): per cluster, the observed
##'         SD against the inner/outer predictive intervals, ordered by
##'         observed SD. Clusters falling outside the outer interval are
##'         highlighted and labelled -- systematic misfit there suggests the
##'         normal likelihood is inadequate (e.g. heavy tails).
##'   \item \code{type = "density"}: the observed density of cluster SDs
##'         overlaid on the replicated densities (one thin line per draw).
##' }
##' Clusters with a single observation are dropped (their SD is undefined).
##' @title Posterior predictive check for ivd models
##' @param object An object of class `ivd`.
##' @param ndraws Number of posterior draws to replicate data from.
##'   Defaults to 100 (capped at the number of stored draws).
##' @param type `"intervals"` (default) or `"density"`; see Details.
##' @param probs Inner and outer probability of the predictive intervals
##'   for `type = "intervals"`. Defaults to `c(0.5, 0.9)`.
##' @param labels Cluster labels for the highlighted points: `"index"`
##'   (default) uses the internal 1..J index; `"original"` the user's own
##'   grouping IDs. Matches [summary.ivd()] and [plot.ivd()].
##' @param seed Optional seed making the replicate noise reproducible.
##' @param ... Not used.
##' @return A `ggplot` object.
##' @examples
##' \dontrun{
##' pp_check(fit)
##' pp_check(fit, type = "density", ndraws = 50)
##' }
##' @author Philippe Rast
##' @export
pp_check.ivd <- function(object, ndraws = 100, type = c("intervals", "density"),
                         probs = c(0.5, 0.9), labels = c("index", "original"),
                         seed = NULL, ...) {
  type <- match.arg(type)
  labels <- match.arg(labels)
  if (identical(labels, "original") && is.null(object$group_labels)) {
    warning("This ivd object predates 'group_labels'; ",
            "points keep the internal index.")
    labels <- "index"
  }
  if (!is.numeric(probs) || length(probs) != 2 || anyNA(probs) ||
      any(probs <= 0) || any(probs >= 1) || probs[1] >= probs[2]) {
    stop("`probs` must be two increasing values in (0, 1), e.g. c(0.5, 0.9).",
         call. = FALSE)
  }

  stats <- .pp_cluster_sds(object, ndraws = ndraws, seed = seed)
  if (length(stats$dropped)) {
    warning(length(stats$dropped),
            " cluster(s) with a single observation were dropped from the check.")
  }

  cluster_label <- if (identical(labels, "original")) {
    object$group_labels[stats$cluster_index]
  } else {
    as.character(stats$cluster_index)
  }

  if (type == "intervals") {
    qs <- apply(stats$rep, 2, quantile,
                probs = c((1 - probs[2]) / 2, (1 - probs[1]) / 2, 0.5,
                          1 - (1 - probs[1]) / 2, 1 - (1 - probs[2]) / 2),
                names = FALSE)
    df <- data.frame(
      cluster_index = stats$cluster_index,
      label = cluster_label,
      obs = as.numeric(stats$observed),
      lo_outer = qs[1, ], lo_inner = qs[2, ], med = qs[3, ],
      hi_inner = qs[4, ], hi_outer = qs[5, ]
    )
    df$outside <- df$obs < df$lo_outer | df$obs > df$hi_outer
    df <- df[order(df$obs), ]
    df$ordered <- seq_len(nrow(df))

    plt <- ggplot(df, aes(x = ordered)) +
      geom_linerange(aes(ymin = lo_outer, ymax = hi_outer),
                     color = "grey65", linewidth = 0.4) +
      geom_linerange(aes(ymin = lo_inner, ymax = hi_inner),
                     color = "grey35", linewidth = 1) +
      geom_point(aes(y = med), shape = 21, fill = "white",
                 color = "grey35", size = 1.6) +
      geom_point(data = subset(df, !outside), aes(y = obs),
                 color = "#0265a5", size = 1.8) +
      geom_point(data = subset(df, outside), aes(y = obs),
                 color = "#B2182B", size = 2.2) +
      labs(
        x = "Clusters (ordered by observed SD)",
        y = "Within-cluster SD",
        title = "Posterior predictive check: within-cluster SD",
        subtitle = sprintf(
          "Observed (points) vs %d%%/%d%% predictive intervals from %d draws; %d cluster(s) outside",
          round(probs[1] * 100), round(probs[2] * 100), nrow(stats$rep),
          sum(df$outside))
      )
    if (any(df$outside)) {
      .require_suggest("ggrepel", "`geom_text_repel()`")
      plt <- plt + ggrepel::geom_text_repel(
        data = subset(df, outside), aes(y = obs, label = label),
        size = 3, color = "#B2182B", point.padding = 0.4
      )
    }
    return(plt)
  }

  ## type == "density": replicated cluster-SD densities vs the observed one
  df_rep <- data.frame(
    rep_id = rep(seq_len(nrow(stats$rep)), times = ncol(stats$rep)),
    value = as.numeric(stats$rep)
  )
  df_obs <- data.frame(value = as.numeric(stats$observed))
  ggplot(df_rep, aes(x = value)) +
    geom_line(aes(group = rep_id), stat = "density",
              color = "grey70", alpha = 0.3, linewidth = 0.3) +
    geom_line(data = df_obs, stat = "density",
              color = "#0265a5", linewidth = 1.2) +
    labs(
      x = "Within-cluster SD",
      y = "Density",
      title = "Posterior predictive check: within-cluster SD",
      subtitle = sprintf("Observed (blue) vs %d replicated densities",
                         nrow(stats$rep))
    )
}
