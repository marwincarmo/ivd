##' Monte Carlo diagnostics for the posterior inclusion probabilities
##'
##' A PIP is the posterior mean of a binary spike-and-slab indicator, so the
##' usual continuous-parameter diagnostics do not transfer: split-Rhat on a
##' 0/1 chain is misleading, and a seemingly precise PIP can still differ
##' markedly between chains when the indicator mixes slowly. This function
##' reports, per cluster and random scale effect:
##' \itemize{
##'   \item the per-chain PIPs (columns `chain1`, `chain2`, ...),
##'   \item `mcse`: the Monte Carlo standard error of the pooled PIP,
##'         estimated from the between-chain spread
##'         (\eqn{\mathrm{sd}(\text{chain PIPs})/\sqrt{C}}),
##'   \item `range`: the largest between-chain difference, and
##'   \item `unstable`: whether the chains disagree on the classification at
##'         `pip_level` (some chains put the cluster at or above the
##'         threshold, others below it).
##' }
##' Clusters flagged `unstable` should not be classified either way from
##' this fit -- run more iterations (or more chains) first.
##' @title Monte Carlo error and chain agreement of the PIPs
##' @param object An object of class `ivd`.
##' @param pip_level PIP threshold whose classification is checked for
##'   between-chain agreement. Defaults to 0.75, matching [plot.ivd()].
##' @return A data frame of class `pip_diagnostics`: one row per cluster x
##'   random scale effect with `scale_var`, `cluster_index`, `cluster_id`,
##'   `pip`, the per-chain PIP columns, `mcse`, `range`, and `unstable`.
##'   `pip_level` is stored as an attribute. Its `print()` method shows a
##'   compact summary and the unstable clusters, if any.
##' @examples
##' \dontrun{
##' diag <- pip_diagnostics(fit)
##' diag
##' subset(diag, unstable)
##' }
##' @author Philippe Rast
##' @export
pip_diagnostics <- function(object, pip_level = 0.75) {
  if (!inherits(object, "ivd")) {
    stop("`object` must be a fitted ivd model.", call. = FALSE)
  }
  if (!is.numeric(pip_level) || length(pip_level) != 1 ||
      is.na(pip_level) || pip_level <= 0 || pip_level >= 1) {
    stop("`pip_level` must be a single number in (0, 1).", call. = FALSE)
  }

  Kr <- object$nimble_constants$Kr
  Sr <- object$nimble_constants$Sr
  J <- object$nimble_constants$J
  scale_vars <- colnames(object$Z_scale)
  labels <- if (!is.null(object$group_labels)) {
    object$group_labels
  } else {
    as.character(seq_len(J))
  }
  labels <- .maybe_numeric(labels)

  chains <- lapply(object$samples, function(chain) chain$samples)
  C <- length(chains)

  out <- do.call(rbind, lapply(seq_len(Sr), function(s) {
    p <- Kr + s
    ss_cols <- paste0("ss[", p, ", ", seq_len(J), "]")
    if (!all(ss_cols %in% colnames(chains[[1]]))) {
      stop("Could not find the ss columns for scale random effect ",
           scale_vars[s], " in the stored samples.")
    }
    ## J x C matrix of per-chain PIPs
    per_chain <- vapply(chains, function(ch) {
      unname(colMeans(ch[, ss_cols, drop = FALSE]))
    }, numeric(J))

    d <- data.frame(
      scale_var = scale_vars[s],
      cluster_index = seq_len(J),
      cluster_id = labels,
      pip = rowMeans(per_chain),
      row.names = NULL
    )
    chain_df <- as.data.frame(per_chain)
    names(chain_df) <- paste0("chain", seq_len(C))
    d <- cbind(d, chain_df)
    d$mcse <- apply(per_chain, 1, stats::sd) / sqrt(C)
    d$range <- apply(per_chain, 1, function(x) max(x) - min(x))
    ## chains disagree about the classification at pip_level
    d$unstable <- apply(per_chain, 1, function(x) {
      any(x >= pip_level) && any(x < pip_level)
    })
    d
  }))

  attr(out, "pip_level") <- pip_level
  class(out) <- c("pip_diagnostics", "data.frame")
  out
}

##' Print method for `pip_diagnostics` objects
##' @title Print PIP Monte Carlo diagnostics
##' @param x An object of class `pip_diagnostics`.
##' @param ... Not used.
##' @return `x`, invisibly.
##' @author Philippe Rast
##' @export
print.pip_diagnostics <- function(x, ...) {
  pip_level <- attr(x, "pip_level")
  chain_cols <- grep("^chain\\d+$", names(x), value = TRUE)
  cat(sprintf("PIP Monte Carlo diagnostics (%d chains, classification at PIP >= %s)\n\n",
              length(chain_cols), format(pip_level)))
  cat(sprintf("  Max. MCSE of a PIP:              %.3f\n", max(x$mcse)))
  cat(sprintf("  Max. between-chain PIP range:    %.3f\n", max(x$range)))

  n_unstable <- sum(x$unstable)
  if (n_unstable == 0) {
    cat(sprintf("  All %d PIPs are classified consistently across chains.\n",
                nrow(x)))
  } else {
    cat(sprintf(paste0("  %d of %d PIPs are classified inconsistently across ",
                       "chains (see below);\n  consider more iterations before ",
                       "classifying these clusters.\n\n"),
                n_unstable, nrow(x)))
    d <- x[x$unstable, c("scale_var", "cluster_id", "pip", chain_cols, "mcse", "range")]
    print.data.frame(d, row.names = FALSE, digits = 3)
  }
  invisible(x)
}
