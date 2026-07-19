##' Core of the Bayesian FDR rule: given PIPs, order clusters by decreasing
##' PIP and select the largest set whose estimated false discovery rate --
##' the running mean of (1 - PIP) among the selected -- stays at or below
##' the target (Newton et al., 2004; Mueller et al., 2006).
##' @param pips Numeric vector of posterior inclusion probabilities.
##' @param fdr Target FDR level.
##' @return List with `order` (indices sorting pips decreasingly), `cum_fdr`
##'   (estimated FDR when selecting up to each ordered cluster), and
##'   `selected` (logical, in the same order as `order`).
##' @keywords internal
##' @noRd
.fdr_select <- function(pips, fdr) {
  ord <- order(pips, decreasing = TRUE)
  cum_fdr <- cumsum(1 - pips[ord]) / seq_along(ord)
  k <- which(cum_fdr <= fdr)
  k_star <- if (length(k)) max(k) else 0L
  list(order = ord,
       cum_fdr = cum_fdr,
       selected = seq_along(ord) <= k_star)
}

##' Bayesian FDR decision rule for the PIPs
##'
##' Instead of flagging clusters by an arbitrary fixed PIP cutoff (e.g. the
##' conventional 0.75), this selects the largest set of clusters whose
##' estimated Bayesian false discovery rate stays at or below `fdr`:
##' clusters are ordered by decreasing PIP and the estimated FDR of
##' selecting the top \eqn{k} is the mean of their `1 - PIP` values
##' (Newton et al., 2004; Mueller, Parmigiani & Rice, 2006). The implied
##' PIP threshold therefore adapts to the data: it is low when many
##' clusters have clear-cut PIPs and high when evidence is diffuse.
##'
##' Selection is done separately for each random scale effect.
##' @title Select deviant clusters by controlling the Bayesian FDR
##' @param object An object of class `ivd`.
##' @param fdr Target false discovery rate, in (0, 1). Defaults to 0.05.
##' @return A data frame of class `pip_fdr` with one row per cluster and
##'   random scale effect, ordered by decreasing PIP within effect:
##'   `scale_var`, `cluster_index`, `cluster_id`, `pip`, `local_fdr`
##'   (`1 - pip`), `cum_fdr` (estimated FDR when selecting all clusters up
##'   to this row), and `selected`. The target level is stored in
##'   `attr(x, "fdr")`. The implied per-effect PIP thresholds -- usable as
##'   `pip_level` in [plot.ivd()] -- are in `attr(x, "thresholds")`.
##' @references
##' Newton, M. A., Noueiry, A., Sarkar, D., & Ahlquist, P. (2004).
##' Detecting differential gene expression with a semiparametric
##' hierarchical mixture method. *Biostatistics*, 5(2), 155-176.
##'
##' Mueller, P., Parmigiani, G., & Rice, K. (2006). FDR and Bayesian
##' multiple comparisons rules. In *Bayesian Statistics 8*. Oxford
##' University Press.
##' @examples
##' \dontrun{
##' sel <- pip_fdr(fit, fdr = 0.05)
##' sel
##' subset(sel, selected)
##' plot(fit, type = "pip", variable = "(Intercept)",
##'      pip_level = attr(sel, "thresholds")[["(Intercept)"]])
##' }
##' @author Philippe Rast
##' @export
pip_fdr <- function(object, fdr = 0.05) {
  if (!inherits(object, "ivd")) {
    stop("`object` must be a fitted ivd model.", call. = FALSE)
  }
  if (!is.numeric(fdr) || length(fdr) != 1 || is.na(fdr) || fdr <= 0 || fdr >= 1) {
    stop("`fdr` must be a single number in (0, 1).", call. = FALSE)
  }

  df <- pip(object)
  parts <- lapply(split(df, df$scale_var), function(d) {
    sel <- .fdr_select(d$pip, fdr)
    d <- d[sel$order, ]
    d$local_fdr <- 1 - d$pip
    d$cum_fdr <- sel$cum_fdr
    d$selected <- sel$selected
    d
  })
  out <- do.call(rbind, parts)
  rownames(out) <- NULL

  thresholds <- vapply(parts, function(d) {
    if (any(d$selected)) min(d$pip[d$selected]) else NA_real_
  }, numeric(1))

  attr(out, "fdr") <- fdr
  attr(out, "thresholds") <- thresholds
  class(out) <- c("pip_fdr", "data.frame")
  out
}

##' Print method for `pip_fdr` objects
##' @title Print a pip_fdr selection
##' @param x An object of class `pip_fdr`.
##' @param ... Not used.
##' @return `x`, invisibly.
##' @author Philippe Rast
##' @export
print.pip_fdr <- function(x, ...) {
  fdr <- attr(x, "fdr")
  thresholds <- attr(x, "thresholds")
  cat(sprintf("Bayesian FDR decision rule (target FDR <= %s)\n\n", format(fdr)))

  for (v in names(thresholds)) {
    d <- x[x$scale_var == v, , drop = FALSE]
    n_sel <- sum(d$selected)
    cat(sprintf("Random scale effect: %s\n", v))
    if (n_sel == 0) {
      cat(sprintf("  No cluster can be selected at FDR <= %s (smallest local fdr: %.3f).\n\n",
                  format(fdr), min(d$local_fdr)))
      next
    }
    cat(sprintf("  Selected %d of %d clusters (implied PIP threshold >= %.3f, expected FDR = %.3f):\n",
                n_sel, nrow(d), thresholds[[v]], d$cum_fdr[n_sel]))
    sel <- d[d$selected, , drop = FALSE]
    print.data.frame(sel[, c("cluster_id", "pip", "local_fdr", "cum_fdr")],
                     row.names = FALSE, digits = 3)
    cat("\n")
  }
  invisible(x)
}
