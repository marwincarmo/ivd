##' Convert a fitted ivd model to a coda `mcmc.list`
##'
##' Returns the posterior samples as a [coda::mcmc.list()] with one element
##' per chain, restricted to the parameters shown by `summary()` and renamed
##' to the same human-readable labels (`Intc`, `scl_*`, `sd_*`, `R[...]`,
##' `pip[...]`). This makes a fit directly usable with the wider MCMC
##' ecosystem, e.g. `bayesplot::mcmc_trace()` or any coda function.
##' @title Convert an ivd fit to an mcmc.list
##' @param x An object of class `ivd`.
##' @param ... Not used.
##' @return A [coda::mcmc.list()] with renamed parameter columns.
##' @examples
##' \dontrun{
##' post <- as.mcmc.list(fit)
##' coda::effectiveSize(post)
##' bayesplot::mcmc_trace(post, pars = "Intc")
##' }
##' @author Philippe Rast
##' @exportS3Method coda::as.mcmc.list
as.mcmc.list.ivd <- function(x, ...) {
  coda::mcmc.list(lapply(.renamed_mcmc_list(x), coda::as.mcmc))
}

##' Convert a fitted ivd model to a posterior `draws` object
##'
##' Method for [posterior::as_draws()]: returns the posterior samples as a
##' `draws_array` (iterations x chains x variables) with the human-readable
##' parameter labels used by `summary()`. From there the whole
##' \pkg{posterior} / \pkg{tidybayes} / \pkg{bayesplot} toolchain applies
##' (e.g. `posterior::summarise_draws()`, `posterior::as_draws_df()`).
##' Requires the suggested \pkg{posterior} package.
##' @title Convert an ivd fit to a posterior draws object
##' @param x An object of class `ivd`.
##' @param ... Not used.
##' @return A [posterior::draws_array()] with renamed variables.
##' @examples
##' \dontrun{
##' posterior::summarise_draws(posterior::as_draws(fit))
##' }
##' @author Philippe Rast
##' @exportS3Method posterior::as_draws
as_draws.ivd <- function(x, ...) {
  .require_suggest("posterior", "`as_draws()`")
  chains <- .renamed_mcmc_list(x)
  vars <- colnames(chains[[1]])
  arr <- array(
    NA_real_,
    dim = c(nrow(chains[[1]]), length(chains), length(vars)),
    dimnames = list(iteration = NULL, chain = NULL, variable = vars)
  )
  for (ch in seq_along(chains)) arr[, ch, ] <- as.matrix(chains[[ch]])
  posterior::as_draws_array(arr)
}
