# ivd (development version)

## New features

* New `chains` argument in `ivd()`, decoupled from `workers`: chains are
  distributed over the worker processes, and each worker compiles the model
  once and reuses it for its chains -- so extra chains cost sampling time
  but no additional compilation. `chains` defaults to `workers`, preserving
  the behaviour (and exact draws) of existing code. `summary()` now labels
  the count "Chains:" instead of "Chains (workers):".
* New `family = "student"` option in `ivd()`: a student-t likelihood with
  estimated degrees of freedom (`nu = 2 + gamma(shape, rate)`, reported as
  `nu` in `summary()`; prior tunable via `priors = list(nu = ...)`). Heavy
  tails can masquerade as variance heterogeneity under the gaussian
  likelihood, inflating PIPs of clusters that merely contain outliers --
  `"student"` absorbs the tails instead. `pp_check()`, WAIC/logLik, plots
  and `simulate_ivd()` (via its own `family`/`df` arguments) all support it.
* New `simulate_ivd()`: simulates data from the intercept-only MELSM with a
  known subset of clusters whose within-cluster SD is inflated (or deflated)
  by a chosen factor. Returns the data ready for `ivd()` plus the ground
  truth per cluster -- for power analysis, teaching, and recovery checks.
* New `pp_check()` method: posterior predictive check of the within-cluster
  SDs -- the quantity the spike-and-slab makes claims about. Default view
  compares each cluster's observed SD with its predictive interval
  (outliers highlighted and labelled); `type = "density"` overlays the
  observed cluster-SD density on replicated ones. Systematic misfit
  suggests the normal likelihood is inadequate (e.g. heavy tails).
* New `priors` argument in `ivd()`: the prior hyperparameters -- location
  intercept and coefficients (normal), scale coefficients (normal),
  random-effect SDs (half-t df/scale), and the LKJ shape -- can now be set
  by the user, e.g. `priors = list(zeta = c(sd = 1), lkj_eta = 2)`. Partial
  specifications are filled with the (unchanged) defaults, and the resolved
  specification is stored on the fit as `$priors`.
* New `pip_diagnostics()`: Monte Carlo diagnostics for the PIPs, which
  binary-indicator Rhat does not cover. Reports per-chain PIPs, the Monte
  Carlo standard error of each pooled PIP, and flags clusters whose
  classification at `pip_level` differs between chains (i.e. PIPs that
  should not be trusted without more iterations).
* New `pip_fdr()`: a Bayesian false-discovery-rate decision rule for
  selecting deviant clusters (Newton et al., 2004). Instead of a fixed PIP
  cutoff, it selects the largest cluster set whose expected FDR stays below
  a target level and reports the implied, data-adaptive PIP threshold
  (usable as `pip_level` in `plot()`).
* Ecosystem bridges: `as.mcmc.list()` (coda) and `as_draws()` (posterior)
  convert a fit to standard formats with the human-readable parameter labels
  used by `summary()`, opening up bayesplot/tidybayes/posterior tooling
  (e.g. `bayesplot::mcmc_trace(as.mcmc.list(fit), pars = "Intc")`).

* `ivd()` accepts arbitrary grouping IDs (character, factor, or numeric with
  gaps, e.g. real school codes). IDs are recoded internally to the 1..J index
  NIMBLE needs -- without reordering the data -- and the original labels are
  stored on the fitted object as `group_labels`.
* `summary()` and `plot()` gained a shared `labels = c("index", "original")`
  argument: the default keeps the compact internal cluster index in both
  displays; `"original"` switches both to the user's own grouping IDs.
* New `print()` method for fitted objects: formulas, data dimensions,
  sampling setup, and a one-line convergence note (max split-Rhat, min
  n_eff).
* New `pip()` extractor returning the posterior inclusion probabilities as a
  tidy data frame (one row per cluster x scale random effect, with the
  original cluster IDs).
* New `fixef()`, `ranef()`, and `coef()` extractors following the lme4/nlme
  conventions (without adding a dependency on either).
* New `pip_sensitivity()`: analytic prior-sensitivity for the PIPs. The
  posterior odds are converted to a prior-independent Bayes factor and
  re-weighted across a grid of prior inclusion probabilities -- no refit
  needed. Comes with a `plot()` method showing per-cluster PIP trajectories
  (prior-sensitive clusters highlighted; subset with `clusters =`).
* `summary()` now returns a structured `summary.ivd` object (`$table`,
  `$waic`, `$lppd`, `$pwaic`, `$chains`) rendered by a separate `print()`
  method; console output is unchanged.

## Bug fixes

* `codaplot()` no longer requires the user to attach coda: plot types such as
  `"traceplot"` are now resolved in the coda namespace, so
  `ivd::codaplot(fit)` works from a plain script (previously it failed with
  `object 'traceplot' not found` unless `library(coda)` had been called).
* Fixed a crash in the default `n_eff = "local"` diagnostics on short or
  strongly autocorrelated chains (Geyer truncation over an empty set
  produced `1:Inf`).
* `.autocorrelation_fft()` now actually zero-pads, computing a linear rather
  than circular autocorrelation; results match `stats::acf()`.
* Fixed the row offset linking PIPs to clusters in `summary()` and
  `codaplot()` for models where the number of location and scale random
  effects differ.
* `DESCRIPTION` now correctly states that the spike-and-slab prior is on the
  *scale* random effects, and cites the published paper:
  Carmo, Williams & Rast (2026) <doi:10.3102/10769986261426004>.

## ivd 1.0.0

Initial CRAN release
