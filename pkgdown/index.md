# Individual Variance Detection <img src="man/figures/logo.png" align="right" height="139" />

[![R-CMD-check](https://github.com/ph-rast/ivd/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ph-rast/ivd/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/consistentlyBetter/ivd/graph/badge.svg?token=SD0PM5BVIL)](https://codecov.io/gh/consistentlyBetter/ivd)

`ivd` implements the **Spike-and-Slab Mixed-Effects Location Scale Model** (SS-MELSM, [Carmo et al. 2026](https://doi.org/10.3102/10769986261426004))
to detect heterogeneous residual variances in hierarchical data.

While standard mixed-effects models assume constant within-group variability, mixed-effects location scale models
explicitly models the residual variance as a function of covariates and random effects. 
In the MELSM approach the covariance matrix (and thus the correlation matrix) of all random effects 
across location and scale is estimated jointly within the likelihood.
SS-MELSM uses a spike-and-slab prior to probabilistically identify specific units (e.g., schools, individuals) 
that exhibit unusually high or low consistency, distinguishing them from the population average.

## Features

- `ivd()` fits the SS-MELSM with `lme4`-style formulas for the location and
  scale sub-models; chains run in parallel, and any grouping IDs (character,
  factor, or non-consecutive numeric codes) are accepted.
- `summary()` and `plot()` report each unit's **posterior inclusion
  probability (PIP)**, with PIP, funnel, and outcome plots; results can be
  labelled with the original grouping IDs (`labels = "original"`).
- `pip()`, `fixef()`, `ranef()`, and `coef()` extract tidy estimates for
  further processing.
- `pip_sensitivity()` shows how the PIPs -- and the resulting
  classifications -- would change under different prior inclusion
  probabilities, computed analytically without refitting.
- `pip_fdr()` selects deviant units by controlling the Bayesian false
  discovery rate instead of an arbitrary fixed PIP cutoff, and
  `pip_diagnostics()` reports each PIP's Monte Carlo error and
  between-chain agreement.
- `pp_check()` runs a posterior predictive check of the within-cluster SDs,
  and `family = "student"` provides a robust t likelihood so that heavy
  tails are not mistaken for variance heterogeneity.
- Priors are tunable (`priors =`), chains are decoupled from workers
  (`chains =`), and `simulate_ivd()` generates data with known deviant
  clusters for power analysis and validation.
- `as.mcmc.list()` and `as_draws()` bridge fits into
  coda/bayesplot/tidybayes/posterior; `codaplot()` provides MCMC diagnostic
  plots with readable parameter names.

## Installation

You can install the development version of `ivd` from GitHub:

```r
# install.packages("devtools")
devtools::install_github("consistentlybetter/ivd")
```

## Acknowledgment

This work was supported by the Tools Competition catalyst award for the
project
[consistentlyBetter](https://tools-competition.org/winner/consistentlybetter/)
to PR. The content is solely the responsibility of the authors and does
not necessarily represent the official views of the funding agency.

## References


<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0" line-spacing="2">

<div id="ref-carmo2026" class="csl-entry">

Carmo, M., Williams, D. R., & Rast, P. (2026). Beyond average scores:
Identification of consistent and inconsistent academic achievement in
grouping units. *Journal of Educational and Behavioral Statistics*.
<https://doi.org/10.3102/10769986261426004>

</div>

<div id="ref-rodriguez2021" class="csl-entry">

Rodriguez, J. E., Williams, D. R., & Rast, P. (2024). Who is and is not
"average"? Random effects selection with spike-and-slab priors.
*Psychological Methods*. <https://doi.org/10.1037/met0000535>

</div>

<div id="ref-williams2022" class="csl-entry">

Williams, D. R., Martin, S. R., & Rast, P. (2022). Putting the
individual into reliability: Bayesian testing of homogeneous
within-person variance in hierarchical models. *Behavior Research
Methods*, *54*(3), 1272–1290.
<https://doi.org/10.3758/s13428-021-01646-x>

</div>

</div>
