## Ecosystem bridges: as.mcmc.list() and as_draws() expose the posterior
## samples with the renamed summary labels. All tests run against the
## committed fixture (2 chains, Kr = Sr = 2, J = 12).

test_that("as.mcmc.list.ivd returns a renamed mcmc.list", {
  skip_if(is.null(ivd_fixture), "fixture missing; run tests/testthat/fixtures/make-ivd-fixture.R")
  post <- coda::as.mcmc.list(ivd_fixture)

  expect_s3_class(post, "mcmc.list")
  expect_length(post, ivd_fixture$workers)
  nm <- colnames(post[[1]])
  ## human-readable labels, no raw NIMBLE names
  expect_true("Intc" %in% nm)
  expect_true(any(grepl("^scl_", nm)))
  expect_true(any(grepl("^sd_", nm)))
  expect_true(any(grepl("^pip\\[", nm)))
  expect_false(any(grepl("^beta\\[|^zeta\\[|^ss\\[", nm)))
  ## per-cluster u and any mu/tau are excluded (as in summary())
  expect_false(any(grepl("^u\\[|^mu\\[|^tau\\[", nm)))
  ## all chains share the same parameters
  expect_true(all(vapply(post, function(ch) identical(colnames(ch), nm), logical(1))))
})

test_that("as.mcmc.list.ivd values match the raw samples", {
  skip_if(is.null(ivd_fixture), "fixture missing; run tests/testthat/fixtures/make-ivd-fixture.R")
  post <- coda::as.mcmc.list(ivd_fixture)

  ## Intc is beta[1]; pooled means must agree with the raw draws
  raw_beta1 <- mean(vapply(ivd_fixture$samples,
                           function(ch) mean(ch$samples[, "beta[1]"]), numeric(1)))
  expect_equal(mean(vapply(post, function(ch) mean(ch[, "Intc"]), numeric(1))),
               raw_beta1, tolerance = 1e-12)
})

test_that("as_draws.ivd returns a draws_array usable by posterior", {
  skip_if(is.null(ivd_fixture), "fixture missing; run tests/testthat/fixtures/make-ivd-fixture.R")
  skip_if_not_installed("posterior")

  dr <- posterior::as_draws(ivd_fixture)
  expect_s3_class(dr, "draws_array")
  expect_equal(posterior::nchains(dr), ivd_fixture$workers)
  expect_true("Intc" %in% posterior::variables(dr))

  ## downstream tooling works and agrees with the mcmc.list bridge
  sm <- posterior::summarise_draws(posterior::subset_draws(dr, variable = "Intc"))
  post <- coda::as.mcmc.list(ivd_fixture)
  expect_equal(sm$mean,
               mean(vapply(post, function(ch) mean(ch[, "Intc"]), numeric(1))),
               tolerance = 1e-12)

  ## draws_df round-trip keeps all variables
  df <- posterior::as_draws_df(dr)
  expect_equal(posterior::nvariables(df), posterior::nvariables(dr))
})

test_that("bridges agree with codaplot's renaming (shared helper)", {
  skip_if(is.null(ivd_fixture), "fixture missing; run tests/testthat/fixtures/make-ivd-fixture.R")
  helper <- ivd:::.renamed_mcmc_list(ivd_fixture)
  post <- coda::as.mcmc.list(ivd_fixture)
  expect_identical(colnames(post[[1]]), colnames(helper[[1]]))
})
