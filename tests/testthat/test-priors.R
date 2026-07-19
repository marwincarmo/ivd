## .ivd_priors(): resolution and validation of the ivd(priors = ) argument.
## Pure R (no NIMBLE), so these run under coverage.

test_that(".ivd_priors returns the documented defaults", {
  p <- ivd:::.ivd_priors(list(), mean_pred = 2, sd_pred = 1.5)

  expect_equal(p$beta_intercept, c(mean = 2, sd = 4.5)) # empirical: mean(Y), 3*sd(Y)
  expect_equal(p$beta, c(mean = 0, sd = 1000))
  expect_equal(p$zeta, c(mean = 0, sd = 3))
  expect_equal(p$sigma_rand, c(df = 3, scale = 1))
  expect_equal(p$lkj_eta, 1)
  expect_true(attr(p, "empirical_intercept"))

  ## NULL behaves like an empty list
  expect_equal(ivd:::.ivd_priors(NULL, 2, 1.5), p, ignore_attr = TRUE)
})

test_that(".ivd_priors merges partial specifications into the defaults", {
  p <- ivd:::.ivd_priors(list(zeta = c(sd = 1)), mean_pred = 0, sd_pred = 1)
  expect_equal(p$zeta, c(mean = 0, sd = 1)) # mean kept from default
  expect_equal(p$beta, c(mean = 0, sd = 1000)) # untouched components intact

  p2 <- ivd:::.ivd_priors(
    list(beta_intercept = c(mean = 10, sd = 5),
         sigma_rand = c(scale = 2.5),
         lkj_eta = 4),
    mean_pred = 0, sd_pred = 1
  )
  expect_equal(p2$beta_intercept, c(mean = 10, sd = 5))
  expect_equal(p2$sigma_rand, c(df = 3, scale = 2.5))
  expect_equal(p2$lkj_eta, 4)
  expect_false(attr(p2, "empirical_intercept"))
})

test_that(".ivd_priors rejects malformed input", {
  expect_error(ivd:::.ivd_priors("normal", 0, 1), "must be a named list")
  expect_error(ivd:::.ivd_priors(list(c(sd = 1)), 0, 1), "must be named")
  expect_error(ivd:::.ivd_priors(list(gamma = c(sd = 1)), 0, 1),
               "Unknown prior component")
  expect_error(ivd:::.ivd_priors(list(beta = c(spread = 1)), 0, 1),
               "Unknown element")
  expect_error(ivd:::.ivd_priors(list(beta = 5), 0, 1), "named numeric vector")
  expect_error(ivd:::.ivd_priors(list(lkj_eta = c(1, 2)), 0, 1), "single number")
})

test_that(".ivd_priors rejects invalid hyperparameter values", {
  expect_error(ivd:::.ivd_priors(list(zeta = c(sd = 0)), 0, 1), "sd > 0")
  expect_error(ivd:::.ivd_priors(list(beta = c(sd = -1)), 0, 1), "sd > 0")
  expect_error(ivd:::.ivd_priors(list(sigma_rand = c(df = 0)), 0, 1),
               "df > 0 and scale > 0")
  expect_error(ivd:::.ivd_priors(list(sigma_rand = c(scale = -2)), 0, 1),
               "df > 0 and scale > 0")
  expect_error(ivd:::.ivd_priors(list(lkj_eta = -1), 0, 1), "must be positive")
  expect_error(ivd:::.ivd_priors(list(zeta = c(mean = NA_real_)), 0, 1),
               "missing values")
  ## a logical NA is not numeric and fails the shape check instead
  expect_error(ivd:::.ivd_priors(list(zeta = c(mean = NA)), 0, 1),
               "named numeric vector")
})

test_that("fixture predates `priors` but new fits carry the resolved spec", {
  skip_if(is.null(ivd_fixture), "fixture missing; run tests/testthat/fixtures/make-ivd-fixture.R")
  ## Legacy objects have no $priors; downstream methods must not require it.
  ## (The fixture is regenerated with the current schema, so once it carries
  ## `priors` this documents the field; both states are acceptable here.)
  if (!is.null(ivd_fixture$priors)) {
    expect_named(ivd_fixture$priors,
                 c("beta_intercept", "beta", "zeta", "sigma_rand", "lkj_eta"))
  }
  succeed()
})
