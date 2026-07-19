## Monte Carlo diagnostics for the PIPs: per-chain PIPs, MCSE from the
## between-chain spread, and a chain-disagreement flag.

test_that("pip_diagnostics returns per-chain PIPs consistent with pip()", {
  skip_if(is.null(ivd_fixture), "fixture missing; run tests/testthat/fixtures/make-ivd-fixture.R")
  diag <- pip_diagnostics(ivd_fixture)

  expect_s3_class(diag, "pip_diagnostics")
  expect_equal(nrow(diag), nrow(pip(ivd_fixture)))
  chain_cols <- grep("^chain\\d+$", names(diag), value = TRUE)
  expect_length(chain_cols, ivd_fixture$workers)

  ## pooled pip = mean of the per-chain pips (equal-length chains), and both
  ## must agree with the pip() extractor
  expect_equal(diag$pip, rowMeans(diag[, chain_cols]))
  expect_equal(diag$pip, pip(ivd_fixture)$pip)

  ## per-chain pips are probabilities; range/mcse are coherent
  expect_true(all(diag[, chain_cols] >= 0 & diag[, chain_cols] <= 1))
  expect_equal(diag$range,
               apply(diag[, chain_cols], 1, max) - apply(diag[, chain_cols], 1, min))
  expect_true(all(diag$mcse >= 0))
  expect_true(all(diag$mcse <= diag$range / sqrt(length(chain_cols)) + 1e-12))
})

test_that("pip_diagnostics hand-check on a doctored two-chain fit", {
  skip_if(is.null(ivd_fixture), "fixture missing; run tests/testthat/fixtures/make-ivd-fixture.R")
  fit <- ivd_fixture
  Kr <- fit$nimble_constants$Kr

  ## force cluster 1 of the first scale effect to disagree: all-1 in chain 1,
  ## all-0 in chain 2 -> pip 0.5, mcse 0.25, range 1, unstable at any level
  col <- paste0("ss[", Kr + 1, ", 1]")
  fit$samples[[1]]$samples[, col] <- 1
  fit$samples[[2]]$samples[, col] <- 0

  diag <- pip_diagnostics(fit, pip_level = 0.75)
  row <- diag[diag$scale_var == colnames(fit$Z_scale)[1] & diag$cluster_index == 1, ]
  expect_equal(row$pip, 0.5)
  expect_equal(row$chain1, 1)
  expect_equal(row$chain2, 0)
  expect_equal(row$mcse, sd(c(1, 0)) / sqrt(2))
  expect_equal(row$range, 1)
  expect_true(row$unstable)

  ## agreement case: both chains all-1 -> no flag, zero mcse
  fit$samples[[2]]$samples[, col] <- 1
  agree <- pip_diagnostics(fit)
  agree_row <- agree[agree$scale_var == colnames(fit$Z_scale)[1] & agree$cluster_index == 1, ]
  expect_equal(agree_row$pip, 1)
  expect_equal(agree_row$mcse, 0)
  expect_false(agree_row$unstable)
})

test_that("unstable flag respects pip_level", {
  skip_if(is.null(ivd_fixture), "fixture missing; run tests/testthat/fixtures/make-ivd-fixture.R")
  fit <- ivd_fixture
  Kr <- fit$nimble_constants$Kr
  col <- paste0("ss[", Kr + 1, ", 2]")
  ## chain pips (0.6, 0.9): disagree at 0.75, agree at 0.5
  n1 <- nrow(fit$samples[[1]]$samples)
  n2 <- nrow(fit$samples[[2]]$samples)
  fit$samples[[1]]$samples[, col] <- c(rep(1, round(0.6 * n1)), rep(0, n1 - round(0.6 * n1)))
  fit$samples[[2]]$samples[, col] <- c(rep(1, round(0.9 * n2)), rep(0, n2 - round(0.9 * n2)))

  d75 <- pip_diagnostics(fit, pip_level = 0.75)
  d50 <- pip_diagnostics(fit, pip_level = 0.5)
  pick <- function(d) d[d$scale_var == colnames(fit$Z_scale)[1] & d$cluster_index == 2, ]
  expect_true(pick(d75)$unstable)
  expect_false(pick(d50)$unstable)
})

test_that("pip_diagnostics validates inputs and prints", {
  skip_if(is.null(ivd_fixture), "fixture missing; run tests/testthat/fixtures/make-ivd-fixture.R")
  expect_error(pip_diagnostics(list()), "fitted ivd model")
  expect_error(pip_diagnostics(ivd_fixture, pip_level = 2), "in \\(0, 1\\)")

  diag <- pip_diagnostics(ivd_fixture)
  expect_output(print(diag), "PIP Monte Carlo diagnostics")
  expect_output(print(diag), "Max. MCSE")
  if (any(diag$unstable)) {
    expect_output(print(diag), "classified inconsistently")
  } else {
    expect_output(print(diag), "classified consistently")
  }
})
