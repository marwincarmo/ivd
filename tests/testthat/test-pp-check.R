## Posterior predictive check of the within-cluster SDs. The statistics core
## is pure R and runs against the fixture plus a deterministic fake fit.

## Minimal hand-buildable fit: 2 clusters x 3 obs, Kr = Sr = 1, one chain of
## constant draws so mu/tau are exactly determined by the set values.
.fake_fit <- function(beta = 1, zeta = log(0.5), u = c(0, 0, 0, 0), n_iter = 4) {
  ## u = (u_loc[1], u_loc[2], u_scl[1], u_scl[2]); P = 2
  group_id <- c(1, 1, 1, 2, 2, 2)
  X <- matrix(1, 6, 1)
  draw <- c("beta[1]" = beta, "zeta[1]" = zeta,
            "u[1, 1]" = u[1], "u[2, 1]" = u[2],
            "u[1, 2]" = u[3], "u[2, 2]" = u[4])
  mat <- matrix(rep(draw, each = n_iter), nrow = n_iter,
                dimnames = list(NULL, names(draw)))
  structure(list(
    samples = list(list(samples = mat)),
    X = X, Z = X, X_scale = X, Z_scale = X,
    Y = data.frame(group_id = group_id, Y = c(0.9, 1.1, 1.0, 3.4, 2.6, 3.0)),
    nimble_constants = list(Kr = 1, Sr = 1, J = 2),
    group_labels = c("g1", "g2"),
    workers = 1
  ), class = c("ivd", "list"))
}

test_that(".pp_cluster_sds reproduces the model's within-cluster SD", {
  ## tau = exp(zeta + u_scl) is tiny -> replicated Y == mu almost exactly,
  ## so the replicated within-cluster SD equals the SD of mu (= 0 here since
  ## mu is constant within cluster).
  fit <- .fake_fit(zeta = log(1e-8))
  stats <- ivd:::.pp_cluster_sds(fit, ndraws = 4, seed = 1)

  expect_equal(dim(stats$rep), c(4, 2))
  expect_true(all(stats$rep < 1e-6))
  ## observed SDs are those of the raw data
  expect_equal(unname(stats$observed),
               c(sd(c(0.9, 1.1, 1.0)), sd(c(3.4, 2.6, 3.0))))

  ## with tau = 2 for cluster 2 only (u_scl[2] = log(4): tau = 0.5*4 = 2),
  ## the mean replicated SD approaches the sampling mean of an SD of n=3
  ## normals: E[S] = sigma * sqrt(2/(n-1)) * gamma(n/2) / gamma((n-1)/2)
  fit2 <- .fake_fit(zeta = log(0.5), u = c(0, 0, 0, log(4)), n_iter = 400)
  stats2 <- ivd:::.pp_cluster_sds(fit2, ndraws = 400, seed = 2)
  e_s <- function(sigma, n) sigma * sqrt(2 / (n - 1)) * gamma(n / 2) / gamma((n - 1) / 2)
  expect_equal(mean(stats2$rep[, 1]), e_s(0.5, 3), tolerance = 0.1)
  expect_equal(mean(stats2$rep[, 2]), e_s(2, 3), tolerance = 0.1)
})

test_that(".pp_cluster_sds is reproducible with a seed and caps ndraws", {
  skip_if(is.null(ivd_fixture), "fixture missing; run tests/testthat/fixtures/make-ivd-fixture.R")
  s1 <- ivd:::.pp_cluster_sds(ivd_fixture, ndraws = 20, seed = 99)
  s2 <- ivd:::.pp_cluster_sds(ivd_fixture, ndraws = 20, seed = 99)
  expect_identical(s1, s2)

  total <- sum(vapply(ivd_fixture$samples, function(ch) nrow(ch$samples), integer(1)))
  s_all <- ivd:::.pp_cluster_sds(ivd_fixture, ndraws = total + 500, seed = 1)
  expect_equal(nrow(s_all$rep), total)
})

test_that(".pp_cluster_sds drops single-observation clusters", {
  fit <- .fake_fit()
  ## make cluster 2 a singleton
  fit$Y <- fit$Y[c(1, 2, 3, 4), ]
  fit$X <- fit$Z <- fit$X_scale <- fit$Z_scale <- matrix(1, 4, 1)
  stats <- ivd:::.pp_cluster_sds(fit, ndraws = 2, seed = 1)
  expect_equal(stats$dropped, 2L)
  expect_equal(ncol(stats$rep), 1)
  expect_length(stats$observed, 1)
})

test_that("pp_check.ivd builds both plot types on the fixture", {
  skip_if(is.null(ivd_fixture), "fixture missing; run tests/testthat/fixtures/make-ivd-fixture.R")
  p_int <- pp_check(ivd_fixture, ndraws = 30, seed = 7)
  expect_s3_class(p_int, "ggplot")
  expect_true(all(c("obs", "med", "lo_outer", "hi_outer", "outside") %in% names(p_int$data)))
  expect_equal(nrow(p_int$data), ivd_fixture$nimble_constants$J)
  ## intervals are ordered coherently
  expect_true(all(p_int$data$lo_outer <= p_int$data$lo_inner))
  expect_true(all(p_int$data$hi_inner <= p_int$data$hi_outer))

  p_dens <- pp_check(ivd_fixture, ndraws = 10, type = "density", seed = 7)
  expect_s3_class(p_dens, "ggplot")
})

test_that("pp_check.ivd honours labels = 'original' and validates probs", {
  skip_if(is.null(ivd_fixture), "fixture missing; run tests/testthat/fixtures/make-ivd-fixture.R")
  fit <- ivd_fixture
  J <- fit$nimble_constants$J
  fit$group_labels <- sprintf("school_%03d", seq_len(J))
  p <- pp_check(fit, ndraws = 10, labels = "original", seed = 1)
  expect_true(all(grepl("^school_", p$data$label)))

  fit$group_labels <- NULL
  expect_warning(pp_check(fit, ndraws = 5, labels = "original"),
                 "predates 'group_labels'")

  expect_error(pp_check(ivd_fixture, probs = c(0.9, 0.5)), "increasing")
  expect_error(pp_check(ivd_fixture, probs = 0.9), "increasing")
})
