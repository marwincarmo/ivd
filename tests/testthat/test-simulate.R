## simulate_ivd(): pure R, fully covr-testable.

test_that("simulate_ivd returns coherent data and truth", {
  sim <- simulate_ivd(J = 30, n_j = 15, n_deviant = 3, seed = 1)

  expect_s3_class(sim, "ivd_sim")
  expect_named(sim, c("data", "truth", "params"))
  expect_equal(nrow(sim$data), 30 * 15)
  expect_equal(sort(unique(sim$data$id)), 1:30)
  expect_equal(as.integer(table(sim$data$id)), rep(15L, 30))
  expect_equal(nrow(sim$truth), 30)
  expect_equal(sum(sim$truth$deviant), 3)
  expect_setequal(sim$truth$id[sim$truth$deviant], sim$params$deviant_clusters)
  ## deviant clusters carry the factor, regular ones don't
  expect_true(all(sim$truth$factor[sim$truth$deviant] == 2))
  expect_true(all(sim$truth$factor[!sim$truth$deviant] == 1))
  expect_equal(sim$truth$sd_within, sim$truth$factor * 1)
})

test_that("simulate_ivd is reproducible and honours explicit clusters/factors", {
  s1 <- simulate_ivd(J = 10, n_j = 5, seed = 42)
  s2 <- simulate_ivd(J = 10, n_j = 5, seed = 42)
  expect_identical(s1, s2)

  sim <- simulate_ivd(J = 10, n_j = 5, deviant_clusters = c(7, 2),
                      deviant_factor = c(3, 0.5), seed = 1)
  ## indices are sorted; factors follow the sorted order
  expect_equal(sim$params$deviant_clusters, c(2L, 7L))
  expect_equal(sim$truth$factor[c(2, 7)], c(3, 0.5))
  expect_equal(sum(sim$truth$deviant), 2)

  ## n_deviant = 0 gives a null model
  s0 <- simulate_ivd(J = 5, n_j = 4, n_deviant = 0, seed = 1)
  expect_false(any(s0$truth$deviant))
})

test_that("simulated SDs match the requested generative values", {
  ## Large n_j: empirical within-cluster SDs concentrate on sd_within.
  sim <- simulate_ivd(J = 12, n_j = 4000, deviant_clusters = c(3, 9),
                      deviant_factor = 2.5, sigma_within = 1.5,
                      tau_loc = 0.8, beta = 10, seed = 7)
  emp_sd <- tapply(sim$data$y, sim$data$id, sd)
  expect_equal(unname(c(emp_sd)), sim$truth$sd_within, tolerance = 0.05)

  ## location side: cluster means ~ beta + u_loc
  emp_mean <- tapply(sim$data$y, sim$data$id, mean)
  expect_equal(unname(c(emp_mean)), 10 + sim$truth$u_loc, tolerance = 0.15)

  ## tau_loc = 0 collapses the location random effects
  s0 <- simulate_ivd(J = 6, n_j = 10, tau_loc = 0, seed = 1)
  expect_true(all(s0$truth$u_loc == 0))
})

test_that("simulate_ivd generates student-t residuals with the right scale", {
  ## t residuals: empirical SD -> scale * sqrt(df/(df-2))
  sim <- simulate_ivd(J = 8, n_j = 5000, n_deviant = 0, sigma_within = 2,
                      tau_loc = 0, family = "student", df = 5, seed = 11)
  emp_sd <- tapply(sim$data$y, sim$data$id, sd)
  expect_equal(unname(c(emp_sd)), rep(2 * sqrt(5 / 3), 8), tolerance = 0.1)
  expect_equal(sim$params$family, "student")
  expect_equal(sim$params$df, 5)

  ## gaussian sims don't carry a df
  expect_null(simulate_ivd(J = 4, n_j = 5, seed = 1)$params$df)

  ## df must exceed 2
  expect_error(simulate_ivd(J = 4, family = "student", df = 2), "> 2")

  ## print mentions the family
  expect_output(print(sim), "student-t residuals, df = 5")
  expect_output(print(sim), "within-cluster scale")
})

test_that("simulate_ivd supports unequal cluster sizes", {
  n_j <- c(5, 10, 15, 20)
  sim <- simulate_ivd(J = 4, n_j = n_j, n_deviant = 1, seed = 3)
  expect_equal(as.integer(table(sim$data$id)), n_j)
  expect_equal(sim$truth$n, n_j)
})

test_that("simulate_ivd validates its inputs", {
  expect_error(simulate_ivd(J = 1), "J")
  expect_error(simulate_ivd(J = 10, n_j = c(5, 5)), "length-J")
  expect_error(simulate_ivd(J = 10, n_deviant = 11), "0..J")
  expect_error(simulate_ivd(J = 10, deviant_clusters = c(2, 2)), "distinct")
  expect_error(simulate_ivd(J = 10, deviant_clusters = 12), "1..J")
  expect_error(simulate_ivd(J = 10, deviant_factor = -1), "positive")
  expect_error(simulate_ivd(J = 10, deviant_clusters = c(1, 2),
                            deviant_factor = c(2, 3, 4)), "one per")
  expect_error(simulate_ivd(J = 10, sigma_within = 0), "positive")
  expect_error(simulate_ivd(J = 10, tau_loc = -1), "non-negative")
})

test_that("print.ivd_sim summarises the design", {
  sim <- simulate_ivd(J = 8, n_j = 5, deviant_clusters = 4,
                      deviant_factor = 3, seed = 1)
  expect_output(print(sim), "8 clusters, 40 observations")
  expect_output(print(sim), "1 deviant cluster\\(s\\): 4 \\(SD factor 3\\)")
  expect_output(print(sim), "Fit with")

  s0 <- simulate_ivd(J = 5, n_j = 4, n_deviant = 0, seed = 1)
  expect_output(print(s0), "No deviant clusters")
})
