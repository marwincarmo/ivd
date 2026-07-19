## Bayesian FDR decision rule. The core is unit-tested on hand-computed
## values; pip_fdr() itself runs against the committed fixture.

test_that(".fdr_select matches a hand computation", {
  ## local fdr: 0, .1, .2, .5 -> cumulative FDR: 0, .05, .1, .2
  pips <- c(1, 0.9, 0.8, 0.5)
  sel <- ivd:::.fdr_select(pips, fdr = 0.10)

  expect_equal(sel$cum_fdr, c(0, 0.05, 0.1, 0.2))
  expect_equal(sel$selected, c(TRUE, TRUE, TRUE, FALSE))

  ## stricter target drops one more; looser keeps all
  expect_equal(ivd:::.fdr_select(pips, 0.06)$selected, c(TRUE, TRUE, FALSE, FALSE))
  expect_equal(ivd:::.fdr_select(pips, 0.25)$selected, rep(TRUE, 4))
})

test_that(".fdr_select handles unordered input and empty selection", {
  pips <- c(0.5, 1, 0.8, 0.9) # unsorted
  sel <- ivd:::.fdr_select(pips, fdr = 0.10)
  expect_equal(pips[sel$order], c(1, 0.9, 0.8, 0.5))
  expect_equal(sum(sel$selected), 3)

  ## nothing selectable when even the best cluster exceeds the target
  none <- ivd:::.fdr_select(c(0.6, 0.5), fdr = 0.05)
  expect_false(any(none$selected))
})

test_that("pip_fdr returns a coherent selection on the fixture", {
  skip_if(is.null(ivd_fixture), "fixture missing; run tests/testthat/fixtures/make-ivd-fixture.R")
  res <- pip_fdr(ivd_fixture, fdr = 0.25)

  expect_s3_class(res, "pip_fdr")
  expect_setequal(unique(res$scale_var), unique(pip(ivd_fixture)$scale_var))
  expect_equal(nrow(res), nrow(pip(ivd_fixture)))
  expect_equal(attr(res, "fdr"), 0.25)

  for (v in unique(res$scale_var)) {
    d <- res[res$scale_var == v, ]
    ## ordered by decreasing pip; selected clusters are a prefix
    expect_true(all(diff(d$pip) <= 0))
    expect_true(all(diff(d$selected) <= 0))
    if (any(d$selected)) {
      k <- sum(d$selected)
      ## achieved FDR within target; threshold attribute consistent
      expect_lte(d$cum_fdr[k], 0.25)
      expect_equal(attr(res, "thresholds")[[v]], min(d$pip[d$selected]))
      ## every selected pip >= every unselected pip
      if (k < nrow(d)) expect_gte(min(d$pip[d$selected]), max(d$pip[!d$selected]))
    }
  }
})

test_that("pip_fdr threshold adapts to the target level", {
  skip_if(is.null(ivd_fixture), "fixture missing; run tests/testthat/fixtures/make-ivd-fixture.R")
  strict <- pip_fdr(ivd_fixture, fdr = 0.01)
  loose <- pip_fdr(ivd_fixture, fdr = 0.5)
  expect_lte(sum(strict$selected), sum(loose$selected))
})

test_that("pip_fdr validates its inputs", {
  skip_if(is.null(ivd_fixture), "fixture missing; run tests/testthat/fixtures/make-ivd-fixture.R")
  expect_error(pip_fdr(list()), "fitted ivd model")
  expect_error(pip_fdr(ivd_fixture, fdr = 0), "in \\(0, 1\\)")
  expect_error(pip_fdr(ivd_fixture, fdr = 1.2), "in \\(0, 1\\)")
})

test_that("print.pip_fdr renders both the selected and the empty case", {
  skip_if(is.null(ivd_fixture), "fixture missing; run tests/testthat/fixtures/make-ivd-fixture.R")
  loose <- pip_fdr(ivd_fixture, fdr = 0.6)
  expect_output(print(loose), "Bayesian FDR decision rule")
  expect_output(print(loose), "Selected \\d+ of \\d+ clusters")

  strict <- pip_fdr(ivd_fixture, fdr = 0.001)
  if (!any(strict$selected)) {
    expect_output(print(strict), "No cluster can be selected")
  }
})
