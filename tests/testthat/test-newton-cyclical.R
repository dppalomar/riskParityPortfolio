context("Newton cyclical")

N <- 5
V <- matrix(rnorm(N^2), N, N)
Sigma <- V %*% t(V)

test_that("Newton and Cyclical Coordinate Descent give the same results", {
  newton <- riskParityPortfolio(Sigma, method_init = "newton")
  cyclical_roncalli <- riskParityPortfolio(Sigma, method_init = "cyclical-roncalli")
  cyclical_spinu <- riskParityPortfolio(Sigma, method_init = "cyclical-spinu")
  expect_true(all(abs(newton$w - cyclical_roncalli$w) < 1e-4))
  expect_true(all(abs(newton$risk_contribution - cyclical_roncalli$risk_contribution) < 1e-4))
  expect_true(all(abs(newton$w - cyclical_spinu$w) < 1e-4), is_true())
  expect_true(all(abs(newton$risk_contribution - cyclical_spinu$risk_contribution) < 1e-4))
})
