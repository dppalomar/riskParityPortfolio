context("Additional checks")

set.seed(123)
N <- 10
V <- matrix(rnorm(N^2), N, N)
Sigma <- cov(V)

test_that("diag formulation throws error when other constraints or terms are
          included", {
  expect_error(riskParityPortfolio(Sigma, w_lb = 0.5, formulation = "diag"))
  expect_error(riskParityPortfolio(Sigma, lmb_var = 0.5, formulation = "diag"))
  expect_error(riskParityPortfolio(Sigma, mu = rep(1/N, N), formulation = "diag"))
})
