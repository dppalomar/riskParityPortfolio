context("test-newton-cyclical.R")
library(riskParityPortfolio)
library(testthat)

N <- 5
V <- matrix(rnorm(N^2), N, N)
Sigma <- V %*% t(V)

test_that("Newton and Cyclical Coordinate Descent give the same results", {
  newton <- riskParityPortfolioNewton(Sigma)
  cyclical <- riskParityPortfolioCyclical(Sigma)
  expect_that(all(abs(newton$w - cyclical$w) < 1e-4), is_true())
  expect_that(all(abs(newton$risk_contribution - cyclical$risk_contribution) < 1e-4),
              is_true())
})
