context("testRiskGrad.R")
library(riskParityPortfolio)
library(testthat)

test_that("risk gradient formulations are equivalent", {
  w <- c(1, 2, 3)
  Sigma <- diag(c(3, 2, 1))
  N <- nrow(Sigma)
  x <- w * (Sigma %*% w)
  risk_grad_alt <- 4 * (Sigma * w + diag(Sigma %*% w)) %*% (N * x - sum(x) * rep(1, N))
  v <- (N * x - sum(x) * rep(1, N))
  risk_grad <- 4 * (Sigma %*% (w * v) + (Sigma %*% w) * v)
  expect_that(all(risk_grad_alt == risk_grad), is_true())
})
