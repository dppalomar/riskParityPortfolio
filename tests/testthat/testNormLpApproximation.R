context("testNormLpApproximation.R")
library(riskParityPortfolio)
library(testthat)

test_that("lp norm approximation", {
  x <- c(-runif(50), runif(50))
  r <- rho(x, p = 1, e = 1e-8)
  expect_that(all((abs(abs(x) - r)) < 1e-1), is_true())
})

test_that("gradient of rho", {
  x <- c(-runif(50), runif(50))
  rg <- rho_grad(x, p = 1, e = 1e-8)
  expect_that(all(sign(rg) == sign(x)), is_true())
})
