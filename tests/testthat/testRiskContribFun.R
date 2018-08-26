context("testRiskContribFun.R")
library(riskParityPortfolio)
library(testthat)

test_that("risk contribution function works", {
  w <- c(1, 2, 3)
  Sigma <- diag(c(3, 2, 1))
  expect_that(all(g(w, Sigma) == c(3, 8, 9)), is_true())
})
