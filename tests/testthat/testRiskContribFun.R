context("testRiskContribFun.R")
library(riskParityPortfolio)
library(testthat)

test_that("risk contribution function works", {
  w <- c(1, 2, 3)
  Sigma <- diag(c(3, 2, 1))
  expect_that(all(g(w, Sigma) == c(3, 8, 9)), is_true())
})

test_that("risk contribution gradient function works", {
  w <- c(2, 3)
  Sigma <- matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE)
  answer <- matrix(c(10, 4, 9, 30), nrow = 2, byrow = TRUE)
  expect_that(all(g_grad(w, Sigma) == answer), is_true())
})
