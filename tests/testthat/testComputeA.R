context("testComputeA.R")
library(riskParityPortfolio)
library(testthat)

test_that("computations of A are equivalent", {
  w <- c(1:2)
  Sigma <- matrix(c(1:4), nrow = 2)
  A1 <- compute_A(w, Sigma)
  A2 <- computeA(w, Sigma)
  expect_that(all(A1 == A2), is_true())
  expect_that(ncol(A1) == ncol(A2), is_true())
  expect_that(nrow(A1) == nrow(A2), is_true())
})
