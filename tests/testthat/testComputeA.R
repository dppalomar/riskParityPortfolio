context("testComputeA.R")
library(riskParityPortfolio)
library(testthat)

test_that("computations of A are equivalent", {
  N <- 2
  w <- runif(N)
  Sigma <- matrix(runif(N ^ 2), nrow = N)
  A1 <- compute_A_double_index_R(w, Sigma, N)
  A2 <- compute_A_double_index(w, Sigma, N)
  expect_that(all(A1 == A2), is_true())
  expect_that(ncol(A1) == ncol(A2), is_true())
  expect_that(nrow(A1) == nrow(A2), is_true())
})
