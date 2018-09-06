context("testComputeA.R")
library(riskParityPortfolio)
library(testthat)

test_that("computations of A are equivalent", {
  N <- 3
  w <- runif(N)
  Sigma <- matrix(runif(N ^ 2), nrow = N)
  A1 <- compute_A(w, Sigma)
  A2 <- computeACpp(w, Sigma)
  expect_that(all(A1 == A2), is_true())
  expect_that(ncol(A1) == ncol(A2), is_true())
  expect_that(nrow(A1) == nrow(A2), is_true())
})
