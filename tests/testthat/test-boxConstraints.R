context("test-boxConstraints.R")
library(riskParityPortfolio)
library(testthat)

# generate a random Sigma and mu
set.seed(123)
N <- 10
V <- matrix(rnorm(N^2), N, N)
Sigma <- cov(V)
w_lb <- .01
w_ub <- .2

test_that("box constraints behave correctly", {
  formulations_list <- c("rc-double-index", "rc-over-b-double-index",
                         "rc-over-var vs b", "rc-over-var",
                         "rc-over-sd vs b-times-sd", "rc vs b-times-var",
                         "rc vs theta", "rc-over-b vs theta")
  for(formulation in formulations_list) {
    rpp <- riskParityPortfolio(Sigma, method = "sca", w_lb = w_lb, w_ub = w_ub,
                               formulation = formulation)
    expect_that(all(rpp$w >= w_lb), is_true())
    expect_that(all(rpp$w <= w_ub), is_true())
  }
})
