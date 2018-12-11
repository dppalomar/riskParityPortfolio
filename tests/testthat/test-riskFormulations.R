context("test-riskFormulations.R")
library(riskParityPortfolio)
library(testthat)

# generate a random Sigma and mu
set.seed(123)
N <- 5
V <- matrix(rnorm(N^2), N, N)
Sigma <- cov(V)

test_that("sca, alabama, and slsq give similar results in low dimensional problems", {
  formulations_list <- c("rc-double-index", "rc-over-b-double-index",
                         "rc-over-var vs b", "rc-over-var",
                         "rc-over-sd vs b-times-sd", "rc vs b-times-var",
                         "rc vs theta", "rc-over-b vs theta")
  for(formulation in formulations_list) {
    rpp_sca <- riskParityPortfolio(Sigma, method = "sca", formulation = formulation)
    rpp_alabama <- riskParityPortfolio(Sigma, method = "alabama", formulation = formulation)
    rpp_slsqp <- riskParityPortfolio(Sigma, method = "slsqp", formulation = formulation)
    expect_that(all(abs(rpp_sca$w - rpp_alabama$w) < 1e-3), is_true())
    expect_that(all(abs(rpp_sca$w - rpp_slsqp$w) < 1e-3), is_true())
  }
})
