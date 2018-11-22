context("test-constraints.R")
library(riskParityPortfolio)
library(patrick)
library(testthat)

# generate a random Sigma and mu
set.seed(123)
N <- 5
V <- matrix(rnorm(N^2), N, N)
Sigma <- cov(V)

with_parameters_test_that("sca, alabama, and slsq give similar results in low dimensional problems", {
  formulations_list <- c("rc-double-index", "rc-over-b-double-index",
                         "rc-over-var vs b", "rc-over-var",
                         "rc-over-sd vs b-times-sd", "rc vs b-times-var",
                         "rc vs theta", "rc-over-b vs theta")
  for(formulation in formulations_list) {
    rpp_sca <- riskParityPortfolio(Sigma, budget = budget, shortselling = shortselling,
                                   method = "sca", formulation = formulation)
    rpp_alabama <- riskParityPortfolio(Sigma, budget = budget, shortselling = shortselling,
                                       method = "alabama", formulation = formulation)
    rpp_slsqp <- riskParityPortfolio(Sigma, budget = budget, shortselling = shortselling,
                                     method = "slsqp", formulation = formulation)
    expect_that(all(abs(rpp_sca$w - rpp_alabama$w) < 1e-3), is_true())
    expect_that(all(abs(rpp_sca$w - rpp_slsqp$w) < 1e-3), is_true())
  }
},
  cases(list(budget = TRUE, shortselling = FALSE),
        list(budget = TRUE, shortselling = TRUE),
        list(budget = FALSE, shortselling = FALSE)
  )
)
