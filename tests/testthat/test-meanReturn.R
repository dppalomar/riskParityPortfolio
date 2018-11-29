context("test-meanReturn.R")
library(riskParityPortfolio)
library(patrick)
library(testthat)

# generate a random Sigma and mu
N <- 5
V <- matrix(rnorm(N^2), N, N)
Sigma <- V %*% t(V)
mu <- runif(N)

with_parameters_test_that("sca and gensolver portfolios are consistent when
                          lambda equals zero or no mu is given", {
  rpp_lambda_zero <- risk_parity_portfolio_foo(Sigma, mu = mu, lmd_mu = 0)
  rpp_no_mu <- risk_parity_portfolio_foo(Sigma)
  expect_that(all.equal(rpp_lambda_zero$w, rpp_no_mu$w), is_true())
  expect_that(all.equal(rpp_lambda_zero$risk_contributions,
                        rpp_no_mu$risk_contributions), is_true())
},
  cases(list(risk_parity_portfolio_foo = riskParityPortfolioSCA),
        list(risk_parity_portfolio_foo = riskParityPortfolioGenSolver))
)
