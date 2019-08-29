context("Mean return")

# generate a random Sigma and mu
N <- 5
V <- matrix(rnorm(100 * N), nrow=N)
Sigma <- V %*% t(V)
mu <- runif(N)

test_that("sca and gensolver portfolios are consistent when lmd_mu equals zero
          or no mu is given", {
  expect_warning(rpp_lmb_zero <- riskParityPortfolio(Sigma, mu = mu, lmd_mu = 0, method = "sca"))
  rpp_no_mu <- riskParityPortfolio(Sigma, method = "sca")
  expect_equal(rpp_lmb_zero$w, rpp_no_mu$w)
  expect_equal(rpp_lmb_zero$relative_risk_contribution, rpp_no_mu$relative_risk_contribution)

  rpp_no_mu <- riskParityPortfolio(Sigma, method = "alabama")
  expect_warning(rpp_lmb_zero <- riskParityPortfolio(Sigma, mu = mu, lmd_mu = 0, method = "alabama"))
  expect_equal(rpp_lmb_zero$w, rpp_no_mu$w)
  expect_equal(rpp_lmb_zero$relative_risk_contribution, rpp_no_mu$relative_risk_contribution)
})

test_that("roncalli's formulation with mu converges to roncalli's formulation
          without mu for sufficiently large mean_volatility_tradeoff", {
  b <- rep(1/N, N)
  mean_volatility_tradeoff <- 1e6
  risk_free_return <- 0
  rpp_mu <- riskParityPortfolio:::active_risk_parity_portfolio_ccd(Sigma, b, mu,
                                                                   mean_volatility_tradeoff,
                                                                   risk_free_return, 1e-6, 50)
  rpp <- riskParityPortfolio:::risk_parity_portfolio_ccd_roncalli(Sigma, b, 1e-6, 50)
  rc <- rpp_mu * (Sigma %*% rpp_mu)
  expect_true(all.equal(rpp_mu, rpp, tolerance = 1e-4))
  expect_true(all.equal(c(round(rc / sum(rc), digits=2)), b, tolerance = 1e-1))
})
