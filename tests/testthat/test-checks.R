context("Additional checks")

set.seed(123)
N <- 10
V <- matrix(rnorm(N^2), N, N)
Sigma <- cov(V)
colnames(Sigma) <- as.character(c(1:10))


test_that("diag formulation throws error when other constraints or terms are
          included", {
  expect_error(riskParityPortfolio(Sigma, w_lb = 0.05, formulation = "diag"))
  expect_error(riskParityPortfolio(Sigma, lmd_var = 0.5, formulation = "diag"))
  expect_error(riskParityPortfolio(Sigma, mu = rep(1/N, N), formulation = "diag"))
})


test_that("an error is raised when method is not supported", {
  expect_error(riskParityPortfolio(Sigma, lmd_var = .1, formulation = "rc-over-var", method = "NA"))
})


test_that("a warning is raised when initial point is not feasible", {
  expect_warning(riskParityPortfolio(Sigma, formulation = "rc-over-var", w0 = rep(1, N)))
})


test_that("stocks names remain consistent with input parameters", {
  portfolio <- riskParityPortfolio(Sigma)
  expect_equal(names(portfolio$w), colnames(Sigma))
  expect_equal(names(portfolio$relative_risk_contribution), colnames(Sigma))
})

