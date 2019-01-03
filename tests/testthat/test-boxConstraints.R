context("Box constraints")

# generate a random Sigma
set.seed(123)
N <- 10
V <- matrix(rnorm(N^2), N, N)
Sigma <- cov(V)

test_that("box constraints behave correctly", {
  formulations_list <- c("rc-double-index", "rc-over-b-double-index",
                         "rc-over-var vs b", "rc-over-var",
                         "rc-over-sd vs b-times-sd", "rc vs b-times-var",
                         "rc vs theta", "rc-over-b vs theta")
  for(formulation in formulations_list) {
    rpp <- riskParityPortfolio(Sigma, method = "sca", w_lb = -1, w_ub = 2, formulation = formulation)
    expect_that(all(rpp$w >= -1), is_true())
    expect_that(all(rpp$w <= 2), is_true())
  }
})

test_that("box constraints feasibility checks work", {
  expect_error(riskParityPortfolio(Sigma, w_lb = 0.5))
  expect_error(riskParityPortfolio(Sigma, w_lb = 2))
  expect_error(riskParityPortfolio(Sigma, w_ub = 0.01))
  expect_error(riskParityPortfolio(Sigma, w_ub = rep(0.01, N)))
  expect_silent(riskParityPortfolio(Sigma, w_lb = 0.05, w_ub = 0.3))
})

test_that("error controls work", {
  expect_warning(riskParityPortfolio(Sigma, w0 = rep(1, N)))
  expect_warning(riskParityPortfolio(Sigma, formulation = "diag", w0 = rep(1, N)))
  expect_warning(riskParityPortfolio(Sigma, formulation = "rc-double-index"))
  expect_warning(riskParityPortfolio(Sigma, formulation = "rc-double-index", mu = rep(0, N), w0 = rep(1, N)))
})
