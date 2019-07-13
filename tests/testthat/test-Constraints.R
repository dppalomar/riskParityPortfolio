context("Constraints")

# generate a random Sigma
set.seed(123)
N <- 10
V <- matrix(rnorm(N^2), N, N)
Sigma <- cov(V)
formulations_list <- c("rc-double-index", "rc-over-b-double-index",
                       "rc-over-var vs b", "rc-over-var",
                       "rc-over-sd vs b-times-sd", "rc vs b-times-var",
                       "rc vs theta", "rc-over-b vs theta")
formulations_list_wo_theta <- c("rc-double-index", "rc-over-b-double-index",
                                "rc-over-var vs b", "rc-over-var",
                                "rc-over-sd vs b-times-sd", "rc vs b-times-var")

test_that("box constraints behave correctly", {
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

test_that("equality constraints behave correctly", {
  Cmat = matrix(0, N, N)
  Cmat[1, ] <- rep(1, N)
  Cmat[2, ] <- c(1, 0, 0, 1, rep(0, 6))
  for(formulation in formulations_list_wo_theta) {
    w_sum <- runif(1)
    w1_sum_w4 <- .2 * w_sum
    rpp <- riskParityPortfolio(Sigma, method = "sca",
                               Cmat = Cmat,
                               cvec = c(w_sum, w1_sum_w4, rep(0, 8)),
                               formulation = formulation)
    expect_that(abs(sum(rpp$w) - w_sum) < 1e-5, is_true())
    expect_that(abs(rpp$w[1] + rpp$w[4] - w1_sum_w4) < 1e-5, is_true())
  }
})

test_that("ineq and eq constraints behave correctly", {
  Cmat = matrix(0, N, N)
  Cmat[1, ] <- rep(1, N)
  Dmat = matrix(0, N, N)
  diag(Dmat) <- rep(-1, N)
  for(formulation in formulations_list_wo_theta) {
    w_sum <- runif(1)
    rpp <- riskParityPortfolio(Sigma, method = "sca",
                               Cmat = Cmat, Dmat = Dmat,
                               cvec = c(w_sum, rep(0, 9)),
                               dvec = c(rep(0, 10)), formulation = formulation)
    expect_that(all(rpp$w > 0), is_true())
    expect_that(abs(sum(rpp$w) - w_sum) < 1e-5, is_true())
  }
})
