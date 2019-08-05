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
    expect_true(all(rpp$w >= -1))
    expect_true(all(rpp$w <= 2))
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
  Cmat <- matrix(1, 1, N)
  for(formulation in formulations_list_wo_theta) {
    w_sum <- runif(1)
    rpp <- riskParityPortfolio(Sigma, method = "sca",
                               Cmat = Cmat,
                               cvec = c(w_sum),
                               formulation = formulation)
    expect_true(abs(sum(rpp$w) - w_sum) < 1e-5)
  }
})

test_that("ineq and eq constraints behave correctly", {
  N <- 10
  V <- matrix(rnorm(1000 * N), nrow=N)
  Sigma <- cov(t(V))
  Cmat <- matrix(1, 1, N)
  Dmat <- matrix(0, N, N)
  diag(Dmat) <- rep(-1, N)
  for(formulation in formulations_list_wo_theta) {
    rpp <- riskParityPortfolio(Sigma, method = "sca",
                               Cmat = Cmat, Dmat = Dmat,
                               cvec = c(1),
                               dvec = c(rep(0, N)), formulation = formulation)
    expect_true(all(rpp$w > 0))
    expect_true(abs(sum(rpp$w) - 1) < 1e-5)
  }
})

test_that("rpp_with_equality_constraints_iteration agree with solve.QP", {
  N <- sample(c(1:10), 1)
  Qk <- diag(c(1:N))
  qk <- rep(1, N)
  Cmat <- matrix(1, 1, N)
  cvec <- c(runif(1))

  c_solution <- riskParityPortfolio:::rpp_equality_constraints_iteration(Cmat, cvec, Qk, qk)
  qp_solution <- quadprog::solve.QP(Qk, -qk, Amat = t(Cmat), bvec = cvec, meq = 1)$solution
  expect_true(all(abs(c_solution - qp_solution) < 1e-5))
})
