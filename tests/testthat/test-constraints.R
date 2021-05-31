context("Constraints")

# generate a random Sigma
set.seed(123)
N <- 10
T <- 10000
V <- matrix(rnorm(T), T/N, N)
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
})

test_that("equality constraints behave correctly", {
  Cmat <- matrix(0, 1, nrow(Sigma))
  Cmat[1] <- 1
  Cmat[2] <- 1
  for(formulation in formulations_list_wo_theta) {
    w12_sum <- runif(1)
    rpp <- riskParityPortfolio(Sigma, method = "sca", Cmat = Cmat, cvec = c(w12_sum),
                               formulation = formulation)
    expect_true(abs(sum(rpp$w) - 1) < 1e-4)
    expect_true(abs(Cmat %*% rpp$w - w12_sum) < 1e-4)
  }
})

test_that("ineq and eq constraints behave correctly", {
  Dmat <- matrix(0, 1, nrow(Sigma))
  Dmat[1] <- 1
  Dmat[2] <- 1
  for(formulation in formulations_list_wo_theta) {
    w12_sum <- runif(1)
    rpp <- riskParityPortfolio(Sigma, method = "sca", Dmat = Dmat, dvec = c(w12_sum),
                               formulation = formulation)
    expect_true(Dmat %*% rpp$w <= (w12_sum + 1e-5))
    expect_true(all(rpp$w > 0))
    expect_true(abs(sum(rpp$w) - 1) < 1e-5)
  }
})

test_that("rpp_with_equality_constraints_iteration agree with solve.QP", {
  N <- sample(c(2:10), 1)
  Qk <- diag(c(1:N))
  qk <- rep(1, N)
  Cmat <- matrix(1, 1, N)
  cvec <- c(runif(1))

  c_solution <- riskParityPortfolio:::rpp_equality_constraints_iteration(Cmat, cvec, Qk, qk)
  qp_solution <- quadprog::solve.QP(Qk, -qk, Amat = t(Cmat), bvec = cvec, meq = 1)$solution
  r_solution <- riskParityPortfolio:::rpp_equality_constraints_iteration_R(Cmat, cvec, Qk, qk)
  expect_true(all(abs(c_solution - qp_solution) < 1e-5))
  expect_true(all(abs(c_solution - r_solution) < 1e-5))
})


# # speed comparison between R and Cpp
# library(microbenchmark)
# op <- microbenchmark(
#   cpp_code = riskParityPortfolio:::rpp_equality_constraints_iteration(Cmat, cvec, Qk, qk),
#   R_code = riskParityPortfolio:::rpp_equality_constraints_iteration_R(Cmat, cvec, Qk, qk),
#   QP_solver = qp_solution <- quadprog::solve.QP(Qk, -qk, Amat = t(Cmat), bvec = cvec, meq = 1),
#   times = 100)
# print(op)
# boxplot(op, main = "Time comparison [milliseconds]",
#         xlab = NULL, ylab = NULL,
#         unit = "ms", outline = FALSE, las = 2)



test_that("rpp_with_ineq_and_eq_constraints_iteration agree with solve.QP", {
  N <- sample(c(2:10), 1)
  Qk <- diag(c(1:N))
  qk <- rep(1, N)
  Cmat <- matrix(1, 1, N)
  cvec <- c(runif(1))
  Dmat <- matrix(0, N, N)
  diag(Dmat) <- rep(-1, N)
  dvec = c(rep(0, N))
  meq <- nrow(Cmat)
  Amat <- t(rbind(Cmat, -Dmat))
  bvec <- c(cvec, -dvec)

  # use the uniform portfolio as initial guess
  wk <- rep(1/N, N)
  c_solution <- riskParityPortfolio:::rpp_eq_and_ineq_constraints_iteration(
                                                                            Cmat, cvec, Dmat, dvec, Qk, qk, wk,
                                                                            rep(0, nrow(Dmat)), rep(0, nrow(Dmat)),
                                                                            rep(0, nrow(Cmat)), rep(0, nrow(Cmat)), 100, 1e-6)[[5]]
  qp_solution <- quadprog::solve.QP(Qk, -qk, Amat = Amat, bvec = bvec, meq = 1)$solution
  expect_true(all(abs(c_solution - qp_solution) < 1e-5))
  r_solution <- riskParityPortfolio:::rpp_eq_and_ineq_constraints_iteration_R(
                                                                            Cmat, cvec, Dmat, dvec, Qk, qk, wk,
                                                                            rep(0, nrow(Dmat)), rep(0, nrow(Dmat)),
                                                                            rep(0, nrow(Cmat)), rep(0, nrow(Cmat)))[[5]]
  expect_true(all(abs(c_solution - r_solution) < 1e-5))
})


# # speed comparison between R and Cpp
# library(microbenchmark)
# op <- microbenchmark(
#   cpp_code = riskParityPortfolio:::rpp_eq_and_ineq_constraints_iteration(Cmat, cvec, Dmat, dvec, Qk, qk, wk,
#                                                                          rep(0, nrow(Dmat)), rep(0, nrow(Dmat)),
#                                                                          rep(0, nrow(Cmat)), rep(0, nrow(Cmat))),
#   R_code = riskParityPortfolio:::rpp_eq_and_ineq_constraints_iteration_R(Cmat, cvec, Dmat, dvec, Qk, qk, wk,
#                                                                          rep(0, nrow(Dmat)), rep(0, nrow(Dmat)),
#                                                                          rep(0, nrow(Cmat)), rep(0, nrow(Cmat))),
#   QP_solver = quadprog::solve.QP(Qk, -qk, Amat = Amat, bvec = bvec, meq = 1),
#   times = 10)
# print(op)
# boxplot(op, main = "Time comparison [milliseconds]",
#         xlab = NULL, ylab = NULL,
#         unit = "ms", outline = FALSE, las = 2)



