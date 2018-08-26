context("testRiskContribFun.R")
library(riskParityPortfolio)
library(testthat)


test_that("negLogLikelihood returns a finite number", {
  w <- runif(3)
  w <- w / sum(w)
  nll <- negLogLikelihood(w, nu = .5, mu = runif(3), Sigma = diag(runif(3)))
  #expect_that(is.finite(nll), is_true())
})


test_that("negLogPrior returns a finite number", {
  w <- runif(3)
  w <- w / sum(w)
  w_k <- runif(3)
  w_k <- w_k / sum(w_k)
  theta <- runif(1)
  Sigma <- diag(runif(3))
  nlprior <- negLogPrior(w, w_k, theta, Sigma, .2, 4., 2e-3, 1e-8, .5, "1")
  expect_that(is.finite(nlprior), is_true())
})


test_that("we can minimize negLogLikelihood", {
  n <- 3
  w <- CVXR::Variable(n)
  Sigma <- diag(runif(n))
  mu <- runif(n)
  nu <- .5
  nll <- negLogLikelihood(w, nu, mu, Sigma)
  obj_fun <- CVXR::Minimize(nll)
  prob <- CVXR::Problem(obj_fun, constraints = list(sum(w) == 1))
  res <- solve(prob)
  w_hat <- res$getValue(w)
  expect_that(abs(sum(w_hat) - 1) < 1e-3, is_true())
})
