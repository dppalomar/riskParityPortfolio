context("Point feasibility")

set.seed(42)
N <- 5
w0 <- 1 + runif(N)
w_lb <- rep(.1, N)
w_ub <- rep(.9, N)

projectQP <- function(w0, w_lb, w_ub) {
  N <- length(w0)
  Dmat <- diag(N)
  dvec <- w0
  Amat <- cbind(rep(1, N), diag(N), -diag(N))
  bvec <- c(1, w_lb, -w_ub)
  meq <- 1
  return(quadprog::solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat,
                            bvec = bvec, meq = meq)$solution)
}

test_that("unitroot and QP solvers agree", {
  expect_equal(projectBudgetLineAndBox(w0, w_lb, w_ub),
               projectQP(w0, w_lb, w_ub))
})


test_that("isFeasiblePortfolio makes sense", {
  Cmat <- c(1, 1)
  cvec <- c(1)
  Dmat <- -diag(2)
  dvec <- c(0, 0)

  w <- c(.5, .5)
  expect_true(isFeasiblePortfolio(w, Cmat, cvec, Dmat, dvec))
  w <- c(-.5, -.3)
  expect_false(isFeasiblePortfolio(w, Cmat, cvec, Dmat, dvec))
})

