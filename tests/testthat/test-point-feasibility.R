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

