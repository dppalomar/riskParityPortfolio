context("Risk Parity with Constraints")
# tests with data from the pyrb package:
# https://github.com/jcrichard/pyrb/blob/master/notebooks/ConstrainedRiskBudgeting.ipynb
vol <- c(0.05, 0.05, 0.07, 0.1, 0.15, 0.15, 0.15, 0.18)
Corr <- rbind(c(100,  80,  60, -20, -10, -20, -20, -20),
              c( 80, 100,  40, -20, -20, -10, -20, -20),
              c( 60,  40, 100,  50,  30,  20,  20,  30),
              c(-20, -20,  50, 100,  60,  60,  50,  60),
              c(-10, -20,  30,  60, 100,  90,  70,  70),
              c(-20, -10,  20,  60,  90, 100,  60,  70),
              c(-20, -20,  20,  50,  70,  60, 100,  70),
              c(-20, -20,  30,  60,  70,  70,  70, 100)) / 100
Sigma <- Corr * (vol %o% vol)

test_that("unconstrained example", {
    answer <- c(26.8306, 28.6769, 11.4095, 9.7985, 5.6135, 5.9029, 6.656, 5.1121) / 100
    rpp <- riskParityPortfolio(Sigma)
    expect_true(max(abs(answer - rpp$w)) < 1e-4)
    expect_true(max(abs(rpp$relative_risk_contribution - 0.125)) < 1e-4)
})


test_that("constrained example", {
    Dmat <- matrix(c(rep(0, 4), rep(-1, 4)), nrow = 1)
    dvec <- c(-0.3)
    rpp <- riskParityPortfolio(Sigma, method = "sca",
                               Dmat = Dmat, dvec = dvec)
    expect_true(abs(sum(rpp$w) - 1) < 1e-4)
    expect_true((Dmat %*% rpp$w - dvec) < 1e-5)
})


test_that("another constrained example", {
    Dmat <- matrix(0, 2, 8)
    Dmat[1, ] <- c(rep(0, 4), rep(-1, 4))
    Dmat[2, ] <- c(1, -1, 0, 0, 1, -1, 0, 0)
    dvec <- c(-0.3, -0.05)
    rpp <- riskParityPortfolio(Sigma, method = "sca",
                               Dmat = Dmat, dvec = dvec)
    expect_true(abs(sum(rpp$w) - 1) < 1e-4)
    expect_true(all(-Dmat %*% rpp$w > -(dvec + 1e-6)))
})
