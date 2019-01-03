context("Formulations")

# generate a random Sigma and mu
set.seed(123)
N <- 5
T <- 1000
V <- matrix(rnorm(N * T), T, N)
Sigma <- cov(V)

test_that("sca, alabama, and slsq give similar results in low dimensional problems", {
  formulations_list <- c("rc-double-index", "rc-over-b-double-index",
                         "rc-over-var vs b", "rc-over-var",
                         "rc-over-sd vs b-times-sd", "rc vs b-times-var",
                         "rc vs theta", "rc-over-b vs theta")
  for(formulation in formulations_list) {
    suppressWarnings(rpp_sca <- riskParityPortfolio(Sigma, method = "sca", formulation = formulation))
    suppressWarnings(rpp_alabama <- riskParityPortfolio(Sigma, method = "alabama", formulation = formulation))
    suppressWarnings(rpp_slsqp <- riskParityPortfolio(Sigma, method = "slsqp", formulation = formulation))
    expect_equal(rpp_sca$w, rpp_alabama$w, tolerance = 1e-3)
    expect_equal(rpp_sca$w, rpp_slsqp$w, tolerance = 1e-3)
  }
})
