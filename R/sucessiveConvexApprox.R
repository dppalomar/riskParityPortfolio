# This module contains the code to implement the SCA
# problem described in FengPalomar-ICASSP2016.pdf

theta_update <- function(theta, x, w, g, gamma) {
  theta_hat <- sum(x * g(w, Sigma))
  return(theta + gamma * (theta_hat - theta))
}

w_update <- function(...) {
  return(NA)
}


d1 <- function(x, p = 0.05, e = 1e-4) {
  d1x <- rep(0, length(x))
  mask <- (abs(x) <= e)
  d1x[mask] <- abs(x) / (e * (p + e))
  d1x[!mask] <- 1 / ((abs(x) + p))
  d1x <- d1x / log(1 + 1 / p)
  return(d1x)
}


d2 <- function(x, p = 0.05, e = 1e-4) {
  d2x <- rep(0, length(x))
  mask <- (abs(x) <= e)
  d2x[mask] <- 1 / (e * (p + e))
  d2x[!mask] <- 1 / (abs(x) * (abs(x) + p))
  d2x <- d2x / (2 * log(1 + 1 / p))
  return(d2x)
}

g <- function(w, Sigma) {
}
