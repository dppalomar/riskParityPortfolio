# This module contains the code to implement the SCA
# problem described in FengPalomar-ICASSP2016.pdf

theta_update <- function(theta, x, w, g, gamma) {
  theta_hat <- sum(x * g(w))
  return(theta + gamma * (theta_hat - theta))
}

w_update <- function(...) {
  return(NA)
}


d1 <- function(x, p = 0.05, e = 1e-4) {
  if (abs(x) < e) {
    return abs(x) / (e * (p + e) * log(1 + 1 / p))
  } else {
    return 1 / ((abs(x) + p) * log(1 + 1 / p))
  }
}


d2 <- function(x, p = 0.05, e = 1e-4) {
  return(NA)
}
