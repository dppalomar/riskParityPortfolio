# This module contains the code to implement the SCA
# problem described in FengPalomar-ICASSP2016.pdf

library(CVXR)

#' Update step for the average RCs of the selected assets
#'
#' @param theta the avg RCs at the previous iteration
#' @param w the portfolio weights at the previous iteration
#' @param g function to represent the risk controbution
#' @param gamma the learning rate
theta_update <- function(theta_k, w_k, g, gamma) {
  rho_sq <- rho(w_k) ^ 2
  x <- rho_sq / sum(rho_sq ^ 2)
  theta_hat <- sum(x * g(w_k, Sigma))
  return(theta_k + gamma * (theta_hat - theta_k))
}

w_update <- function(w_k, theta_k, nu, gamma, l1, l2,
                     mu, Sigma, tau, type = "1") {
  w <- Variable(length(w_k))
  F <- negLogLikelihood(...)


  obj_fun <- Minimize(F + l1 * D + l2 * P + tau * sum((w - w_k) ^ 2))

  return(w + gamma * (w_hat - w))
}


negLogLikelihood <- function(w, nu, mu, Sigma) {
  return(t(w) %*% (Sigma %*% w - nu * mu))
}


negLogPrior <- function(l1 = .1, l2 = 4, type = "1") {
  if (type == "1") {
    D <- norm(d1(w_k) * w, type = type)
  } else if (type == "2") {
    D <- sum(d2(w_k) * (w ^ 2))
  } else {
    stop("type is not implemented")
  }
}


d1 <- function(x, p = 2e-3, e = 1e-8) {
  d1x <- rep(0, length(x))
  abs_x <- abs(x)
  mask <- (abs_x <= e)
  not_mask <- !mask
  d1x[mask] <- abs_x[mask] / (e * (p + e))
  d1x[not_mask] <- 1 / (abs_x[not_mask] + p)
  d1x <- d1x / log(1 + 1 / p)
  return (d1x)
}


d2 <- function(x, p = 2e-3, e = 1e-8) {
  d2x <- rep(1 / (e * (p + e)), length(x))
  abs_x <- abs(x)
  mask <- (abs_x > e)
  d2x[mask] <- 1 / (abs_x[mask] * (abs_x[mask] + p))
  d2x <- d2x / (2 * log(1 + 1 / p))
  return (d2x)
}

g <- function(w, Sigma) {
  return (w * (Sigma %*% w))
}

g_grad <- function(w, Sigma) {
  return (2 * Sigma %*% w)
}


#' Approximation for the lp-norm, 0 < p < 1
rho <- function(x, p, e) {
  val <- rep(0, length(x))
  abs_x <- abs(x)
  mask <- (abs_x <= e)
  not_mask <- !mask
  val[mask] <- x[mask] ^ 2 / (2 * e * (p + e))
  val[not_mask] <- log(1 + abs_x[not_mask] / p) - log(1 + e / p) + e / (2 * (p + e))
  val <- val / log(1 + 1 / p)
  return (val)
}

#' Gradient of the approximation for the lp-norm, 0 < p < 1, w.r.t. to x
rho_grad <- function(x, p, e) {
  val <- rep(0, length(x))
  abs_x <- abs(x)
  mask <- (abs_x <= e)
  not_mask <- !mask
  val[mask] <- x[mask] / (e * (p + e))
  val[not_mask] <- sign(x[not_mask]) / (p + abs_x[not_mask])
  val <- val / log(1 + 1 / p)
  return (val)
}
