# This module contains the code to implement the SCA
# problem described in FengPalomar-ICASSP2016.pdf

#' Update step for the average RCs of the selected assets
#'
#' @param theta the avg RCs at the previous iteration
#' @param w the portfolio weights at the previous iteration
#' @param g function to represent the risk controbution
#' @param gamma the learning rate
#'
#' @export
theta_update <- function(theta_k, w_k, Sigma, g, p, e, gamma) {
  rho_sq <- rho(w_k, p, e) ^ 2
  x <- rho_sq / sum(rho_sq ^ 2)
  theta_hat <- sum(x * g(w_k, Sigma))
  return(theta_k + gamma * (theta_hat - theta_k))
}


#' @export
w_update <- function(w_k, theta_k, nu, gamma, l1, l2,
                     mu, Sigma, tau, p, e, type) {
  w <- CVXR::Variable(length(w_k))
  nll <- negLogLikelihoodCVX(w, nu, mu, Sigma)
  nlprior <- negLogPrior(w, w_k, theta_k, Sigma, l1, l2, p, e, tau, type)
  obj_fun <- CVXR::Minimize(nll + nlprior)
  prob <- CVXR::Problem(obj_fun, constraints = list(w >= 0, sum(w) == 1))
  result <- solve(prob)
  w_hat <- result$getValue(w)
  return(w_k + gamma * (w_hat - w_k))
}


#' @export
negLogLikelihoodCVX <- function(w, nu, mu, Sigma) {
  return (CVXR::quad_form(w, Sigma) - nu * t(w) %*% mu)
}


#' @export
negLogLikelihood <- function(w, nu, mu, Sigma) {
  return (t(w) %*% (Sigma %*% w - nu * mu))
}


#' @export
negLogPrior <- function(w, w_k, theta, Sigma, l1, l2, p, e, tau, type) {
  if (type == "1") {
    D <- sum(abs(d1(w_k, p, e) * w))
  } else if (type == "2") {
    D <- sum(d2(w_k, p, e) * (w * w))
  } else {
    stop("type is not implemented")
  }

  P <- sum((g_tilde(w_k, theta, p, e, Sigma) +
            sum(g_tilde_grad(w_k, theta, p, e, Sigma) * (w - w_k))) ^ 2)

  return (l1 * D + l2 * P + tau * sum((w - w_k) ^ 2))
}


#' @export
d1 <- function(x, p, e) {
  d1x <- rep(0, length(x))
  abs_x <- abs(x)
  mask <- (abs_x <= e)
  not_mask <- !mask
  d1x[mask] <- abs_x[mask] / (e * (p + e))
  d1x[not_mask] <- 1 / (abs_x[not_mask] + p)
  d1x <- d1x / log(1 + 1 / p)
  return (d1x)
}


#' @export
d2 <- function(x, p, e) {
  d2x <- rep(1 / (e * (p + e)), length(x))
  abs_x <- abs(x)
  mask <- (abs_x > e)
  d2x[mask] <- 1 / (abs_x[mask] * (abs_x[mask] + p))
  d2x <- d2x / (2 * log(1 + 1 / p))
  return (d2x)
}


#' @export
g <- function(w, Sigma) {
  return (w * (Sigma %*% w))
}


#' @export
g_grad <- function(w, Sigma) {
  # Jacobian of the Hadamard product
  return (diag(as.vector(Sigma %*% w)) + diag(as.vector(w)) %*% Sigma)
}


#' @export
g_tilde <- function(w, theta, p, e, Sigma) {
  return ((g(w, Sigma) - theta) * rho(w, p, e))
}


#' @export
g_tilde_grad <- function(w, theta, p, e, Sigma) {
  return(rho(w, p, e) * g_grad(w, Sigma) +
         (g(w, Sigma) - theta) * rho_grad(w, p, e))
}


#' Approximation for the lp-norm, 0 < p < 1
#' @export
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
#' @export
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
