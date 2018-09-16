# This file is meant to have the implementation of various risk concentration
# functions and its quantities such as gradient, Jacobian, etc, used by the
# successive convex approximation algorithm and by a general solver available
# in the riskParityPortfolio package.
# We use the following definitions of functions in this file:
#   - functions starting with the letter g are assumed to compute the risk
#     concentration vector.
#   - functions starting with the letter R followed by the name of the
#     formulation are assumed to compute the risk concentration. Most of them
#     are computed as sum(g ^ 2) unless a more efficient way is available.
#   - functions starting with "R_grad" are assumed to compute the gradient of
#     the risk concentration wrt the portfolio weights.
#   - functions starting with the letter A are assumed to compute the Jacobian
#     of the risk concentration vector wrt the portfolio weights.

#############################################################
# Compute g, R, and A for the formulation "rc-double-index" #
#############################################################

#' Risk concentration vector for the formulation rc-double-index
#'
#' @param w portfolio weights
#' @param Sigma covariance or correlation matrix
#' @param N number of stocks (length of w)
#' @param r the quantity w * (Sigma %*% w)
#' @return g risk concentration vector
g_rc_double_index <- function(w, Sigma, N, r) {
  return (rep(r, times = N) - rep(r, each = N))
}

#' Objective function for the formulation "rc-double-index"
#'
#' @param w portfolio weights
#' @param Sigma covariance or correlation matrix
#' @param N number of stocks (length of w)
#' @return R scalar valued objective function
R_rc_double_index <- function(w, Sigma, N) {
  r <- w * (Sigma %*% w)
  return (2 * (N * sum(r ^ 2) - sum(r) ^ 2))
}

#' Gradient of the objective function wrt the portfolio weights for the
#' formulation "rc-double-index"
#'
#' @param w portfolio weights
#' @param Sigma covariance or correlation matrix
#' @param N number of stocks (length of w)
#' @return grad_R the gradient of the obejctive function wrt to w
R_grad_rc_double_index <- function(w, Sigma, N) {
  Sigma_w <- Sigma %*% w
  r <- w * Sigma_w
  v <- N * r - sum(r)
  return (4 * (Sigma %*% (w * v) + Sigma_w * v))
}

#' Jacobian of the risk concentration vector wrt the portfolio weights for the
#' formulation "rc-double-index"
#'
#' @param w portfolio weights
#' @param Sigma covariance or correlation matrix
#' @param N number of stocks (length of w)
#' @param r the quantity w * Sigma_w
#' @param Sigma_w the quantity Sigma %*% w
#' @return Jg the jacobian of the risk concentration vector wrt to w
A_rc_double_index <- function(w, Sigma, N, r, Sigma_w) {
  Ut <- diag(Sigma_w) + Sigma * w
  return(matrix(rep(t(Ut), N), ncol = N, byrow = TRUE) -
         matrix(rep(Ut, each = N), ncol = N))
}


##############################################################
# Compute g, R, and A for the formulation "rc-over-var-vs-b" #
##############################################################

#' Risk vector for the formulation "rc-over-var-vs-b"
#'
#' @param w portfolio weights
#' @param Sigma covariance or correlation matrix
#' @param N number of stocks (length of w)
#' @param r the quantity w * (Sigma %*% w)
#' @param b budget vector
#' @return g risk concentration vector
g_rc_over_var_vs_b <- function(w, Sigma, r, b) {
  sum_r <- sum(r)
  return (r / sum_r - b)
}

#' Objective function for the formulation "rc-over-var-vs-b"
#'
#' @param w portfolio weights
#' @param Sigma covariance or correlation matrix
#' @param N number of stocks (length of w)
#' @param b budget vector
#' @return R scalar valued objective function
R_rc_over_var_vs_b <- function(w, Sigma, N, b) {
  r <- w * (Sigma %*% w)
  return (sum((g_rc_over_var_vs_b(w, Sigma, r, b)) ^ 2))
}

#' Gradient of the objective function wrt the portfolio weights for the
#' formulation "rc-over-var-vs-b"
#'
#' @param w portfolio weights
#' @param Sigma covariance or correlation matrix
#' @param N number of stocks (length of w)
#' @param b budget vector
#' @return grad_R the gradient of the obejctive function wrt to w
R_grad_rc_over_var_vs_b <- function(w, Sigma, N, b) {
  Sigma_w <- Sigma %*% w
  r <- w * Sigma_w
  sum_r <- sum(r)
  r_b <- r / sum_r - b
  v <- r_b - sum(r_b * r) / sum_r
  return ((2 / sum_r) * (Sigma %*% (w * v) + Sigma_w * v))
}

#' Jacobian of the risk concentration vector wrt the portfolio weights for the
#' formulation "rc-over-var-vs-b"
#'
#' @param w portfolio weights
#' @param Sigma covariance or correlation matrix
#' @param N number of stocks (length of w)
#' @param r the quantity w * Sigma_w
#' @param Sigma_w the quantity Sigma %*% w
#' @return Jg the jacobian of the risk concentration vector wrt to w
A_rc_over_var_vs_b <- function(w, Sigma, N, r, Sigma_w) {
  sum_r <- sum(r)
  Ut <- diag(Sigma_w) + Sigma * w
  retun(Ut / sum_r - 2 / (sum_r^2) * r %o% Sigma_w)
}

######################################################################
# Compute g, R, and A for the formulation "rc-over-sd-vs-b-times-sd" #
######################################################################

#' Risk vector for the formulation "rc-over-sd-vs-b-times-sd"
#'
#' @param w portfolio weights
#' @param Sigma covariance or correlation matrix
#' @param N number of stocks (length of w)
#' @param r the quantity w * (Sigma %*% w)
#' @param b budget vector
#' @return g risk concentration vector
g_rc_over_sd_vs_b_times_sd <- function(w, Sigma, r, b) {
  sqrt_sum_r <- sqrt(sum(r))
  return (r / sqrt_sum_r - b * sqrt_sum_r)
}

#' Objective function for the formulation "rc-over-sd-vs-b-times-sd"
#'
#' @param w portfolio weights
#' @param Sigma covariance or correlation matrix
#' @param N number of stocks (length of w)
#' @param b budget vector
#' @return R scalar valued objective function
R_rc_over_sd_vs_b_times_sd <- function(w, Sigma, N, b) {
  r <- w * (Sigma %*% w)
  return(sum((g_rc_over_sd_vs_b_times_sd(w, Sigma, r, b)) ^ 2))
}

#' Gradient of the objective function wrt the portfolio weights for the
#' formulation "rc-over-sd-vs-b-times-sd"
#'
#' @param w portfolio weights
#' @param Sigma covariance or correlation matrix
#' @param N number of stocks (length of w)
#' @param b budget vector
#' @return grad_R the gradient of the obejctive function wrt to w
R_grad_rc_over_sd_vs_b_times_sd <- function(w, Sigma, N, b) {
  Sigma_w <- Sigma %*% w
  r <- w * Sigma_w
  sum_r <- sum(r)
  r_b <- r / sum_r - b
  v <- 2 * r_b - sum(r_b ^ 2)
  return (Sigma %*% (w * v) + Sigma_w * v)
}

#' Jacobian of the risk concentration vector wrt the portfolio weights for the
#' formulation "rc-over-sd-vs-b-times-sd"
#'
#' @param w portfolio weights
#' @param Sigma covariance or correlation matrix
#' @param N number of stocks (length of w)
#' @param r the quantity w * Sigma_w
#' @param Sigma_w the quantity Sigma %*% w
#' @param b budget vector
#' @return Jg the jacobian of the risk concentration vector wrt to w
A_rc_over_sd_vs_b_times_sd <- function(w, Sigma, N, r, Sigma_w, b) {
  sum_r <- sum(r)
  Ut <- diag(Sigma_w) + Sigma * w
  A <- (Ut - (r/sum_r + b) %o% Sigma_w) / sqrt(sum_r)
}
