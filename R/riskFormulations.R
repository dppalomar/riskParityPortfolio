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

R_rc_double_index <- function(w, Sigma, b = NA, r = w*(Sigma %*% w)) {
  N <- length(w)
  return (2*(N*sum(r^2) - sum(r)^2))
}

R_grad_rc_double_index <- function(w, Sigma, b = NA, Sigma_w = Sigma %*% w, r = w*Sigma_w) {
  N <- length(w)
  v <- N*r - sum(r)
  return (as.vector(4*(Sigma %*% (w*v) + Sigma_w*v)))
}

g_rc_double_index <- function(w, Sigma, b = NA, r = w*(Sigma %*% w)) {
  N <- length(w)
  return (rep(r, times = N) - rep(r, each = N))
}

A_rc_double_index <- function(w, Sigma, b = NA, Sigma_w = Sigma %*% w) {
  N <- length(w)
  Sigma_w <- as.vector(Sigma_w)
  Ut <- diag(Sigma_w) + Sigma*w
  return (matrix(rep(t(Ut), N), ncol = N, byrow = TRUE) - matrix(rep(Ut, each = N), ncol = N))
}


####################################################################
# Compute g, R, and A for the formulation "rc-over-b-double-index" #
####################################################################

R_rc_over_b_double_index <- function(w, Sigma, b, r = w*(Sigma %*% w)) {
  N <- length(w)
  rb <- r/b
  return (2*(N*sum(rb^2) - sum(rb)^2))
}

R_grad_rc_over_b_double_index <- function(w, Sigma, b, Sigma_w = Sigma %*% w, r = w*Sigma_w) {
  N <- length(w)
  rb <- r/b
  v <- N*rb - sum(rb)
  return (as.vector(4*(Sigma %*% (w*v) + Sigma_w*v)))
}

g_rc_over_b_double_index <- function(w, Sigma, b, r = w*(Sigma %*% w)) {
  N <- length(w)
  rb <- r/b
  return (rep(rb, times = N) - rep(rb, each = N))
}

A_rc_over_b_double_index <- function(w, Sigma, b, Sigma_w = Sigma %*% w) {
  N <- length(w)
  Sigma_w <- as.vector(Sigma_w)
  Ut <- diag(Sigma_w) + Sigma*w
  Utb <- Ut / b
  return (matrix(rep(t(Utb), N), ncol = N, byrow = TRUE) - matrix(rep(Utb, each = N), ncol = N))
}


##############################################################
# Compute g, R, and A for the formulation "rc-over-var-vs-b" #
##############################################################

g_rc_over_var_vs_b <- function(w, Sigma, b, r = w*(Sigma %*% w)) {
  return (as.vector(r/sum(r) - b))
}

R_rc_over_var_vs_b <- function(w, Sigma, b, r = w*(Sigma %*% w)) {
  return (sum((g_rc_over_var_vs_b(w, Sigma, b, r = r))^2))
}

R_grad_rc_over_var_vs_b <- function(w, Sigma, b, Sigma_w = Sigma %*% w, r = w*Sigma_w) {
  sum_r <- sum(r)
  r_sumr_b <- r/sum_r - b
  v <- r_sumr_b - sum(r_sumr_b*r)/sum_r
  return ((2/sum_r) * as.vector(Sigma %*% (w*v) + Sigma_w*v))
}

A_rc_over_var_vs_b <- function(w, Sigma, b = NA, Sigma_w = Sigma %*% w, r = w*Sigma_w) {
  sum_r <- sum(r)
  Sigma_w <- as.vector(Sigma_w)
  r <- as.vector(r)
  Ut <- diag(Sigma_w) + Sigma*w
  return (Ut/sum_r - 2/(sum_r^2) * r %o% Sigma_w)
}


#########################################################
# Compute g, R, and A for the formulation "rc-over-var" #
#########################################################

g_rc_over_var <- function(w, Sigma, b = NA, r = w*(Sigma %*% w)) {
  return (as.vector(r/sum(r)))
}

R_rc_over_var <- function(w, Sigma, b = NA, r = w*(Sigma %*% w)) {
  return (sum((g_rc_over_var(w, Sigma, r = r))^2))
}

R_grad_rc_over_var <- function(w, Sigma, b, Sigma_w = Sigma %*% w, r = w*Sigma_w) {
  sum_r <- sum(r)
  r_sumr <- r/sum_r
  v <- r_sumr - sum(r_sumr^2)
  return ((2/sum_r) * as.vector(Sigma %*% (w*v) + Sigma_w*v))
}

A_rc_over_var <- A_rc_over_var_vs_b


######################################################################
# Compute g, R, and A for the formulation "rc-over-sd-vs-b-times-sd" #
######################################################################

g_rc_over_sd_vs_b_times_sd <- function(w, Sigma, b, r = w*(Sigma %*% w)) {
  sqrt_sum_r <- sqrt(sum(r))
  return (r/sqrt_sum_r - b*sqrt_sum_r)
}

R_rc_over_sd_vs_b_times_sd <- function(w, Sigma, b, r = w*(Sigma %*% w)) {
  return (sum((g_rc_over_sd_vs_b_times_sd(w, Sigma, b, r = r))^2))
}

R_grad_rc_over_sd_vs_b_times_sd <- function(w, Sigma, b, Sigma_w = Sigma %*% w, r = w*Sigma_w) {
  sum_r <- sum(r)
  r_sumr_b <- r/sum_r - b
  v <- 2*r_sumr_b - sum(r_sumr_b^2)
  return (as.vector(Sigma %*% (w*v) + Sigma_w*v))
}

A_rc_over_sd_vs_b_times_sd <- function(w, Sigma, b, Sigma_w = Sigma %*% w, r = w*Sigma_w) {
  sum_r <- sum(r)
  Sigma_w <- as.vector(Sigma_w)
  r <- as.vector(r)
  Ut <- diag(Sigma_w) + Sigma * w
  A <- (Ut - (r/sum_r + b) %o% Sigma_w) / sqrt(sum_r)
}
