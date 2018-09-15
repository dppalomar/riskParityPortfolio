#########################################################
# Compute g and R for the formulation "rc double-index" #
# The matrix A is computed in C++                       #
#########################################################
g_rc_double_index <- function(w, Sigma, N, r) {
  return (rep(r, times = N) - rep(r, each = N))
}

R_rc_double_index <- function(w, Sigma, N) {
  r <- w * (Sigma %*% w)
  return (2 * (N * sum(r ^ 2) - sum(r) ^ 2))
}

R_grad_rc_double_index <- function(w, Sigma, N) {
  Sigma_w <- Sigma %*% w
  r <- w * Sigma_w
  v <- N * r - sum(r)
  return (4 * (Sigma %*% (w * v) + Sigma_w * v))
}

##############################################################
# Compute g, R, and A for the formulation "rc-over-var vs b" #
##############################################################
g_rc_over_var_vs_b <- function(w, Sigma, r, b) {
  sum_r <- sum(r)
  return (r / sum_r - b)
}

R_rc_over_var_vs_b <- function(w, Sigma, N, b) {
  r <- w * (Sigma %*% w)
  return (sum((g_rc_over_var_vs_b(w, Sigma, r, b)) ^ 2))
}

R_grad_rc_over_var_vs_b <- function(w, Sigma, N, b) {
  Sigma_w <- Sigma %*% w
  r <- w * Sigma_w
  sum_r <- sum(r)
  r_b <- r / sum_r - b
  v <- r_b - sum(r_b * r) / sum_r
  return ((2 / sum_r) * (Sigma %*% (w * v) + Sigma_w * v))
}

A_rc_over_var_vs_b <- function(w, Sigma, N, r) {
  sum_r <- sum(r)
  Mat <- t(Sigma * w) + diag(as.vector(Sigma %*% w))
  inv_sum_r <- 1 / sum_r
  return (inv_sum_r * (Mat - inv_sum_r * matrix(t(r) %*% Mat, N, N, byrow = TRUE)))
}

######################################################################
# Compute g, R, and A for the formulation "rc-over-sd vs b-times-sd" #
######################################################################
g_rc_over_sd_vs_b_times_sd <- function(w, Sigma, r, b) {
  sqrt_sum_r <- sqrt(sum(r))
  return (r / sqrt_sum_r - b * sqrt_sum_r)
}

R_rc_over_sd_vs_b_times_sd <- function(w, Sigma, N, b) {
  r <- w * (Sigma %*% w)
  return(sum((g_rc_over_sd_vs_b_times_sd(w, Sigma, r, b)) ^ 2))
}

R_grad_rc_over_sd_vs_b_times_sd <- function(w, Sigma, N, b) {
  Sigma_w <- Sigma %*% w
  r <- w * Sigma_w
  sum_r <- sum(r)
  r_b <- r / sum_r - b
  v <- 2 * r_b - sum(r_b ^ 2)
  return (Sigma %*% (w * v) + Sigma_w * v)
}

A_rc_over_sd_vs_b_times_sd <- function(w, Sigma, N, r, b) {
  sum_r <- sum(r)
  inv_sum_r <- 1 / sum_r
  Mat <- t(Sigma * w) + diag(as.vector(Sigma %*% w))
  return(sqrt(inv_sum_r) * (Mat - .5 * matrix(t(r / inv_sum_r + b) %*% Mat,
                                              N, N, byrow = TRUE)))
}
