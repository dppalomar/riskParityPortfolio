# Compute g and R for the formulation "rc double-index"
# The matrix A is computed in C++
g_rc_double_index <- function(w, Sigma, N, ...) {
  kwargs <- list(...)
  r <- kwargs$r
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
  risk_grad <- 4 * (Sigma %*% (w * v) + Sigma_w * v)
  return (risk_grad)
}

# Compute g, R, and A for the formulation "rc-over-var vs b"
g_rc_over_var_vs_b <- function(w, Sigma, N, ...) {
  kwargs <- list(...)
  r <- kwargs$r
  sum_r <- sum(r)
  return (r / sum_r - 1 / N)
}

R_rc_over_var_vs_b <- function(w, Sigma, N) {
  r <- w * (Sigma %*% w)
  return (sum((r / sum(r) - 1 / N) ^ 2))
}

R_grad_rc_over_var_vs_b <- function(w, Sigma, N) {
  Sigma_w <- Sigma %*% w
  r <- w * Sigma_w
  sum_r <- sum(r)
  r_b <- r / sum_r - 1 / N
  v <- r_b - sum(r_b * r) / sum_r
  risk_grad <- (2 / sum_r) * (Sigma %*% (w * v) + Sigma_w * v)
  return (risk_grad)
}

A_rc_over_var_vs_b <- function(w, Sigma, N, ...) {
  kwargs <- list(...)
  r <- kwargs$r
  sum_r <- sum(r)
  Mat <- t(Sigma * w) + diag(as.vector(Sigma %*% w))
  inv_sum_r <- 1 / sum_r
  A <- inv_sum_r * (Mat - inv_sum_r) * matrix(t(r) %*% Mat, N, N, byrow = TRUE)
  return (A)
}
