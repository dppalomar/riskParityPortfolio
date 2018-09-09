#' @export
compute_A_double_index_R <- function(w, Sigma, N) {
  A <- matrix(NA, N^2, N)
  for (i in 1:N) {
    Mi <- matrix(0, N, N)
    Mi[i, ] <- Sigma[i, ]
    for (j in 1:N) {
      Mj <- matrix(0, N, N)
      Mj[j, ] <- Sigma[j, ]
      A[i + (j-1)*N, ] <- (Mi + t(Mi) - Mj - t(Mj)) %*% w
    }
  }
  return (A)
}

#' @export
compute_A_single_index_R <- function(w, Sigma, N) {
  wSw <- w * (Sigma %*% w)
  sum_wSw <- sum(wSw)
  Mat <- t(Sigma * w) + diag(as.vector(Sigma %*% w))
  inv_sum_wSw <- 1 / sum_wSw
  A <- inv_sum_wSw * (Mat - inv_sum_wSw) * matrix(t(wSw) %*% Mat, N, N, byrow = TRUE)
  return (A)
}
