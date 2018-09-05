#' @export
computeA <- function(w, Sigma) {
  N <- length(w)
  diagSw <- Sigma * w
  Sdiagw <- Sigma %*% diag(w)
  A <- matrix(NA, N^2, N)
  for (i in 1:N) {
    for (j in 1:N) {
      A[i + (j-1)*N, ] <-  diagSw[, i] + Sdiagw[, i] - (diagSw[, j] + Sdiagw[, j])
    }
  }
  return(A)
}

#' @export
compute_A <- function(w, Sigma) {
  N <- length(w)
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
  return(A)
}
