#' function used for sanity-check / unit testing
compute_A_dan_version <- function(w, Sigma, N) {
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

#' Computes the A matrix for the formulation "rc-over-var vs b"
#' Also for sanity-check / unit testing
#' @export
compute_A_double_index_R <- function(w, Sigma, N) {
  A <- Matrix::Matrix(0, N ^ 2, N, sparse = TRUE)
  for (i in 1:(N-1)) {
    Mi <- Matrix::Matrix(0, N, N, sparse = TRUE)
    Mi[i, ] <- Sigma[i, ]
    for (j in (i+1):N) {
      Mj <- Matrix::Matrix(0, N, N, sparse = TRUE)
      Mj[j, ] <- Sigma[j, ]
      A[i + (j - 1) * N, ] <- (Mi + t(Mi) - Mj - t(Mj)) %*% w
      A[j + (i - 1) * N, ] <- - A[i + (j - 1) * N, ]
    }
  }
  return (A)
}
