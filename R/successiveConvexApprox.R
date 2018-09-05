#' @export
g_16 <- function(w, Sigma) {
  N <- length(w)
  risks <-  w * (Sigma %*% w)
  return (rep(risks, times = N) - rep(risks, each = N))
}

#' @export
compute_A <- function(w, Sigma) {
  N <- length(w)
  g <- rep(NA, N^2)
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
