# This module contains the code to implement the SCA
# solution for the risk parity portfolio design problem
# as presented in FengPalomar-TSP2015.pdf

#' @export
riskParityPortfolioCVX <- function(mu, Sigma, nu = 0, shortselling = FALSE,
                                   w0 = NA, gamma = .9, zeta = 1e-7, tau = NA,
                                   lambda = .5, maxiter = 500, ftol = 1e-5,
                                   wtol = 1e-5) {
  if (any(is.na(w0))) {
    wk <- rep(1 / nrow(Sigma), nrow(Sigma))
  } else {
    wk <- w0
  }
  if (is.na(tau)) {
    tau <- .05 * sum(diag(Sigma)) / (2 * nrow(Sigma))
  }
  # build constraints
  w <- CVXR::Variable(nrow(Sigma))
  constraints <- list(sum(w) == 1)
  if (!shortselling) {
    constraints <- c(constraints, w >= 0)
  }

  nll_seq <- c()
  fun_seq <- c()
  funk <- Inf
  for (k in 1:maxiter) {
    # auxiliary quantities
    AkT <- t(g_grad(wk, Sigma))
    gk <- t(g(wk, Sigma))
    Qk <- 2 * AkT %*% t(AkT) + tau * diag(nrow(AkT))
    qk <- 2 * AkT %*% g(wk, Sigma) - Qk %*% wk
    # build and solve problem (39) as in Feng & Palomar TSP2015
    F <- negLogLikelihoodCVX(w, nu, mu, Sigma)
    P <- .5 * CVXR::quad_form(w, Qk) + t(w) %*% qk
    obj_fun <- CVXR::Minimize(P + lambda * F)
    prob <- CVXR::Problem(obj_fun, constraints = constraints)
    result <- solve(prob)
    w_hat <- result$getValue(w)
    w_next <- wk + gamma * (w_hat - wk)
    # save likelihood and objective function values
    nll_seq <- c(nll_seq, negLogLikelihood(w_next, nu, mu, Sigma))
    fun_next <- result$value
    fun_seq <- c(fun_seq, fun_next)
    # check convergence on parameters
    werr <- norm(w_next - wk, "2") / max(1., norm(w_next, "2"))
    if (werr < wtol) {
      wk <- w_next
      break
    }
    # check convergence on objective function
    ferr <- abs(fun_next - funk) / max(1., abs(fun_next))
    if (ferr < ftol) {
      wk <- w_next
      break
    }
    # update variables
    wk <- w_next
    funk <- fun_next
    gamma <- gamma * (1 - zeta * gamma)
  }

  return(list(portfolio_weights = w_next,
              negloglike = nll_seq, obj_fun = fun_seq))
}
