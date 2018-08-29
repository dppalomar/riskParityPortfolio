# This module contains the code to implement the SCA
# solution for the risk parity portfolio design problem
# as presented in FengPalomar-TSP2015.pdf

riskParityPortfolioCVX <- function(mu, Sigma, nu = 0, shortselling = FALSE,
                                   w0 = NA, gamma = .9, zeta = .1, tau = .5,
                                   lambda = .5, maxiter = 500) {
  if (any(is.na(w0))) {
    wk <- 1 / diag(Sigma)
    wk <- wk / sum(wk)
  } else {
    wk <- w0
  }
  constraints <- list(sum(w) == 1)
  if (!shortselling) {
    constraints <- list(constraints, w >= 0)
  }

  w <- CVXR::Variable(length(mu))
  for k in (1:maxiter) {
    # auxiliary quantities
    AkT <- g_grad(wk, Sigma)
    gk <- t(g(wk, Sigma))
    Qk <- 2 * AkT %*% Ak + tau * I
    qk <- 2 * AkT %*% g(wk) - Qk %*% wk
    # build problem (39) as in Feng & Palomar TSP2015
    F <- negLogLikelihoodCVX(w, nu, mu, Sigma)
    P <- .5 * CVXR::quad_form(w, Qk) + t(w) %*% qk
    w_hat <- solve(CVXR::Problem(CVXR::Minimize(P + lambda * F),
                                 constraints = constraints))$getValue(w)
    w_next <- wk + gamma * (w_hat - wk)
    # update variables
    wk <- w_next
    gamma <- gamma * (1 - zeta * gamma)
  }

  return(portfolio_weights = w_next, negloglike = nll_seq, obj_fun = fun_seq)
}
