# This module contains the code to implement the SCA
# solution for the risk parity portfolio design problem
# as presented in FengPalomar-TSP2015.pdf

#' @export
riskParityPortfolioCVX <- function(mu, Sigma, nu = 0, budget = TRUE,
                                   shortselling = FALSE, w0 = NA, gamma = .9,
                                   zeta = 1e-7, tau = NA, lambda = .5,
                                   maxiter = 500, ftol = 1e-5, wtol = 1e-5) {
  N <- nrow(Sigma)
  if (any(is.na(w0))) {
    wk <- rep(1 / N, N)
  } else {
    wk <- w0
  }
  if (is.na(tau)) {
    tau <- .05 * sum(diag(Sigma)) / (2 * N)
  }
  # build constraints
  w <- CVXR::Variable(N)
  constraints <- list()
  if (budget) {
    constraints <- c(constraints, sum(w) == 1)
  }
  if (!shortselling) {
    constraints <- c(constraints, w >= 0)
  }

  nll_seq <- c()
  fun_seq <- c()
  funk <- Inf
  for (k in 1:maxiter) {
    # auxiliary quantities
    AkT <- t(g_grad(wk, Sigma))
    Qk <- 2 * AkT %*% t(AkT) + tau * diag(N)
    qk <- 2 * AkT %*% g(wk, Sigma) - Qk %*% wk
    # build and solve problem (39) as in Feng & Palomar TSP2015
    F <- CVXR::quad_form(w, Sigma) - nu * t(w) %*% mu
    P <- .5 * CVXR::quad_form(w, Qk) + t(w) %*% qk
    obj_fun <- CVXR::Minimize(P + lambda * F)
    prob <- CVXR::Problem(obj_fun, constraints = constraints)
    result <- solve(prob)
    w_hat <- result$getValue(w)
    w_next <- wk + gamma * (w_hat - wk)
    # save likelihood and objective function values
    nll_seq <- c(nll_seq, t(w_next) %*% Sigma %*% w_next - nu * t(w_next) %*% mu)
    fun_next <- result$value
    fun_seq <- c(fun_seq, fun_next)
    # check convergence on parameters
    werr <- norm(w_next - wk, "2") / max(1., norm(w_next, "2"))
    if (werr < wtol) {
      break
    }
    # check convergence on objective function
    ferr <- abs(fun_next - funk) / max(1., abs(fun_next))
    if (ferr < ftol) {
      break
    }
    # update variables
    wk <- w_next
    funk <- fun_next
    gamma <- gamma * (1 - zeta * gamma)
  }

  return(list(portfolio_weights = w_next, risk_contrib = w_next * (Sigma %*% w_next),
              negloglike = nll_seq, obj_fun = fun_seq))
}
