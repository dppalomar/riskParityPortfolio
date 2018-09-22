# Implements the risk parity portfolio using SCA and CVXR
# this implementation is basically used for the sake of unit
# testing.
riskParityPortfolioCVX <- function(Sigma, w0 = NA, budget = TRUE,
                                   shortselling = FALSE, gamma = .9,
                                   zeta = 1e-7, tau = NA, maxiter = 500,
                                   ftol = 1e-5, wtol = 1e-5) {
  N <- nrow(Sigma)
  if (any(is.na(w0))) {
    wk <- 1 / sqrt(diag(Sigma))
    wk <- wk / sum(wk)
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

  funk <- Inf
  for (k in 1:maxiter) {
    # auxiliary quantities
    risks <-  wk * (Sigma %*% wk)
    g_wk <- rep(risks, times = N) - rep(risks, each = N)
    Ak <- compute_A(wk, Sigma)
    Qk <- 2 * crossprod(Ak) + tau * diag(N)
    qk <- 2 * t(Ak) %*% g_wk - Qk %*% wk
    # build and solve problem (39) as in Feng & Palomar TSP2015
    P <- .5 * CVXR::quad_form(w, Qk) + t(w) %*% qk
    obj_fun <- CVXR::Minimize(P)
    prob <- CVXR::Problem(obj_fun, constraints = constraints)
    result <- solve(prob)
    w_hat <- result$getValue(w)
    w_next <- wk + gamma * (w_hat - wk)
    # save the objective function values
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

  return(list(w = w_next, risk_contributions = w_next * (Sigma %*% w_next),
              obj_fun = fun_seq))
}
