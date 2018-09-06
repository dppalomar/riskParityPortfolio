#' Implements the risk parity portfolio using SCA and CVXR
#' @export
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

  fun_seq <- c()
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

  return(list(portfolio_weights = w_next,
              risk_contributions = w_next * (Sigma %*% w_next),
              obj_fun = fun_seq))
}

#' Implements the risk parity portfolio using SCA and CVXR
#' @export
riskParityPortfolioQP <- function(Sigma, w0 = NA, gamma = .9,
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

  fun_seq <- c()
  funk <- Inf
  for (k in 1:maxiter) {
    # auxiliary quantities
    Ak <- computeACpp(wk, Sigma)
    risks <-  wk * (Sigma %*% wk)
    g_wk <- rep(risks, times = N) - rep(risks, each = N)
    Qk <- 2 * crossprod(Ak) + tau * diag(N)
    qk <- 2 * t(Ak) %*% g_wk - Qk %*% wk
    # build and solve problem (39) as in Feng & Palomar TSP2015
    w_hat <- quadprog::solve.QP(Qk, -qk, matrix(1, N, 1), 1, meq = 1)$solution
    w_next <- wk + gamma * (w_hat - wk)
    # save objective function values
    fun_next <- sum(g_wk * g_wk)
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

  return(list(portfolio_weights = w_next,
              risk_contributions = w_next * (Sigma %*% w_next),
              obj_fun = fun_seq))
}

#' Implements the risk parity portfolio using a general constrained
#' solver from the alabama package
#' @export
riskParityPortfolioGenSolver <- function(Sigma, w0 = NA, budget = TRUE,
                                         shortselling = FALSE) {
  N <- nrow(Sigma)

  if (any(is.na(w0))) {
    w0 <- 1 / sqrt(diag(Sigma))
    w0 <- w0 / sum(w0)
  }

  if (budget) {
    budget <- function(w, ...) {
      return (sum(w) - 1)
    }

    budget.jac <- function(w, ...) {
      return (matrix(1, 1, N))
    }
  } else {
    budget <- NULL
    budget.jac <- NULL
  }

  if (!shortselling) {
    shortselling <- function(w, ...) {
      return (w)
    }

    shortselling.jac <- function(w, ...) {
      return (diag(N))
    }
  } else {
    shortselling <- NULL
    shortselling.jac <- NULL
  }

  fn <- function(w, Sigma) {
    wSw <- w * (Sigma %*% w)
    risk <- 2 * (N * sum(wSw^2) - sum(wSw)^2)
    return (risk)
  }

  fn_grad <- function(w, Sigma) {
    wSw <- w * (Sigma %*% w)
    v <- (N * wSw - sum(wSw) * rep(1, N))
    risk_grad <- 4 * (Sigma %*% (w * v) + (Sigma %*% w) * v)
    return (risk_grad)
  }

  res <- alabama::constrOptim.nl(w0, fn, fn_grad,
                                 hin = shortselling, hin.jac = shortselling.jac,
                                 heq = budget, heq.jac = budget.jac,
                                 Sigma = Sigma, control.outer = list(trace = FALSE))
  wopt <- res$par
  return(list(portfolio_weights = wopt,
              risk_contributions = wopt * (Sigma %*% wopt),
              init_portfolio_weights = w0))
}
