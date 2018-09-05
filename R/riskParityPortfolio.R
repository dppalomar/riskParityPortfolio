library(stats)
#' Implements the risk parity portfolio using SCA and CVXR
#' @export
riskParityPortfolioCVX <- function(mu, Sigma, nu = 0, budget = TRUE,
                                   shortselling = FALSE, w0 = NA, gamma = .9,
                                   zeta = 1e-7, tau = NA, lambda = .5,
                                   maxiter = 500, ftol = 1e-5, wtol = 1e-5) {
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

  nll_seq <- c()
  fun_seq <- c()
  funk <- Inf
  for (k in 1:maxiter) {
    # auxiliary quantities
    AkT <- t(compute_A(wk, Sigma))
    Qk <- 2 * AkT %*% t(AkT) + tau * diag(N)
    qk <- 2 * AkT %*% g_16(wk, Sigma) - Qk %*% wk
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

  return(list(portfolio_weights = w_next,
              risk_contributions = w_next * (Sigma %*% w_next),
              negloglike = nll_seq, obj_fun = fun_seq))
}

#' @export
riskParityPortfolioGenSolver <- function(Sigma, w0 = NA) {
  if (any(is.na(w0))) {
    w0 <- 1 / sqrt(diag(Sigma))
    w0 <- w0 / sum(w0)
  }

  N <- nrow(Sigma)
  fn <- function(w, Sigma) {
    wSw <- w * (Sigma %*% w)
    risk <- 2 * (N * sum(wSw^2) - sum(wSw)^2)
    return (risk)
  }

  fn_grad <- function(w, Sigma) {
    wSw <- w * (Sigma %*% w)
    v <- (N * wSw - sum(wSw) * rep(1, N))
    risk_grad <- 4 * (Sigma %*% (w * v) + (Sigma %*% w) * v)
    return(risk_grad)
  }

  res <- constrOptim(w0, fn, fn_grad, ui = diag(N), ci = rep(0, N),
                     Sigma = Sigma, method = "BFGS")
  wopt <- res$par
  return(list(init_portfolio_weights = w0,
              portfolio_weights = wopt,
              risk_contributions = wopt * (Sigma %*% wopt)))
}
