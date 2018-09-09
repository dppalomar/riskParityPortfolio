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

#' Implements the risk parity portfolio using SCA and QP solver
#' @export
riskParityPortfolioQP <- function(Sigma, w0 = NA, budget = TRUE,
                                  shortselling = FALSE, gamma = .9,
                                  zeta = 1e-7, tau = NA, maxiter = 500,
                                  ftol = 1e-9, wtol = 1e-9) {
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

  if (budget & (!shortselling)) {
    Amat <- cbind(matrix(1, N, 1), diag(N))
    bvec <- c(1, rep(0, N))
    meq <- 1
  } else if (budget) {
    Amat <- matrix(1, N, 1)
    bvec <- 1
    meq <- 1
  } else if (!shortselling) {
    Amat <- diag(N)
    bvec <- rep(0, N)
    meq <- 0
  }
  # compute and store objective function at the initial value
  fn <- function(w, Sigma, N) {
    wSw <- w * (Sigma %*% w)
    return(2 * (N * sum(wSw^2) - sum(wSw)^2))
  }
  fun_k <- fn(wk, Sigma, N)
  fun_seq <- c(fun_k)
  time_seq <- c(0)
  start_time <- Sys.time()
  for (k in 1:maxiter) {
    # auxiliary quantities
    Ak <- computeACpp(wk, Sigma)
    risks <-  wk * (Sigma %*% wk)
    g_wk <- rep(risks, times = N) - rep(risks, each = N)
    Qk <- 2 * crossprod(Ak) + tau * diag(N)
    qk <- 2 * t(Ak) %*% g_wk - Qk %*% wk
    # build and solve problem (39) as in Feng & Palomar TSP2015
    w_hat <- quadprog::solve.QP(Qk, -qk, Amat = Amat,
                                bvec = bvec, meq = meq)$solution
    w_next <- wk + gamma * (w_hat - wk)
    # save objective function values and elapsed time
    end_time <- Sys.time()
    time_seq <- c(time_seq, end_time - start_time)
    fun_next <- fn(w_next, Sigma, N)
    fun_seq <- c(fun_seq, fun_next)
    # check convergence on parameters
    werr <- norm(w_next - wk, "2") / max(1., norm(w_next, "2"))
    if (werr < wtol) {
      break
    }
    # check convergence on objective function
    ferr <- abs(fun_next - fun_k) / max(1., abs(fun_next))
    if (ferr < ftol) {
      break
    }
    # update variables
    wk <- w_next
    fun_k <- fun_next
    gamma <- gamma * (1 - zeta * gamma)
  }

  return(list(w = w_next, risk_contributions = w_next * (Sigma %*% w_next),
              obj_fun = fun_seq, elapsed_time = time_seq))
}

#' Implements the risk parity portfolio using a general constrained
#' solver from the alabama package
#' @export
riskParityPortfolioGenSolver <- function(Sigma, w0 = NA, budget = TRUE,
                                         shortselling = FALSE, maxiter = 500,
                                         ftol = 1e-9, wtol = 1e-9) {
  N <- nrow(Sigma)

  if (any(is.na(w0))) {
    wk <- 1 / sqrt(diag(Sigma))
    wk <- wk / sum(wk)
  } else {
    wk <- w0
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

  fn <- function(w, Sigma, N) {
    wSw <- w * (Sigma %*% w)
    return(2 * (N * sum(wSw^2) - sum(wSw)^2))
  }

  fn_grad <- function(w, Sigma, N) {
    wSw <- w * (Sigma %*% w)
    v <- (N * wSw - sum(wSw) * rep(1, N))
    risk_grad <- 4 * (Sigma %*% (w * v) + (Sigma %*% w) * v)
    return (risk_grad)
  }

  fun_k <- Inf
  fun_seq <- c(fn(wk, Sigma, N))
  time_seq <- c(0)
  start_time <- Sys.time()
  for (i in 1:maxiter) {
    res <- alabama::constrOptim.nl(wk, fn, fn_grad, hin = shortselling,
                                   hin.jac = shortselling.jac,
                                   heq = budget, heq.jac = budget.jac,
                                   Sigma = Sigma, N = N,
                                   control.outer = list(trace = FALSE, itmax = 1))
    # save objective value and elapsed time
    end_time <- Sys.time()
    time_seq <- c(time_seq, end_time - start_time)
    fun_next <- res$value
    fun_seq <- c(fun_seq, fun_next)
    # check convergence on parameters
    w_next <- res$par
    werr <- norm(w_next - wk, "2") / max(1., norm(w_next, "2"))
    if (werr < wtol) {
      break
    }
    # check convergence on objective function
    ferr <- abs(fun_next - fun_k) / max(1., abs(fun_next))
    if (ferr < ftol) {
      break
    }
    # update variables
    wk <- w_next
    fun_k <- fun_next
  }
  return(list(w = w_next, risk_contributions = w_next * (Sigma %*% w_next),
              obj_fun = fun_seq, elapsed_time = time_seq))
}
