#' Implements the risk parity portfolio for the case of diagonal Sigma
#' that satisfies the constraints sum(w) = 1 and w >= 0.
#'
#' @export
riskParityPortfolioDiagSigma <- function(Sigma) {
  w <- 1 / sqrt(diag(Sigma))
  w <- w / sum(w)
  return (list(w = w, r = w * (Sigma %*% w)))
}

#' Implements the risk parity portfolio using SCA and QP solver
#' @export
riskParityPortfolioSCA <- function(Sigma, w0 = NA, budget = TRUE,
                                   shortselling = FALSE,
                                   formulation = "double-index", gamma = .9,
                                   zeta = 1e-7, tau = NA, maxiter = 500,
                                   ftol = 1e-9, wtol = 1e-6) {
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

  if (formulation == "double-index") {
    # define the risk function for a double-index formulation
    fn <- function(w, Sigma, N) {
      wSw <- w * (Sigma %*% w)
      return (2 * (N * sum(wSw^2) - sum(wSw)^2))
    }
    g <- function(w, Sigma, N) {
      wSw <-  w * (Sigma %*% w)
      return (rep(wSw, times = N) - rep(wSw, each = N))
    }
    computeA <- compute_A_double_index
  } else if (formulation == "single-index") {
    # define the risk function for a single-index formulation
    fn <- function(w, Sigma, N) {
      wSw <- w * (Sigma %*% w)
      # return (sum((wSw / sum(wSw) - 1 / N) ^ 2))
      return (sum(wSw ^ 2) / (sum(wSw) ^ 2) - 1 / N)
    }
    g <- function(w, Sigma, N) {
      wSw <- w * (Sigma %*% w)
      return (wSw / sum(wSw) - 1 / N)
    }
    computeA <- compute_A_single_index_R
  }
  # compute and store objective function at the initial value
  fun_k <- fn(wk, Sigma, N)
  fun_seq <- c(fun_k)
  time_seq <- c(0)
  start_time <- Sys.time()
  for (k in 1:maxiter) {
    # auxiliary quantities
    Ak <- computeA(wk, Sigma, N)
    g_wk <- g(wk, Sigma, N)
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
    if ((werr < wtol) & k > 1) {
      break
    }
    # check convergence on objective function
    ferr <- abs(fun_next - fun_k) / max(1., abs(fun_next))
    if ((ferr < ftol) & k > 1) {
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
                                         shortselling = FALSE, use_gradient = TRUE,
                                         formulation = "double-index", method = "slsqp",
                                         maxiter = 500, ftol = 1e-9, wtol = 1e-6) {
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

  if (formulation == "double-index") {
    # define the risk for a double-index formulation
    fn <- function(w, Sigma, N) {
      wSw <- w * (Sigma %*% w)
      return (2 * (N * sum(wSw^2) - sum(wSw)^2))
    }

    if (use_gradient) {
      # define the gradient of the risk for a double-index formulation
      fn_grad <- function(w, Sigma, N) {
        wSw <- w * (Sigma %*% w)
        v <- N * wSw - sum(wSw)
        risk_grad <- 4 * (Sigma %*% (w * v) + (Sigma %*% w) * v)
        return (risk_grad)
      }
    } else {
      fn_grad <- NULL
    }
  } else if (formulation == "single-index") {
    # define the risk for a single-index formulation
    fn <- function(w, Sigma, N) {
      wSw <- w * (Sigma %*% w)
      # return (sum((wSw / sum(wSw) - 1 / N) ^ 2))
      return (sum(wSw ^ 2) / (sum(wSw) ^ 2) - 1 / N)
    }

    if (use_gradient) {
      # define the gradient of the risk for a single-index formulation
      fn_grad <- function(w, Sigma, N) {
        wSw <- w * (Sigma %*% w)
        sum_wSw <- sum(wSw)
        wSw_b <- wSw / sum_wSw - 1 / N
        v <- wSw_b - sum(wSw_b * wSw) / (sum_wSw ^ 2)
        risk_grad <- (2 / sum_wSw) * (Sigma %*% (w * v) + (Sigma %*% w) * v)
        return (risk_grad)
      }
    } else {
      fn_grad <- NULL
    }
  } else {
    stop("formulation ", formulation, " is not included.")
  }

  fun_seq <- c(fn(w0, Sigma, N))
  time_seq <- c(0)
  if (method == "alabama") {
    start_time <- Sys.time()
    res <- alabama::constrOptim.nl(w0, fn, fn_grad, hin = shortselling,
                                   hin.jac = shortselling.jac,
                                   heq = budget, heq.jac = budget.jac,
                                   Sigma = Sigma, N = N,
                                   control.outer = list(trace = FALSE,
                                                        itmax = maxiter))
    end_time <- Sys.time()
  } else if (method == "slsqp") {
    start_time <- Sys.time()
    res <- nloptr::slsqp(w0, fn, fn_grad, hin = shortselling,
                         hinjac = shortselling.jac,
                         heq = budget, heqjac = budget.jac,
                         Sigma = Sigma, N = N, control = list(xtol_rel = wtol,
                                                              ftol_rel = ftol))
    end_time <- Sys.time()
  }
  # save objective value and elapsed time
  time_seq <- c(time_seq, end_time - start_time)
  fun_seq <- c(fun_seq, res$value)
  w <- res$par
  return(list(w = w, r = w * (Sigma %*% w),
              obj_fun = fun_seq, elapsed_time = time_seq))
}
