#' Risk parity portfolio optimization for the case of diagonal Sigma
#' that satisfies the constraints sum(w) = 1 and w >= 0.
#'
#' @param Sigma covariance or correlation matrix
#' @param b budget vector
#' @return w optimal portfolio vector
#' @return risk_contribution the risk contribution of every asset
#' @export
riskParityPortfolioDiagSigma <- function(Sigma, b = rep(1/nrow(Sigma), nrow(Sigma))) {
  w <- sqrt(b) / sqrt(diag(Sigma))
  w <- w / sum(w)
  return (list(w = w,
               risk_contribution = as.vector(w * (Sigma %*% w))))
}


#' Risk parity portfolio optimization using successive convex approximation (SCA)
#' and a quadratic programming (QP) solver.
#'
#' @param Sigma covariance or correlation matrix
#' @param b budget vector
#' @param budget boolean indicating whether to consider sum(w) = 1 as a
#'        constraint
#' @param shortselling boolean indicating whether to allow short-selling, i.e.,
#'        w < 0
#' @param formulation string indicating the formulation to be used for the risk
#'        parity optimization problem. It must be one of: "rc-double-index",
#'        "rc-over-b-double-index", "rc-over-var vs b", "rc-over-var",
#'        "rc-over-sd vs b-times-sd", or "rc vs b-times-var"
#' @export
riskParityPortfolioSCA <- function(Sigma, b = rep(1/nrow(Sigma), nrow(Sigma)),
                                   mu = NA, lambda = 1e-4, budget = TRUE,
                                   shortselling = FALSE,
                                   formulation = "rc-over-var vs b", w0 = NA,
                                   theta0 = NA, gamma = .9, zeta = 1e-7, tau = NA,
                                   maxiter = 500, ftol = 1e-9, wtol = 1e-6) {
  N <- nrow(Sigma)

  if (anyNA(w0))
    w0 <- riskParityPortfolioDiagSigma(Sigma, b)$w

  has_theta <- grepl("theta", formulation)
  if (has_theta) {
    if (is.na(theta0))
      theta0 <- mean(w0 * (Sigma %*% w0))
    w0 <- as.vector(c(w0, theta0))
  }

  if (is.na(tau))
    tau <- .05 * sum(diag(Sigma)) / (2*N)

  if (has_theta) {
    if (budget & (!shortselling)) {
      Amat <- cbind(rbind(matrix(1, N, 1), 0), diag(c(rep(1, N), 0)))
      bvec <- c(1, rep(0, N+1))
      meq <- 1
    } else if (budget) {
      Amat <- rbind(matrix(1, N, 1), 0)
      bvec <- 1
      meq <- 1
    } else if (!shortselling) {
      Amat <- diag(c(rep(1, N), 0))
      bvec <- rep(0, N+1)
      meq <- 0
    }
  } else {
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
  }

  switch(formulation,
         "rc-double-index" = {
           R <- R_rc_double_index
           g <- g_rc_double_index
           A <- A_rc_double_index
         },
         "rc-over-b-double-index" = {
           R <- R_rc_over_b_double_index
           g <- g_rc_over_b_double_index
           A <- A_rc_over_b_double_index
         },
         "rc-over-var vs b" = {
           R <- R_rc_over_var_vs_b
           g <- g_rc_over_var_vs_b
           A <- A_rc_over_var_vs_b
         },
         "rc-over-var" = {
           R <- R_rc_over_var
           g <- g_rc_over_var
           A <- A_rc_over_var
         },
         "rc-over-sd vs b-times-sd" = {
           R <- R_rc_over_sd_vs_b_times_sd
           g <- g_rc_over_sd_vs_b_times_sd
           A <- A_rc_over_sd_vs_b_times_sd
         },
         "rc vs b-times-var"  = {
           R <- R_rc_vs_b_times_var
           g <- g_rc_vs_b_times_var
           A <- A_rc_vs_b_times_var
         },
         "rc vs theta" = {
           R <- R_rc_vs_theta
           g <- g_rc_vs_theta
           A <- A_rc_vs_theta
         },
         "rc-over-b vs theta" = {
           R <- R_rc_over_b_vs_theta
           g <- g_rc_over_b_vs_theta
           A <- A_rc_over_b_vs_theta
         },
         stop("formulation ", formulation, " is not included.")
  )

  # compute and store objective function at the initial value
  wk <- w0
  fun_k <- R(wk, Sigma, b)
  fun_seq <- c(fun_k)
  time_seq <- c(0)

  if (has_theta)
    tauI <- diag(rep(tau, N + 1))
  else
    tauI <- diag(rep(tau, N))

  start_time <- proc.time()[3]
  for (k in 1:maxiter) {
    # auxiliary quantities
    Sigma_wk <- Sigma %*% wk[1:N] # slicing is slow
    rk <- wk[1:N] * Sigma_wk
    Ak <- A(wk, Sigma, b, Sigma_w = Sigma_wk)
    g_wk <- g(wk, Sigma, b, r = rk)
    Qk <- 2 * crossprod(Ak) + tauI
    qk <- 2 * t(Ak) %*% g_wk - Qk %*% wk
    # build and solve problem (39) as in Feng & Palomar TSP2015
    w_hat <- quadprog::solve.QP(Qk, -qk, Amat = Amat,
                                bvec = bvec, meq = meq)$solution
    w_next <- wk + gamma * (w_hat - wk)
    # save objective function values and elapsed time
    time_seq <- c(time_seq, proc.time()[3] - start_time)
    fun_next <- R(w_next, Sigma, b)
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

  if (!has_theta)
    return(list(w = w_next,
           risk_contribution = as.vector(w_next * (Sigma %*% w_next)),
           obj_fun = fun_seq,
           elapsed_time = time_seq,
           convergence = sum(!(k == maxiter))))
  else
    return(list(w = w_next[1:N],
           risk_contribution = as.vector(w_next[1:N] * (Sigma %*% w_next[1:N])),
           obj_fun = fun_seq,
           elapsed_time = time_seq,
           convergence = sum(!(k == maxiter)),
           theta = w_next[N+1]))
}


#' Implements the risk parity portfolio using a general constrained
#' solver from the alabama package
#'
#' @export
riskParityPortfolioGenSolver <- function(Sigma, b = rep(1/nrow(Sigma), nrow(Sigma)),
                                         budget = TRUE, shortselling = FALSE,
                                         formulation = "rc-over-var vs b",
                                         method = "slsqp", use_gradient = TRUE,
                                         w0 = NA, theta0 = NA, maxiter = 500, ftol = 1e-9,
                                         wtol = 1e-6) {
  N <- nrow(Sigma)

  if (anyNA(w0))
    w0 <- riskParityPortfolioDiagSigma(Sigma, b)$w

  has_theta <- grepl("theta", formulation)
  if (has_theta) {
    if (is.na(theta0)) {
      r0 <- w0 * (Sigma %*% w0)
      theta0 <- mean(r0 / b)
    }
    w0 <- as.vector(c(w0, theta0))
  }

  if (budget) {
    if (has_theta) {
      budget <- function(w, ...)
        return(sum(w[1:N]) - 1)  # slicing is slow
      budget.jac <- function(w, ...)
        return(cbind(matrix(1, 1, N), 0))
    } else {
      budget <- function(w, ...)
        return(sum(w) - 1)
      budget.jac <- function(w, ...)
        return(matrix(1, 1, N))
    }
  } else {
    budget <- NULL
    budget.jac <- NULL
  }

  if (!shortselling) {
    if (has_theta) {
      shortselling <- function(w, ...)
        return(w[1:N])
      shortselling.jac <- function(w, ...)
        return(cbind(diag(N), rep(0, N)))
    } else {
      shortselling <- function(w, ...)
        return(w)
      shortselling.jac <- function(w, ...)
        return(diag(N))
    }
  } else {
    shortselling <- NULL
    shortselling.jac <- NULL
  }

  switch(formulation,
         "rc-double-index" = {
           R <- R_rc_double_index
           R_grad <- R_grad_rc_double_index
         },
         "rc-over-b-double-index" = {
           R <- R_rc_over_b_double_index
           R_grad <- R_grad_rc_over_b_double_index
         },
         "rc-over-var vs b" = {
           R <- R_rc_over_var_vs_b
           R_grad <- R_grad_rc_over_var_vs_b
         },
         "rc-over-var" = {
           R <- R_rc_over_var
           R_grad <- R_grad_rc_over_var
         },
         "rc-over-sd vs b-times-sd" = {
           R <- R_rc_over_sd_vs_b_times_sd
           R_grad <- R_grad_rc_over_sd_vs_b_times_sd
         },
         "rc vs b-times-var"  = {
           R <- R_rc_vs_b_times_var
           R_grad <- R_rc_vs_b_times_var
         },
         "rc vs theta" = {
           R <- R_rc_vs_theta
           R_grad <- R_grad_rc_vs_theta
         },
         "rc-over-b vs theta" = {
           R <- R_rc_over_b_vs_theta
           R_grad <- R_grad_rc_over_b_vs_theta
         },
         stop("formulation ", formulation, " is not included.")
  )
  if (!use_gradient) R_grad <- NULL

  fun_seq <- c(R(w0, Sigma, b))
  time_seq <- c(0)
  if (method == "alabama") {
    start_time <- proc.time()[3]
    res <- alabama::constrOptim.nl(w0, R, R_grad,
                                   hin = shortselling, hin.jac = shortselling.jac,
                                   heq = budget, heq.jac = budget.jac,
                                   Sigma = Sigma, b = b,
                                   control.outer = list(trace = FALSE, itmax = maxiter))
    end_time <- proc.time()[3]
  } else if (method == "slsqp") {
    start_time <- proc.time()[3]
    res <- nloptr::slsqp(w0, R, R_grad,
                         hin = shortselling, hinjac = shortselling.jac,
                         heq = budget, heqjac = budget.jac,
                         Sigma = Sigma, b = b, control = list(xtol_rel = wtol, ftol_rel = ftol))
    end_time <- proc.time()[3]
  }
  # save objective value and elapsed time
  time_seq <- c(time_seq, end_time - start_time)
  fun_seq <- c(fun_seq, res$value)
  w <- res$par

  if (!has_theta)
    return(list(w = w,
                risk_contribution = as.vector(w * (Sigma %*% w)),
                obj_fun = fun_seq,
                elapsed_time = time_seq,
                convergence = res$convergence))
  else
    return(list(w = w[1:N],
                risk_contribution = as.vector(w[1:N] * (Sigma %*% w[1:N])),
                obj_fun = fun_seq,
                elapsed_time = time_seq,
                convergence = res$convergence,
                theta = w[N+1]))
}
