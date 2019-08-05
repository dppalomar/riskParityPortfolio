riskParityPortfolioDiagSigma <- function(Sigma, b = rep(1/nrow(Sigma), nrow(Sigma))) {
  w <- sqrt(b) / sqrt(diag(Sigma))
  w <- w / sum(w)
  return (list(w = w, risk_contribution = as.vector(w * (Sigma %*% w))))
}


riskParityPortfolioSCA <- function(Sigma, b = rep(1/nrow(Sigma), nrow(Sigma)),
                                   mu = NULL, lmd_mu = 0, lmd_var = 0,
                                   w_lb = rep(0, nrow(Sigma)), w_ub = rep(1, nrow(Sigma)),
                                   Cmat = NULL, cvec = NULL, Dmat = NULL, dvec = NULL,
                                   formulation = c("rc-over-b-double-index",
                                                   "rc-double-index",
                                                   "rc-over-var vs b",
                                                   "rc-over-var",
                                                   "rc-over-sd vs b-times-sd",
                                                   "rc vs b-times-var",
                                                   "rc vs theta",
                                                   "rc-over-b vs theta"),
                                   w0 = NULL, theta0 = NULL, gamma = .9, zeta = 1e-7,
                                   tau = NULL, maxiter = 500, ftol = 1e-6, wtol = 1e-6,
                                   use_qp_solver = FALSE) {
  N <- nrow(Sigma)
  formulation <- match.arg(formulation)
  if (is.null(w0)) w0 <- riskParityPortfolioDiagSigma(Sigma, b)$w
  # define boolean cases for easier reading of code
  has_eq_constraints <- !(is.null(Cmat) || is.null(cvec))
  has_ineq_constraints <- !(is.null(Dmat) || is.null(dvec))
  #has_eq_and_ineq_constraints <- !(is.null(Dmat) || is.null(dvec) || is.null(Cmat) || is.null(cvec)) 
  #TODO{Vinicius}: please make sure the next line is correct
  has_eq_and_ineq_constraints <- has_eq_constraints & has_ineq_constraints
  has_mu <- !is.null(mu)
  has_var <- lmd_var > 0  
  has_theta <- grepl("theta", formulation)
  
  # computation of feasible initial point
  if (has_eq_and_ineq_constraints) {
    w0 <- project_onto_eq_and_ineq_constraint_set(w0, Cmat, cvec, Dmat, dvec)
    xi <- xi_prev <- rep(0, length(cvec))
    chi <- chi_prev <- rep(0, length(dvec))
  }
  #TODO{Vinicius}: why do we need a separate projection? We should have a single one. Let's skype
  else if (has_equality_constraints) {
    w0 <- project_onto_equality_constraint_set(w0, Cmat, cvec)
  } else {
    w0 <- projectBudgetLineAndBox(w0, w_lb, w_ub)
  }
  #TODO{Vinicius}: Hehe. The projection above is wrong because when you project on the linear constraints you forgot the other constraints. Let's skype
  if (has_theta) {
    if (is.null(theta0)) {
      r0 <- w0 * (Sigma %*% w0)
      theta0 <- switch(formulation,
                       "rc vs theta" = mean(r0),
                       "rc-over-b vs theta" = mean(r0 / b),
                       stop("Oops... Formulation name contains theta but it is not supported."))
      }
    w0 <- as.vector(c(w0, theta0))
  }
  
  # parameters tau for SCA method
  if (is.null(tau))
    tau <- .05 * sum(diag(Sigma)) / (2*N)
  tauI <- diag(rep(tau, length(w0)))
  
  # packing linear constrains
  if (has_theta) {  
    #TODO{Vinicius}: Hey, here we can remove the last column of Amat and last element of bvec!!
    Amat <- cbind(c(rep(1, N), 0), diag(rep(1, N+1)), -diag(c(rep(1, N), 0)))
    bvec <- c(1, c(w_lb, 0), c(-w_ub, 0))
  } else {
    Amat <- cbind(rep(1, N), diag(N), -diag(N))
    bvec <- c(1, w_lb, -w_ub)
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
  if (has_mu) {
    if (lmd_mu == 0) stop("Ooppss... You specified mu but chose lmd_mu == 0...")
    if (has_theta) fun_k <- fun_k - lmd_mu * t(mu) %*% wk[1:N]
    else fun_k <- fun_k - lmd_mu * t(mu) %*% wk
  }
  fun_seq <- c(fun_k)
  time_seq <- c(0)

  # SCA outer loop
  start_time <- proc.time()[3]
  for (k in 1:maxiter) {
    # auxiliary quantities
    if (has_theta) {
      Sigma_wk <- Sigma %*% wk[1:N]
      rk <- wk[1:N] * Sigma_wk
    } else {
      Sigma_wk <- Sigma %*% wk
      rk <- wk * Sigma_wk
    }
    Ak <- A(wk, Sigma, b, Sigma_w = Sigma_wk)
    g_wk <- g(wk, Sigma, b, r = rk)
    Qk <- 2 * crossprod(Ak) + tauI
    if (has_var)
      if (has_theta) Qk <- Qk + lmd_var * cbind(rbind(Sigma, rep(0, N)), rep(0, N+1))
      else Qk <- Qk + lmd_var * Sigma
    qk <- 2 * t(Ak) %*% g_wk - Qk %*% wk
    if (has_mu)
      if (has_theta) qk <- qk - lmd_mu * c(mu, 0)
      else qk <- qk - lmd_mu * mu
    # build and solve problem (39) as in Feng & Palomar TSP2015
    if (has_eq_and_ineq_constraints || has_equality_constraints) {
      if (!use_qp_solver){
        if (has_eq_and_ineq_constraints) {
          params <- rpp_eq_and_ineq_constraints_iteration(Cmat, cvec, Dmat, dvec, Qk,
                                                          qk, wk, chi, chi_prev, xi,
                                                          xi_prev)
          chi_prev <- params[[1]]
          chi <- params[[2]]
          xi_prev <- params[[3]]
          xi <- params[[4]]
          w_hat <- params[[5]]
        } else if (has_equality_constraints) {
          w_hat <- rpp_equality_constraints_iteration(Cmat, cvec, Qk, qk)
        }
      } else {
        if (has_equality_constraints) {
          meq <- nrow(Cmat)
          w_hat <- quadprog::solve.QP(Qk, -qk, Amat = t(Cmat), bvec = cvec,
                                      meq = meq)$solution
        } else if (has_eq_and_ineq_constraints) {
          meq <- nrow(Cmat)
          Amat <- rbind(Cmat, -Dmat)
          bvec <- c(cvec, -dvec)
          w_hat <- quadprog::solve.QP(Qk, -qk, Amat = t(Amat), bvec = bvec,
                                      meq = meq)$solution
        }
      }
    } else {
      w_hat <- quadprog::solve.QP(Qk, -qk, Amat = Amat, bvec = bvec,
                                  meq = 1)$solution
    }
    w_next <- wk + gamma * (w_hat - wk)
    # save objective function value and elapsed time
    time_seq <- c(time_seq, proc.time()[3] - start_time)
    fun_next <- R(w_next, Sigma, b)
    if (has_mu)
      if (has_theta) fun_next <- fun_next - lmd_mu * t(mu) %*% w_next[1:N]
      else fun_next <- fun_next - lmd_mu * t(mu) %*% w_next
    if (has_var)
      if (has_theta) fun_next <- fun_next + lmd_var * (t(w_next[1:N]) %*% Sigma %*% w_next[1:N])
      else fun_next <- fun_next + lmd_var * (t(w_next) %*% Sigma %*% w_next)
    fun_seq <- c(fun_seq, fun_next)
    # check convergence on parameters and objective function
    has_w_converged <- all((abs(w_next - wk) <= .5 * wtol * (abs(wk) + abs(w_next))))
    has_fun_converged <- abs(fun_next - fun_k) <= .5 * ftol * (abs(fun_k) + abs(fun_next))
    if (k > 1 && (has_w_converged && has_fun_converged)) break
    # update variables
    wk <- w_next
    fun_k <- fun_next
    gamma <- gamma * (1 - zeta * gamma)
  }

  portfolio_results <- list()
  portfolio_results$risk_parity <- fun_seq[length(fun_seq)]
  if (!has_theta) {
    portfolio_results$w <- w_next
    portfolio_results$risk_contribution <- as.vector(w_next * (Sigma %*% w_next))
  } else {
    portfolio_results$w <- w_next[1:N]
    portfolio_results$theta <- w_next[N+1]
    portfolio_results$risk_contribution <- as.vector(w_next[1:N] * (Sigma %*% w_next[1:N]))
  }
  if (has_mu) {
    portfolio_results$mean_return <- t(mu) %*% portfolio_results$w
    portfolio_results$risk_parity <- portfolio_results$risk_parity + lmd_mu * portfolio_results$mean_return
  }
  if (has_var) {
    portfolio_results$variance <- t(portfolio_results$w) %*% Sigma %*% portfolio_results$w
    portfolio_results$risk_parity <- portfolio_results$risk_parity - lmd_var * portfolio_results$variance
  }
  portfolio_results$obj_fun <- fun_seq
  portfolio_results$elapsed_time <- time_seq
  portfolio_results$convergence <- !(k == maxiter)
  return(portfolio_results)
}


riskParityPortfolioGenSolver <- function(Sigma, b = rep(1/nrow(Sigma), nrow(Sigma)),
                                         mu = NULL, lmd_mu = 0,
                                         formulation = c("rc-over-b-double-index",
                                                         "rc-double-index",
                                                         "rc-over-var vs b",
                                                         "rc-over-var",
                                                         "rc-over-sd vs b-times-sd",
                                                         "rc vs b-times-var",
                                                         "rc vs theta",
                                                         "rc-over-b vs theta"),
                                         method = c("slsqp", "alabama"),
                                         use_gradient = TRUE, w0 = NULL, theta0 = NULL,
                                         maxiter = 500, ftol = 1e-6, wtol = 1e-6) {
  N <- nrow(Sigma)
  if (is.null(w0)) w0 <- riskParityPortfolioDiagSigma(Sigma, b)$w
  formulation <- match.arg(formulation)
  # set initial value for theta
  has_theta <- grepl("theta", formulation)
  if (has_theta) {
    if (is.null(theta0)) {
      r0 <- w0 * (Sigma %*% w0)
      if (formulation == "rc vs theta")
        theta0 <- mean(r0)
      else if (formulation == "rc-over-b vs theta")
        theta0 <- mean(r0 / b)
    }
    w0 <- as.vector(c(w0, theta0))
  }
  if (has_theta) {
    budget <- function(w, ...) {
      N <- length(w) - 1
      return(sum(w[1:N]) - 1)
    }
    budget.jac <- function(w, ...) {
      N <- length(w) - 1
      return(cbind(matrix(1, 1, N), 0))
    }
    shortselling <- function(w, ...) {
      N <- length(w) - 1
      return(w[1:N])
    }
    shortselling.jac <- function(w, ...) {
      N <- length(w) - 1
      return(cbind(diag(N), rep(0, N)))
    }
  } else {
    budget <- function(w, ...)
      return(sum(w) - 1)
    budget.jac <- function(w, ...) {
      N <- length(w)
      return(matrix(1, 1, N))
    }
    shortselling <- function(w, ...)
      return(w)
    shortselling.jac <- function(w, ...)
      return(diag(length(w)))
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
           R_grad <- R_grad_rc_vs_b_times_var
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
  has_mu <- !is.null(mu)
  if (has_mu) {
    wrap_R <- function(R, lmd_mu, mu, has_theta, N) {
      if (has_theta) {
        func <- function(...) {
          kwargs <- list(...)
          w_ <- kwargs[[1]]
          return(R(...) - lmd_mu * t(mu) %*% w_[1:N])
        }
      } else {
        func <- function(...) {
          kwargs <- list(...)
          return(R(...) - lmd_mu * t(mu) %*% kwargs[[1]])
        }
      }
      return(func)
    }
    wrap_R_grad <- function(R_grad, lmd_mu, mu) {
      grad <- function(...) {
        return(R_grad(...) - lmd_mu * mu)
      }
      return(grad)
    }
    R_ <- wrap_R(R, lmd_mu, mu, has_theta, N)
    R_grad_ <- wrap_R_grad(R_grad, lmd_mu, mu)
  } else {
    R_ <- R
    R_grad_ <- R_grad
  }
  fun_seq <- c(R(w0, Sigma, b))
  time_seq <- c(0)
  switch(match.arg(method),
         "alabama" = {
           start_time <- proc.time()[3]
           res <- alabama::constrOptim.nl(w0, R_, R_grad_,
                                          hin = shortselling,
                                          hin.jac = shortselling.jac,
                                          heq = budget,
                                          heq.jac = budget.jac,
                                          Sigma = Sigma, b = b,
                                          control.outer = list(trace = FALSE,
                                                               itmax = maxiter))
           end_time <- proc.time()[3]
         },
         "slsqp" = {
           start_time <- proc.time()[3]
           res <- nloptr::slsqp(w0, R_, R_grad_,
                                hin = shortselling, hinjac = shortselling.jac,
                                heq = budget, heqjac = budget.jac,
                                Sigma = Sigma, b = b,
                                control = list(xtol_rel = wtol, ftol_rel = ftol))
           end_time <- proc.time()[3]
         },
         stop("method ", method, "is not included.")
  )
  # save objective value and elapsed time
  time_seq <- c(time_seq, end_time - start_time)
  fun_seq <- c(fun_seq, res$value)
  w <- res$par

  portfolio_results <- list()
  portfolio_results$risk_parity <- fun_seq[length(fun_seq)]
  if (!has_theta) {
    portfolio_results$w <- w
    portfolio_results$risk_contribution <- as.vector(w * (Sigma %*% w))
  } else {
    portfolio_results$w <- w[1:N]
    portfolio_results$theta <- w[N+1]
    portfolio_results$risk_contribution <- as.vector(w[1:N] * (Sigma %*% w[1:N]))
  }
  if (has_mu) {
    portfolio_results$mean_return <- t(mu) %*% portfolio_results$w
    portfolio_results$risk_parity <- portfolio_results$risk_parity + lmd_mu * portfolio_results$mean_return
  }
  portfolio_results$obj_fun <- fun_seq
  portfolio_results$elapsed_time <- time_seq
  portfolio_results$convergence <- res$convergence
  return(portfolio_results)
}


riskParityPortfolioNewton <- function(Sigma, b = rep(1/nrow(Sigma), nrow(Sigma)),
                                      maxiter = 50, ftol = 1e-8) {
  w <- risk_parity_portfolio_nn(Sigma, b, ftol, maxiter)
  return(list(w = w, risk_contribution = c(w * (Sigma %*% w)),
              obj_fun = obj_function_spinu(Sigma, w, b)))
}


riskParityPortfolioCyclicalRoncalli <- function(Sigma, b = rep(1/nrow(Sigma), nrow(Sigma)),
                                        maxiter = 50, ftol = 1e-8) {
  w <- risk_parity_portfolio_ccd_roncalli(Sigma, b, ftol, maxiter)
  return(list(w = w, risk_contribution = c(w * (Sigma %*% w)),
              obj_fun = obj_function_roncalli(Sigma, w, b)))
}


riskParityPortfolioCyclicalSpinu <- function(Sigma, b = rep(1/nrow(Sigma), nrow(Sigma)),
                                             maxiter = 50, ftol = 1e-8) {
  w <- risk_parity_portfolio_ccd_spinu(Sigma, b, ftol, maxiter)
  return(list(w = w, risk_contribution = c(w * (Sigma %*% w))))
}

# minimize ||w - w0||^2
# s.t.     sum(w) = 1
#          w_lb <= w <= w_ub
projectBudgetLineAndBox <- function(w0, w_lb, w_ub) {
  if (sum(w_lb) > 1 || sum(w_ub) < 1) stop("Problem infeasible: relax the bounds!")

  obj_fun <- function(mu, w0) {
    sum(pmax(pmin(w0 - mu, w_ub), w_lb)) - 1
  }
  mu_ub <- max(w0 - w_lb)
  mu_lb <- min(w0 - w_ub)
  mu <- stats::uniroot(obj_fun, interval = c(mu_lb, mu_ub), w0, tol = 1e-8)$root
  w <- pmax(pmin(w0 - mu, w_ub), w_lb)
  return(w)
}


# minimize ||w - w0||^2
# s.t.     Cw = c
#          Dw <= d
project_onto_eq_and_ineq_constraint_set <- function(w0, Cmat, cvec, Dmat, dvec) {
  meq <- nrow(Cmat)
  I <- diag(length(w0))
  Amat <- -t(rbind(Cmat, Dmat))
  b0 <- -c(cvec, dvec)
  return(quadprog::solve.QP(Dmat = I, dvec = w0, Amat = Amat, bvec = b0, meq = meq)$solution)
}


#' @title Design of risk parity portfolios
#'
#' @description This function designs risk parity portfolios to equalize/distribute
#' the risk contributions of the different assets, which is missing if we simply
#' consider the overall volatility of the portfolio as in the mean-variance
#' Markowitz portfolio. In addition to the vanilla formulation, where the risk
#' contributions are perfectly equalized subject to no shortselling and budget
#' constraints, many other formulations are considered that allow for box
#' constraints, as well as the inclusion of additional objectives like the
#' expected return and overall variance. In short, this function solves the
#' following problem:
#'
#'       \code{minimize    R(w) - lmd_mu * t(w) \%*\% mu + lmd_var * t(w) \%*\% Sigma \%*\% w}
#'
#'       \code{subject to  sum(w) = 1, w_lb <= w <= w_ub},
#'       \code{            Cmat \%*\% w = cvec, Dmat \%*\% w <= dvec},
#'
#' where \code{R(w)} denotes the risk concentration,
#' \code{t(w) \%*\% mu} is the expected return, \code{t(w) \%*\% Sigma \%*\% w} is the
#' overall variance, \code{lmd_mu} and \code{lmd_var} are the trade-off weights
#' for the expected return and the variance terms, respectively, \code{w_lb} and
#' \code{w_ub} are the lower and upper bound vector values for the portfolio vector \code{w},
#' \code{Cmat \%*\% w = cvec} denotes arbitrary linear equality constrains, and
#' \code{Dmat \%*\% w = dvec} denotes arbitrary linear inequality constrains.
#'
#' @details By default, the problem considered is the vanilla risk parity portfolio:
#' \code{w >= 0, sum(w) = 1}, with no expected return term, and no variance term. In this case,
#' the problem formulation is convex and the optimal solution is guaranteed to be achieved with
#' a perfect risk concentration, i.e., \code{R(w) = 0}. By default, we use the formulation by
#' Spinu (2013) (\code{method_init = "cyclical-spinu"}), but the user can also select the formulation
#' by Roncalli et al. (2013) (\code{method_init = "cyclical-roncalli"}).
#'
#' In case of additional box constraints, expected return term, or variance term,
#' then the problem is nonconvex and the global optimal solution cannot be
#' guaranteed, just a local optimal. We use the efficient sucessive
#' convex approximation (SCA) method proposed in Feng & Palomar (2015),
#' where the user can choose among many different risk concentration
#' terms (through the argument \code{formulation}), namely:
#' \itemize{
#' \item{\code{formulation = "rc-double-index":} }{\code{sum_{i,j} (r_i - r_j)^2}}
#' \item{\code{formulation = "rc-vs-theta":} }{\code{sum_{i} (r_i - theta)^2}}
#' \item{\code{formulation = "rc-over-var-vs-b":} }{\code{sum_{i} (r_i/r - b_i)^2}}
#' \item{\code{formulation = "rc-over-b double-index":} }{\code{sum_{i,j} (r_i/b_i - r_j/b_j)^2}}
#' \item{\code{formulation = "rc-vs-b-times-var":} }{\code{sum_{i} (r_i - b_i*r)^2}}
#' \item{\code{formulation = "rc-over-sd vs b-times-sd":} }{\code{sum_{i} (r_i/sqrt(r) - b_i*sqrt(r))^2}}
#' \item{\code{formulation = "rc-over-b vs theta":} }{\code{sum_{i} (r_i/b_i - theta)^2}}
#' \item{\code{formulation = "rc-over-var":} }{\code{sum_{i} (r_i/r)^2}}}
#' where \code{r_i = w_i*(Sigma\%*\%w)_i} is the risk contribution and
#' \code{r = t(w)\%*\%Sigma\%*\%w} is the overall risk (i.e., variance).
#'
#' For more details, please check the vignette.
#'
#' @param Sigma covariance or correlation matrix (this is the only mandatory argument)
#' @param b budget vector, i.e., the risk budgeting targets. The default is the
#'        uniform 1/N vector.
#' @param mu vector of expected returns (only needed if the expected return term
#'        is desired in the objective)
#' @param lmd_mu scalar weight to control the importance of the expected return term
#' @param lmd_var scalar weight to control the importance of the variance term
#'        (only currently available for the SCA method)
#' @param w_lb lower bound (either a vector or a scalar) on the value of each
#'        portfolio weight
#' @param w_ub upper bound (either a vector or a scalar) on the value of each
#'        portfolio weight
#' @param Cmat TBD
#' @param cvec TBD
#' @param Dmat TBD
#' @param dvec TBD
#' @param method_init method to compute the vanilla solution. In case of
#'        additional constraints or objective terms, this solution is used as
#'        the initial point for the subsequent method. The default is
#'        \code{"cyclical-spinu"}. See details below.
#' @param method method to solve the non-vanilla formulation. The default is \code{"sca"}.
#'        See details below. (DEPRECATED)
#' @param formulation string indicating the risk concentration formulation to be used.
#'        It must be one of: "diag", "rc-double-index",
#'        "rc-over-b-double-index", "rc-over-var vs b",
#'        "rc-over-var", "rc-over-sd vs b-times-sd",
#'        "rc vs b-times-var", "rc vs theta", or
#'        "rc-over-b vs theta". The default is "rc-over-b-double-index".
#'        If \code{formulation} is not provided and no additional terms or
#'        constraints are set, such as expected return or shortselling, then the
#'        vanilla risk parity portfolio will be returned. If formulation is
#'        "diag" then the analytical solution of the risk parity optimization for
#'        for a diagonal covariance matrix will be returned. See details below.
#' @param w0 initial value for the portfolio weights. Default is a convex
#'        combination of the risk parity portfolio, the (uncorrelated) minimum variance
#'        portfolio, and the maximum return portfolio.
#' @param theta0 initial value for theta (in case formulation uses theta). If not provided,
#'        the optimum solution for a fixed vector of portfolio weights will be used.
#' @param gamma learning rate for the SCA method
#' @param zeta factor used to decrease the learning rate at each iteration for the SCA method
#' @param tau regularization factor
#' @param maxiter maximum number of iterations for the SCA loop
#' @param ftol convergence tolerance on the objective function
#' @param wtol convergence tolerance on the values of the portfolio weights
#' @param use_gradient this parameter is meaningful only if method is either
#'        \code{"alabama"} or \code{"slsqp"}. If \code{TRUE} (default value), analytical gradients of the
#'        objective function will be used (strongly recommended to achieve faster results).
#' @return A list containing possibly the following elements:
#' \item{\code{w}}{optimal portfolio vector}
#' \item{\code{risk_contribution}}{the risk contribution of every asset}
#' \item{\code{theta}}{the optimal value for theta (in case that it is part of
#'                     the chosen formulation)}
#' \item{\code{obj_fun}}{the sequence of values of the objective function at
#'                       each iteration}
#' \item{\code{risk_parity}}{the risk concentration term of the portfolio \code{R(w)}}
#' \item{\code{mean_return}}{the expected return term of the portoflio \code{t(w)\%*\%mu},
#'                           if the term is included in the optimization}
#' \item{\code{variance}}{the variance term of the portfolio \code{t(w)\%*\%Sigma\%*\%w},
#'                        if the term is included in the optimization}
#' \item{\code{elapsed_time}}{elapsed time recorded at every iteration}
#' \item{\code{convergence}}{boolean flag to indicate whether or not the optimization converged}
#'
#' @examples
#' library(riskParityPortfolio)
#'
#' # create covariance matrix
#' N <- 5
#' V <- matrix(rnorm(N^2), ncol = N)
#' Sigma <- cov(V)
#'
#' # risk parity portfolio
#' res <- riskParityPortfolio(Sigma)
#' names(res)
#' #> [1] "w"                 "risk_contribution"
#'
#' res$w
#' #> [1] 0.04142886 0.38873465 0.34916787 0.09124019 0.12942842
#'
#' res$risk_contribution
#' #> [1] 0.007361995 0.007361995 0.007361995 0.007361995 0.007361995
#'
#' c(res$w * (Sigma %*% res$w))
#' #> [1] 0.007361995 0.007361995 0.007361995 0.007361995 0.007361995
#'
#' # risk budggeting portfolio
#' res <- riskParityPortfolio(Sigma, b = c(0.4, 0.4, 0.1, 0.05, 0.05))
#' res$risk_contribution/sum(res$risk_contribution)
#' #> [1] 0.40 0.40 0.10 0.05 0.05
#'
#' @references
#' Y. Feng, and D. P. Palomar (2015). SCRIP: Successive Convex Optimization Methods
#' for Risk Parity Portfolio Design. \emph{IEEE Trans. on Signal Processing},
#' vol. 63, no. 19, pp. 5285-5300. <https://doi.org/10.1109/TSP.2015.2452219>
#'
#' F. Spinu (2013). An Algorithm for Computing Risk Parity Weights.
#' <https://dx.doi.org/10.2139/ssrn.2297383>
#'
#' T. Griveau-Billion, J. Richard, and T. Roncalli (2013). A fast algorithm for computing High-dimensional
#' risk parity portfolios. <https://arxiv.org/pdf/1311.4057.pdf>
#'
#' @author Ze Vinicius and Daniel P. Palomar
#'
#' @export
riskParityPortfolio <- function(Sigma, b = NULL, mu = NULL,
                                lmd_mu = 0, lmd_var = 0,
                                w_lb = 0, w_ub = 1,
                                Cmat = NULL, cvec = NULL, Dmat = NULL, dvec = NULL,
                                method_init = c("cyclical-spinu", "cyclical-roncalli", "newton"),
                                method = c("sca", "alabama", "slsqp"),
                                formulation = NULL, w0 = NULL, theta0 = NULL,
                                gamma = .9, zeta = 1e-7, tau = NULL,
                                maxiter = 500, ftol = 1e-8, wtol = 1e-6,
                                use_gradient = TRUE, use_qp_solver = FALSE) {
  # default values
  stocks_names <- colnames(Sigma)
  N <- nrow(Sigma)
  if (is.null(b)) b <- rep(1/N, N)
  if (length(w_ub) == 1) w_ub <- rep(w_ub, N)
  if (length(w_lb) == 1) w_lb <- rep(w_lb, N)
  # define boolean cases for easier reading of code
  has_mu <- !is.null(mu)
  has_theta <- !is.null(theta0)
  has_var <- lmd_var > 0
  has_formulation <- !is.null(formulation)
  has_fancy_box <- any(w_lb != 0) || any(w_ub != 1)
  has_equality_constraints <- !(is.null(Cmat) || is.null(cvec))
  has_inequality_constraints <- !(is.null(Dmat) || is.null(dvec))
  has_initial_point <- !is.null(w0)
  is_vanilla_formulation <- !(has_mu || has_theta || has_var || has_fancy_box ||
                              has_equality_constraints || has_inequality_constraints)

  # check wrong parameters
  if (sum(w_lb) > 1) stop("Problem infeasible: relax the lower bounds.")
  if (sum(w_ub) < 1) stop("Problem infeasible: relax the upper bounds.")
  if (length(b) != N) stop("Shape mismatch: b has to have nrow(Sigma) number of elements.")
  if (has_mu && (length(mu) != N))
    stop("Shape mismatch: mu has to have nrow(Sigma) number of elements")
  if (has_initial_point && (length(w0) != N))
    stop("Shape mismatch: w0 has to have nrow(Sigma) number of elements")
  # if diag, then call diagonal solver
  if (has_formulation && formulation == "diag") {
    if (!is_vanilla_formulation)
      stop("Additional constraints (box-constraints or other linear constrains) or",
           " terms (expected return or variance) are not supported by the 'diag' formulation.")
    if (has_initial_point)
      warning("The problem is a naive (diagonal) risk parity portfolio, but an initial",
              " point has been provided: The initial point is being ignored.")
    return(riskParityPortfolioDiagSigma(Sigma, b))
  }
  # if vanilla, then call the vanilla solver
  if (is_vanilla_formulation && has_formulation) {
    warning("The problem is a vanilla risk parity portofolio, but a nonconvex",
            " formulation has been chosen. Consider not specifying the formulation",
            " argument in order to use the convex formulation and get a guaranteed",
            " global solution.")
    is_vanilla_formulation <- FALSE  # so the nonconvex formulation will be used
  }
  if (is_vanilla_formulation) {
    if (has_initial_point)
      warning("The problem is a vanilla risk parity portfolio, but an initial",
              " point has been provided: The initial point is being ignored.")
    portfolio <- switch(match.arg(method_init),
                        "newton" = riskParityPortfolioNewton(Sigma, b, maxiter, ftol),
                        "cyclical-spinu" = riskParityPortfolioCyclicalSpinu(Sigma, b, maxiter, ftol),
                        "cyclical-roncalli" = riskParityPortfolioCyclicalRoncalli(Sigma, b, maxiter, ftol),
                        stop("method_init ", method_init, " is not supported."))
  } else {  # nonconvex solver
    if (!has_initial_point) {
      w_rc <- switch(match.arg(method_init),
                     "newton" = riskParityPortfolioNewton(Sigma, b, maxiter, ftol)$w,
                     "cyclical-spinu" = riskParityPortfolioCyclicalSpinu(Sigma, b, maxiter, ftol)$w,
                     "cyclical-roncalli" = riskParityPortfolioCyclicalRoncalli(Sigma, b, maxiter, ftol)$w,
                   stop("method_init ", method_init, " is not supported."))
      # create fancy initial point for the case of additional objectives
      w_gmvp <- 1 / diag(Sigma)
      w_gmvp <- w_gmvp / sum(w_gmvp)
      if (has_mu) {
        w_mu <- as.numeric(max(mu) == mu)
        w_mu <- w_mu / sum(w_mu)
      } else w_mu <- 0
      theta_rc <- 1 / (1 + lmd_var + lmd_mu*sum(has_mu))
      theta_mu <- lmd_mu*sum(has_mu) / (1 + lmd_var + lmd_mu*sum(has_mu))
      theta_var <- lmd_var / (1 + lmd_var + lmd_mu*sum(has_mu))
      w0 <- theta_rc*w_rc + theta_mu*w_mu + theta_var*w_gmvp
    }
    # TODO: Here we should do different projections depending on whether there are or not linear constraints
    # Also, should we incorporate the other constraints in the general matrices? (sum(w)==1 and box constraints)
    # make w0 feasible (excluding the general linear constraints)
    if (sum(w0) != 1 || any(w0 < w_lb) || any(w0 > w_ub)) {
      if (has_initial_point) warning("Initial point is infeasible. Projecting it onto the feasible set.")
      w0 <- projectBudgetLineAndBox(w0, w_lb, w_ub)
    }
    # solve nonconvex formulation
    switch(match.arg(method),
           "sca" = portfolio <- riskParityPortfolioSCA(Sigma = Sigma, b = b, mu = mu, lmd_mu = lmd_mu,
                                                       lmd_var = lmd_var, w_lb = w_lb, w_ub = w_ub,
                                                       Cmat = Cmat, cvec = cvec, Dmat = Dmat, dvec = dvec,
                                                       formulation = formulation, w0 = w0, theta0 = theta0,
                                                       gamma = gamma, zeta = zeta, tau = tau,
                                                       maxiter = maxiter, ftol = ftol, wtol = wtol,
                                                       use_qp_solver = use_qp_solver),
           "slsqp" = ,
           "alabama" = {
             warning("Methods 'slsqp' and 'alabama' are deprecated and will be removed in the next release.
                      We strongly recommend the more robust and scalable method = 'sca'.")
             # reminder: remove argument use_gradient after decrecating...
             if (has_fancy_box) stop("Box constraints are not supported for method ", method)
             if (has_var) stop("Variance term is not supported for method ", method)
             if (has_equality_constraints || has_inequality_constraints) 
               stop("General linear constraints not supported for method ", method)
             portfolio <- riskParityPortfolioGenSolver(Sigma = Sigma, b = b, mu = mu, lmd_mu = lmd_mu,
                                                       formulation = formulation, method = method,
                                                       use_gradient = use_gradient, w0 = w0, theta0 = theta0,
                                                       maxiter = maxiter, ftol = ftol, wtol = wtol)
           },
           stop("method ", method, " is not included.")
    )
  }
  names(portfolio$w) <- names(portfolio$risk_contribution) <- stocks_names
  return(portfolio)
}
