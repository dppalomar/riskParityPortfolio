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
                                         maxiter = 100, ftol = 1e-8, wtol = .5e-6) {
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
  portfolio_results$risk_concentration <- fun_seq[length(fun_seq)]
  if (!has_theta) {
    portfolio_results$w <- w
    w_Sigmaw <- as.vector(w * (Sigma %*% w))
    portfolio_results$relative_risk_contribution <- w_Sigmaw / sum(w_Sigmaw)
  } else {
    portfolio_results$w <- w[1:N]
    portfolio_results$theta <- w[N+1]
    w_Sigmaw <- as.vector(w[1:N] * (Sigma %*% w[1:N]))
    portfolio_results$relative_risk_contribution <- w_Sigmaw / sum(w_Sigmaw)
  }
  if (has_mu) {
    portfolio_results$mean_return <- t(mu) %*% portfolio_results$w
    portfolio_results$risk_concentration <- portfolio_results$risk_concentration + lmd_mu * portfolio_results$mean_return
  }
  portfolio_results$obj_fun <- fun_seq
  portfolio_results$elapsed_time <- time_seq
  portfolio_results$convergence <- res$convergence
  return(portfolio_results)
}
