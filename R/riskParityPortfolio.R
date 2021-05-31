riskParityPortfolioDiagSigma <- function(Sigma, b = rep(1/nrow(Sigma), nrow(Sigma))) {
  w <- sqrt(b) / sqrt(diag(Sigma))
  w <- w / sum(w)
  w_Sigmaw <- as.vector(w * (Sigma %*% w))
  return (list(w = w, relative_risk_contribution = w_Sigmaw / sum(w_Sigmaw)))
}


isFeasiblePortfolio <- function(w, Cmat, cvec, Dmat, dvec, tol = 1e-6) {
  equality_feasibility <- all(abs(Cmat %*% w - cvec) < tol)
  if (is.null(Dmat))
    return (equality_feasibility)
  inequality_feasibility <- all(Dmat %*% w - dvec <= tol)
  return (equality_feasibility && inequality_feasibility)
}


riskParityPortfolioSCA <- function(Sigma, w0, b = rep(1/nrow(Sigma), nrow(Sigma)),
                                   mu = NULL, lmd_mu = 0, lmd_var = 0,
                                   w_lb = rep(0, nrow(Sigma)), w_ub = rep(1, nrow(Sigma)),
                                   Cmat = matrix(1, 1, nrow(Sigma)), cvec = c(1),
                                   Dmat = rbind(diag(nrow(Sigma)), -diag(nrow(Sigma))),
                                   dvec = c(w_ub, -w_lb),
                                   formulation = c("rc-over-var vs b",
                                                   "rc-over-b-double-index",
                                                   "rc-double-index",
                                                   "rc-over-var",
                                                   "rc-over-sd vs b-times-sd",
                                                   "rc vs b-times-var",
                                                   "rc vs theta",
                                                   "rc-over-b vs theta"),
                                   theta0 = NULL, gamma = .9, zeta = 1e-7,
                                   tau = NULL, maxiter = 1000, ftol = 1e-8, wtol = .5e-6,
                                   use_qp_solver = TRUE) {
  N <- nrow(Sigma)
  formulation <- match.arg(formulation)
  # define boolean cases for easier reading of code
  has_mu <- !is.null(mu)
  has_var <- lmd_var > 0
  has_theta <- grepl("theta", formulation)
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
    Cmat <- cbind(Cmat, 0)
    Dmat <- cbind(Dmat, 0)
  }
  # check the type of constraints
  has_only_equality_constraints <- all(w_lb == (-Inf), w_ub == Inf)
  # initiliaze some variables depending on the type of solver
  if (use_qp_solver) {
    meq <- nrow(Cmat)
    Amat <- t(rbind(Cmat, -Dmat))
    bvec <- c(cvec, -dvec)
  } else {
    dual_lmd_0 <- dual_lmd_minus_1 <- rep(0, length(cvec))
    dual_mu_0 <- dual_mu_minus_1 <- rep(0, length(dvec))
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
    if (lmd_mu == 0) warning("The mean return vector has been given, but lmd_mu = 0.")
    if (has_theta) fun_k <- fun_k - lmd_mu * t(mu) %*% wk[1:N]
    else fun_k <- fun_k - lmd_mu * t(mu) %*% wk
  }
  if (has_var) fun_k <- fun_k + lmd_var * t(wk) %*% Sigma %*% wk
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
    # solve problem (39) as in Feng & Palomar TSP2015
    if (use_qp_solver) {
      w_hat <- quadprog::solve.QP(Qk, -qk, Amat = Amat, bvec = bvec, meq = meq)$solution
    } else if (has_only_equality_constraints) {
      w_hat <- rpp_equality_constraints_iteration(Cmat, cvec, Qk, qk)
    } else {
      params <- rpp_eq_and_ineq_constraints_iteration(Cmat, cvec, Dmat, dvec, Qk, qk, wk,
                                                      dual_mu_0, dual_mu_minus_1, dual_lmd_0,
                                                      dual_lmd_minus_1, maxiter, wtol)
      dual_mu_minus_1 <- params[[1]]
      dual_mu_0 <- params[[2]]
      dual_lmd_minus_1 <- params[[3]]
      dual_lmd_0 <- params[[4]]
      w_hat <- params[[5]]
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
    # check convergence on parameters and objective function
    has_w_converged <- all((abs(w_next - wk) <= .5 * wtol * (abs(wk) + abs(w_next))))
    has_fun_converged <- (abs(fun_next - fun_k) <= .5 * ftol * (abs(fun_k) + abs(fun_next)))
    if (k > 1 && (has_w_converged || has_fun_converged)) break
    # update variables
    wk <- w_next
    fun_seq <- c(fun_seq, fun_next)
    fun_k <- fun_next
    gamma <- gamma * (1 - zeta * gamma)
  }

  portfolio_results <- list()
  portfolio_results$risk_concentration <- fun_seq[length(fun_seq)]
  if (!has_theta) {
    portfolio_results$w <- w_next
    w_Sigmaw <- as.vector(w_next * (Sigma %*% w_next))
    portfolio_results$relative_risk_contribution <- w_Sigmaw / sum(w_Sigmaw)
  } else {
    portfolio_results$w <- w_next[1:N]
    portfolio_results$theta <- w_next[N+1]
    w_Sigmaw <- as.vector(w_next[1:N] * (Sigma %*% w_next[1:N]))
    portfolio_results$relative_risk_contribution <- w_Sigmaw / sum(w_Sigmaw)
  }
  if (has_mu) {
    portfolio_results$mean_return <- t(mu) %*% portfolio_results$w
    portfolio_results$risk_concentration <- portfolio_results$risk_concentration + lmd_mu * portfolio_results$mean_return
  }
  if (has_var) {
    portfolio_results$variance <- t(portfolio_results$w) %*% Sigma %*% portfolio_results$w
    portfolio_results$risk_concentration <- portfolio_results$risk_concentration - lmd_var * portfolio_results$variance
  }
  portfolio_results$obj_fun <- fun_seq
  portfolio_results$elapsed_time <- time_seq
  portfolio_results$convergence <- !(k == maxiter)
  return(portfolio_results)
}


riskParityPortfolioNewton <- function(Sigma, b = rep(1/nrow(Sigma), nrow(Sigma)),
                                      maxiter = 50, ftol = 1e-8) {
  w <- risk_parity_portfolio_nn(Sigma, b, ftol, maxiter)
  w_Sigmaw <- c(w * (Sigma %*% w))
  return(list(w = w, relative_risk_contribution = w_Sigmaw / sum(w_Sigmaw),
              obj_fun = obj_function_spinu(Sigma, w, b)))
}


riskParityPortfolioCyclicalRoncalli <- function(Sigma, b = rep(1/nrow(Sigma), nrow(Sigma)),
                                        maxiter = 50, ftol = 1e-8) {
  w <- risk_parity_portfolio_ccd_roncalli(Sigma, b, ftol, maxiter)
  w_Sigmaw <- c(w * (Sigma %*% w))
  return(list(w = w, relative_risk_contribution = w_Sigmaw / sum(w_Sigmaw),
              obj_fun = obj_function_roncalli(Sigma, w, b)))
}


riskParityPortfolioCyclicalSpinu <- function(Sigma, b = rep(1/nrow(Sigma), nrow(Sigma)),
                                             maxiter = 50, ftol = 1e-8) {
  w <- risk_parity_portfolio_ccd_spinu(Sigma, b, ftol, maxiter)
  w_Sigmaw <- c(w * (Sigma %*% w))
  return(list(
              w = w,
              relative_risk_contribution = w_Sigmaw / sum(w_Sigmaw),
              obj_fun = obj_function_spinu(Sigma, w, b)
        )
  )
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
#' @param Sigma Covariance or correlation matrix (this is the only mandatory argument).
#' @param b Budget vector, i.e., the risk budgeting targets. The default is the
#'        uniform 1/N vector.
#' @param mu Vector of expected returns (only needed if the expected return term
#'        is desired in the objective).
#' @param lmd_mu Scalar weight to control the importance of the expected return term.
#' @param lmd_var Scalar weight to control the importance of the variance term
#'        (only currently available for the SCA method).
#' @param w_lb Lower bound (either a vector or a scalar) on the value of each
#'        portfolio weight.
#' @param w_ub Upper bound (either a vector or a scalar) on the value of each
#'        portfolio weight.
#' @param Cmat Equality constraints matrix.
#' @param cvec Equality constraints vector.
#' @param Dmat Inequality constraints matrix.
#' @param dvec Inequality constraints vector.
#' @param method_init Method to compute the vanilla solution. In case of
#'        additional constraints or objective terms, this solution is used as
#'        the initial point for the subsequent method. The default is
#'        \code{"cyclical-spinu"}. See details below.
#' @param method Method to solve the non-vanilla formulation. The default is \code{"sca"}.
#'        See details below. (DEPRECATED)
#' @param formulation String indicating the risk concentration formulation to be used.
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
#' @param w0 Initial value for the portfolio weights. Default is a convex
#'        combination of the risk parity portfolio, the (uncorrelated) minimum variance
#'        portfolio, and the maximum return portfolio.
#' @param theta0 Initial value for theta (in case formulation uses theta). If not provided,
#'        the optimum solution for a fixed vector of portfolio weights will be used.
#' @param gamma Learning rate for the SCA method.
#' @param zeta Factor used to decrease the learning rate at each iteration for the SCA method.
#' @param tau Regularization factor.
#' @param maxiter Maximum number of iterations for the SCA loop.
#' @param ftol Convergence tolerance on the objective function.
#' @param wtol Convergence tolerance on the values of the portfolio weights.
#' @param use_gradient This parameter is meaningful only if method is either
#'        \code{"alabama"} or \code{"slsqp"}. If \code{TRUE} (default value), analytical gradients of the
#'        objective function will be used (strongly recommended to achieve faster results).
#' @param use_qp_solver Whether or not to use the general QP solver from
#'        quadprog to solve each iteration of the SCA algorithm. Default is TRUE.
#' @return A list containing possibly the following elements:
#' \item{\code{w}}{Optimal portfolio vector.}
#' \item{\code{relative_risk_contribution}}{The relative risk contribution of every asset.}
#' \item{\code{theta}}{Optimal value for theta (in case that it is part of
#'                     the chosen formulation.)}
#' \item{\code{obj_fun}}{Sequence of values of the objective function at each iteration.}
#' \item{\code{risk_concentration}}{Risk concentration term of the portfolio \code{R(w)}.}
#' \item{\code{mean_return}}{Expected return term of the portoflio \code{t(w)\%*\%mu},
#'                           if the term is included in the optimization.}
#' \item{\code{variance}}{Variance term of the portfolio \code{t(w)\%*\%Sigma\%*\%w},
#'                        if the term is included in the optimization.}
#' \item{\code{elapsed_time}}{Elapsed time recorded at every iteration.}
#' \item{\code{convergence}}{Boolean flag to indicate whether or not the optimization converged.}
#' \item{\code{is_feasible}}{Boolean flag to indicate whether or not the computed portfolio respects the linear constraints.}
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
#' res$relative_risk_contribution
#' #> [1] 0.2 0.2 0.2 0.2 0.2
#'
#' # risk budggeting portfolio
#' res <- riskParityPortfolio(Sigma, b = c(0.4, 0.4, 0.1, 0.05, 0.05))
#' res$relative_risk_contribution
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
                                Cmat = NULL, cvec = NULL,
                                Dmat = NULL, dvec = NULL,
                                method_init = c("cyclical-spinu", "cyclical-roncalli", "newton"),
                                method = c("sca", "alabama", "slsqp"),
                                formulation = NULL, w0 = NULL, theta0 = NULL,
                                gamma = .9, zeta = 1e-7, tau = NULL,
                                maxiter = 1000, ftol = 1e-8, wtol = .5e-6,
                                use_gradient = TRUE, use_qp_solver = TRUE) {
  N <- nrow(Sigma)
  stocks_names <- colnames(Sigma)
  # check that constraints are consistent
  is_Cmat_null <- is.null(Cmat)
  is_cvec_null <- is.null(cvec)
  if ((!is_Cmat_null) && is_cvec_null) stop("Matrix Cmat has been given, but vector cvec is NULL.")
  if (is_Cmat_null && (!is_cvec_null)) stop("Vector cvec has been given, but matrix Cmat is NULL.")
  if (is_Cmat_null)
    Cmat <- matrix(1, 1, N)
  else
    Cmat <- rbind(Cmat, matrix(1, 1, N))
  if (is_cvec_null)
    cvec <- c(1)
  else
    cvec <- c(cvec, 1)
  if (nrow(Cmat) > 1)
    if (Matrix::rankMatrix(Cmat) != nrow(Cmat)) stop("Cmat contains linearly dependent rows.")
  if (nrow(Cmat) != length(cvec)) stop("Shapes of Cmat and cvec are inconsistent.")
  is_Dmat_null <- is.null(Dmat)
  is_dvec_null <- is.null(dvec)
  if ((!is_Dmat_null) && is_dvec_null) stop("Matrix Dmat has been given, but vector dvec is NULL.")
  if (is_Dmat_null && (!is_dvec_null)) stop("Vector dvec has been given, but matrix Dmat is NULL.")
  # default values
  if (is.null(b))
    b <- rep(1/N, N)
  else if (sum(b) != 1)
    warning("Budget vector b does not sum up to 1.")

  if (length(w_ub) == 1) w_ub <- rep(w_ub, N)
  if (length(w_lb) == 1) w_lb <- rep(w_lb, N)
  has_only_equality_constraints <- all(w_lb == (-Inf), w_ub == Inf)
  if (!has_only_equality_constraints) {
    In <- diag(nrow(Sigma))
    if (is.null(Dmat))
      Dmat <- rbind(-In, In)
    else
      Dmat <- rbind(Dmat, -In, In)
    if (is.null(dvec))
      dvec <- c(-w_lb, w_ub)
    else
      dvec <- c(dvec, -w_lb, w_ub)
    if (nrow(Dmat) != length(dvec)) stop("Shapes of Dmat and dvec are inconsistent.")
  }
  # define boolean cases for easier reading of code
  has_mu <- !is.null(mu)
  has_theta <- !is.null(theta0)
  has_var <- lmd_var > 0
  has_formulation <- !is.null(formulation)
  has_fancy_box <- any(w_lb != 0) || any(w_ub != 1)
  has_initial_point <- !is.null(w0)
  has_equality_constraints <- nrow(Cmat) > 1
  has_inequality_constraints <- !is_Dmat_null
  has_std_constraints <- !has_only_equality_constraints && (nrow(Cmat) == 1) && is_Dmat_null
  is_vanilla_formulation <- !(has_mu || has_theta || has_var || has_fancy_box ||
                              has_only_equality_constraints || has_equality_constraints ||
                              has_inequality_constraints)

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
    if (!isFeasiblePortfolio(w0, Cmat, cvec, Dmat, dvec)) {
      if (has_initial_point) warning("Initial portfolio is unfeasible, projecting it onto the feasible set.")
      if (has_only_equality_constraints) {
        w0 <- project_onto_eq_and_ineq_constraint_set(w0 = w0, Cmat = Cmat,
                                                      cvec = cvec, Dmat = Dmat,
                                                      dvec = dvec)
      } else if(has_inequality_constraints) {
        w0 <- project_onto_eq_and_ineq_constraint_set(w0 = w0, Cmat = Cmat, cvec = cvec,
                                                      Dmat = Dmat, dvec = dvec)
      } else {
        w0 <- projectBudgetLineAndBox(w0 = w0, w_lb = w_lb, w_ub = w_ub)
      }
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
             if (has_only_equality_constraints || has_inequality_constraints)
               stop("General linear constraints not supported for method ", method)
             portfolio <- riskParityPortfolioGenSolver(Sigma = Sigma, b = b, mu = mu, lmd_mu = lmd_mu,
                                                       formulation = formulation, method = method,
                                                       use_gradient = use_gradient, w0 = w0, theta0 = theta0,
                                                       maxiter = maxiter, ftol = ftol, wtol = wtol)
           },
           stop("method ", method, " is not included.")
    )
  }
  names(portfolio$w) <- names(portfolio$relative_risk_contribution) <- stocks_names
  portfolio$is_feasible <- isFeasiblePortfolio(portfolio$w, Cmat, cvec, Dmat, dvec)
  return(portfolio)
}
