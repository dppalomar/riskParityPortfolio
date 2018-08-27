#' Jointly estimates a portfolio with risk-parity strategy and asset selection
#'
#' @param mu mean vector of the random returns of n assets
#' @param Sigma positive definite covariance matrix of the random returns
#'        of n assets
#' @param w0 initial estimate of the portfolio weights
#' @param risk_contrib function that takes the portfolio weights and the
#'        covariance matrix of the returns as inputs and outputs the risk
#'        contribution of each asset
#' @param gamma learning rate for the optimization
#' @param zeta how much to decrease the learning rate per iteration such that
#'        gamma <- gamma * (1 - zeta * gamma)
#' @param nu regularization parameter
#' @param l1 regularization parameter
#' @param l2 regularization parameter
#' @param tau regularization parameter
#' @param type approximation type
#' @param maxiter maximum number of iterations
#' @param w_tol convergence tolerance on the portfolio weights
#' @param theta_tol convergence tolerance on the average of the RCs for the
#'        selected assets
#' @param ftol convergence tolerance on the objective function
#' @export
riskParityPortfolio <- function(mu, Sigma, w0 = NA, risk_contrib = NA,
                                gamma = 1e-1, zeta = 1e-2, nu = 5e-1,
                                l1 = 1e-1, l2 = 4, tau = 1e-3, type = "1",
                                p = 2e-3, e = 1e-8, maxiter = 5000,
                                w_tol = 1e-4, theta_tol = 1e-4, ftol = 1e-4) {

  if (is.na(risk_contrib)) {
    # g is defined in successiveConvexApprox.R
    risk_contrib = g
  }

  if (any(is.na(w0))) {
    w0 <- 1 / diag(Sigma)
    w0 <- w0 / sum(w0)
  }

  # initialization
  rho_sq <- rho(w0, p, e) ^ 2
  x <- rho_sq / sum(rho_sq ^ 2)
  theta0 <- sum(x * risk_contrib(w0, Sigma))

  nll0 <- negLogLikelihood(w0, nu, mu, Sigma)
  nll_seq <- c(nll0)
  fun0 <- nll0 + negLogPrior(w0, w0, theta0, Sigma, l1, l2, p, e, tau, type)
  fun_seq <- c(fun0)

  for (k in 1:maxiter) {
    # update theta
    theta <- theta_update(theta0, w0, Sigma, risk_contrib, p, e, gamma)
    # update w
    w <- w_update(w0, theta, nu, gamma, l1, l2, mu, Sigma, tau, p, e, type)
    # check tolerance on parameters
    w_err <- norm(w - w0, "2") / max(1., norm(w, "2"))
    theta_err <- abs(theta - theta0) / max(1., theta)
    if ((w_err < w_tol) & (theta_err < theta_tol))
      break
    # check tolerance on objective function
    nll <- negLogLikelihood(w, nu, mu, Sigma)
    nll_seq <- c(nll_seq, nll)
    fun <- nll + negLogPrior(w, w0, theta, Sigma, l1, l2, p, e, tau, type)
    fun_seq <- c(fun_seq, fun)
    ferr <- abs(fun - fun0) / max(1., abs(fun))
    if (ferr < ftol)
      break

    gamma <- gamma * (1 - zeta * gamma)
    w0 <- w
    theta0 <- theta
    fun0 <- fun
  }

  return(list(portfolio_weights = w, negloglike = nll_seq, obj_fun = fun_seq))
}
