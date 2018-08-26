riskParityPortfolio <- function (mu, Sigma, w0 = NA, g = NA, gamma = 1e-1,
                                 zeta = 1e-2, l1 = 1e-1, l2 = 4, tau = 1e-3,
                                 type = "1", maxiter = 5000, w_tol = 1e-4,
                                 theta_tol = 1e-4, ftol = 1e-4) {

  if (is.na(g)) {
    g = riskParityPortfolio::g
  }

  # initialization
  rho_sq <- rho(w0) ^ 2
  x <- rho_sq / sum(rho_sq ^ 2)
  theta <- sum(x * g(w0, Sigma))
  w <- w0

  w_seq <- c()
  fun_seq <- c()
  nll_seq <- c()

  for (i in 1:maxiter) {
    # update theta
    theta <- theta_update(theta, w, g, gamma = gamma)
    # update w
    w <- w_update(w, theta, nu = nu, gamma = gamma, l1 = l1, l2 = l2, mu = mu,
                  Sigma = Sigma, tau = tau, type = type)
    gamma <- gamma * (1 - zeta * gamma)
  }

  return(w)
}
