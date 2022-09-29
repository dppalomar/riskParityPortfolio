#include "cyclical_coordinate_descent.h"
// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace std;

// cyclical coordinate descent algo by Choi & Chen 2022
// ref: https://arxiv.org/pdf/2203.00148.pdf
// [[Rcpp::export]]
Eigen::VectorXd risk_parity_portfolio_ccd_choi(const Eigen::VectorXd& cov,
                                        const Eigen::VectorXd& b,
                                        const double tol = 1E-4,
                                        const unsigned int maxiter = 100) {
  const unsigned int n = b.size();
  Eigen::VectorXd a(n);
  Eigen::VectorXd vol = cov.diagonal().array().sqrt();
  Eigen::VectorXd invvol = (1 / vol.array()).matrix();
  Eigen::MatrixXd corr = cov.array().colwise() * invvol.array();
  corr = corr.array().rowwise() * invvol.transpose().array();
  Eigen::MatrixXd adj = corr;
  adj.diagonal().array() = 0;
  Eigen::VectorXd wk = Eigen::VectorXd::Ones(n);
  wk = (wk.array() / std::sqrt(corr.sum())).matrix();
  for (unsigned int k = 0; k < maxiter; ++k) {
    // compute portfolio weights
    a = 0.5 * adj * wk;
    wk = ((a.array() * a.array() + b.array()).sqrt() - a.array()).matrix();
    if ((wk.array() * (corr * wk).array() - b.array()).abs().maxCoeff() < tol)
      break;
  }
  Eigen::VectorXd w = wk.array() / vol.array();
  return (w / w.sum()).matrix();
}

// Cyclical coordinate descent for Spinu's formulation
// of the risk parity portfolio problem
// [[Rcpp::export]]
Eigen::VectorXd risk_parity_portfolio_ccd_spinu(const Eigen::MatrixXd& Sigma,
                                                const Eigen::VectorXd& b,
                                                const double tol,
                                                const unsigned int maxiter) {
  double aux, x_diff, xk_sum;
  const unsigned int n = b.size();
  Eigen::VectorXd xk = Eigen::VectorXd::Constant(n, 1);
  Eigen::VectorXd x_star(n), Sigma_xk(n), rc(n);
  xk = (1 / Sigma.sum()) * xk;
  Sigma_xk = Sigma * xk;
  for (unsigned int k = 0; k < maxiter; ++k) {
    for (unsigned int i = 0; i < n; ++i) {
      // compute update for the portfolio weights x
      aux = xk(i) * Sigma(i, i) - Sigma_xk(i);
      x_star(i) = (.5 / Sigma(i, i)) * (aux + std::sqrt(aux * aux + 4 * Sigma(i, i) * b(i)));
      // update auxiliary terms
      x_diff = x_star(i) - xk(i);
      Sigma_xk += (Sigma.col(i).array() * x_diff).matrix();
      xk(i) = x_star(i);
    }
    xk_sum = xk.sum();
    rc = (xk.array() * (Sigma_xk).array() / (xk_sum * xk_sum)).matrix();
    if ((rc.array() / rc.sum() - b.array()).abs().maxCoeff() < tol)
      break;
  }
  return x_star / xk_sum;
}

// Cyclical coordinate descent for Roncalli's square-root formulation
// of the risk parity portfolio problem
// [[Rcpp::export]]
Eigen::VectorXd risk_parity_portfolio_ccd_roncalli(const Eigen::MatrixXd& Sigma,
                                                   const Eigen::VectorXd& b,
                                                   const double tol,
                                                   const unsigned int maxiter) {
  double aux, sigma, x_diff, xk_sum;
  const unsigned int n = b.size();
  Eigen::VectorXd xk = Eigen::VectorXd::Constant(n, 1);
  Eigen::VectorXd x_star(n), Sigma_xk(n), rc(n);
  xk = (1 / Sigma.sum()) * xk;
  Sigma_xk = Sigma * xk;
  sigma = std::sqrt(xk.transpose() * Sigma * xk);
  for (unsigned int k = 0; k < maxiter; ++k) {
    for (unsigned int i = 0; i < n; ++i) {
      // compute update for the portfolio weights x
      aux = xk(i) * Sigma(i, i) - Sigma_xk(i);
      x_star(i) = (.5 / Sigma(i, i)) * (aux + std::sqrt(aux * aux + 4 * Sigma(i, i) * b(i) * sigma));
      // update auxiliary terms
      x_diff = x_star(i) - xk(i);
      Sigma_xk += (Sigma.col(i).array() * x_diff).matrix();
      sigma = std::sqrt(sigma * sigma + 2 * x_diff * (Sigma.row(i).array() * xk.transpose().array()).sum() +
                        Sigma(i, i) * x_diff * x_diff);
      xk(i) = x_star(i);
    }
    xk_sum = xk.sum();
    rc = (xk.array() * (Sigma_xk).array() / (xk_sum * xk_sum)).matrix();
    if ((rc.array() / rc.sum() - b.array()).abs().maxCoeff() < tol)
      break;
  }
  return x_star / xk_sum;
}


// Cyclical coordinate descent for the active risk parity portfolio
// proposed by Roncalli
// [[Rcpp::export]]
Eigen::VectorXd active_risk_parity_portfolio_ccd(const Eigen::MatrixXd& Sigma,
                                                 const Eigen::VectorXd& b,
                                                 const Eigen::VectorXd& mu,
                                                 const double c, const double r,
                                                 const double tol,
                                                 const unsigned int maxiter) {
  double aux, sigma, x_diff, xk_sum;
  const unsigned int n = b.size();
  Eigen::VectorXd xk = Eigen::VectorXd::Constant(n, 1);
  Eigen::VectorXd x_star(n), Sigma_xk(n), pi(n), rc(n);
  xk = (1 / Sigma.sum()) * xk;
  Sigma_xk = Sigma * xk;
  pi.array() = mu.array() - r;
  sigma = std::sqrt(xk.transpose() * Sigma * xk);
  for (unsigned int k = 0; k < maxiter; ++k) {
    for (unsigned int i = 0; i < n; ++i) {
      // compute update for the portfolio weights x
      aux = c * (xk(i) * Sigma(i, i) - Sigma_xk(i)) + pi(i) * sigma;
      x_star(i) = (aux + std::sqrt(aux * aux + 4 * c * Sigma(i, i) * b(i) * sigma)) / (2 * c * Sigma(i, i));
      x_diff = x_star(i) - xk(i);
      // update auxiliary terms
      Sigma_xk += (Sigma.col(i).array() * x_diff).matrix();
      sigma = std::sqrt(sigma * sigma + 2 * x_diff * (Sigma.row(i).array() * xk.transpose().array()).sum() +
                        Sigma(i, i) * x_diff * x_diff);
      xk(i) = x_star(i);
    }
    xk_sum = xk.sum();
    rc = (xk.array() * (Sigma_xk).array() / (xk_sum * xk_sum)).matrix();
    if ((rc.array() / rc.sum() - b.array()).abs().maxCoeff() < tol)
      break;
  }
  return x_star / xk_sum;
}
