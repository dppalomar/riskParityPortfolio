#include "cyclical_coordinate_descent.h"
// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace std;

// Cyclical coordinate descent for Spinu's formulation
// of the risk parity portfolio problem
// [[Rcpp::export]]
Eigen::VectorXd risk_parity_portfolio_ccd_spinu(const Eigen::MatrixXd& Sigma,
                                                const Eigen::VectorXd& b,
                                                const double tol,
                                                const unsigned int maxiter) {
  double aux, x_diff, xk_sum;
  const unsigned int N = b.size();
  Eigen::VectorXd xk = Eigen::VectorXd::Constant(N, 1);
  Eigen::VectorXd x_star(N), Sigma_xk(N), rc(N);
  xk = (1 / Sigma.sum()) * xk;
  Sigma_xk = Sigma * xk;
  for (unsigned int k = 0; k < maxiter; ++k) {
    for (unsigned int i = 0; i < N; ++i) {
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
  const unsigned int N = b.size();
  Eigen::VectorXd xk = Eigen::VectorXd::Constant(N, 1);
  Eigen::VectorXd x_star(N), Sigma_xk(N), rc(N);
  xk = (1 / Sigma.sum()) * xk;
  Sigma_xk = Sigma * xk;
  sigma = std::sqrt(xk.transpose() * Sigma * xk);
  for (unsigned int k = 0; k < maxiter; ++k) {
    for (unsigned int i = 0; i < N; ++i) {
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
