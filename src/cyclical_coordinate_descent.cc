#include "cyclical_coordinate_descent.hh"
// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace std;

// [[Rcpp::export]]
Eigen::VectorXd risk_parity_portfolio_ccd(const Eigen::MatrixXd& Sigma,
                                          const Eigen::VectorXd& b,
                                          const double tol,
                                          const int maxiter) {
  double aux, sigma, x_diff;
  int N = b.size();
  Eigen::VectorXd xk = Eigen::VectorXd::Constant(N, 1);
  Eigen::VectorXd x_star(N);
  Eigen::VectorXd Sigma_xk(N);
  xk = (1 / Sigma.sum()) * xk;
  Sigma_xk = Sigma * xk;
  sigma = std::sqrt(xk.transpose() * Sigma * xk);
  for (int k = 0; k < maxiter; ++k) {
    for (int i = 0; i < N; ++i) {
      // compute update for the portfolio weights x
      aux = xk(i) * Sigma(i, i) - Sigma_xk(i);
      x_star(i) = (.5 / Sigma(i, i)) * (aux + std::sqrt(aux * aux + 4 * Sigma(i, i) * b(i) * sigma));
      // update auxiliary terms
      x_diff = x_star(i) - xk(i);
      Sigma_xk += (Sigma.col(i).array() * x_diff).matrix();
      sigma = std::sqrt(sigma * sigma + (2 * (x_star(i) - xk(i)) *
                                         (Sigma.row(i).array() * xk.transpose().array()).sum())
                        + Sigma(i, i) * (x_diff * x_diff));
      xk(i) = x_star(i);
    }
    if ((Sigma_xk.array() - b.array()).abs().maxCoeff() < tol)
      break;
  }
  return x_star / x_star.sum();
}


double obj_function(const Eigen::MatrixXd& Sigma,
                    const Eigen::VectorXd& x,
                    const Eigen::VectorXd& b) {
  return std::sqrt((x.transpose() * Sigma * x).sum()) - (b.array() * (x.array().log())).sum();
}
