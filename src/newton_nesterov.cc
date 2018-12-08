#include "newton_nesterov.hh"
// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace std;

// [[Rcpp::export]]
Eigen::VectorXd risk_parity_portfolio_nn(const Eigen::MatrixXd& Sigma,
                                         const Eigen::VectorXd& b,
                                         const double tol,
                                         const int maxiter) {
  int N = b.size();
  Eigen::VectorXd xk = Eigen::VectorXd::Constant(N, 1);
  Eigen::VectorXd uk(N), d(N);
  Eigen::MatrixXd Hk(N, N);
  double dx, lambdak, lambda_star = 0.3628676;

  // initial guess
  xk = std::sqrt(b.sum() / Sigma.sum()) * xk;
  // damped phase
  for (int i = 0; i < maxiter; ++i) {
    uk = gradient_log_formulation(Sigma, xk, b);
    Hk = hessian_log_formulation(Sigma, xk, b);
    d = Hk.llt().solve(uk);
    dx = (d.array() / xk.array()).maxCoeff();
    lambdak = std::sqrt(uk.dot(d));
    xk = xk - d / (1 + dx);
    if (lambdak < lambda_star)
      break;
  }
  // quadratic phase
  for(int i = 0; i < maxiter; ++i) {
    uk = gradient_log_formulation(Sigma, xk, b);
    Hk = hessian_log_formulation(Sigma, xk, b);
    d = Hk.llt().solve(uk);
    lambdak = std::sqrt(uk.dot(d));
    xk = xk - d;
    if (lambdak < tol)
      break;
  }
  return xk / xk.sum();
}


Eigen::VectorXd gradient_log_formulation(const Eigen::MatrixXd& Sigma,
                                         const Eigen::VectorXd& xk,
                                         const Eigen::VectorXd& b) {
  return Sigma * xk - (b.array() / xk.array()).matrix();
}


Eigen::MatrixXd hessian_log_formulation(const Eigen::MatrixXd& Sigma,
                                        const Eigen::VectorXd& xk,
                                        const Eigen::VectorXd& b) {
  Eigen::MatrixXd M = (b.array() / (xk.array() * xk.array())).matrix().asDiagonal();
  return Sigma + M;
}


double obj_function_spinu(const Eigen::MatrixXd& Sigma,
                    const Eigen::VectorXd& xk,
                    const Eigen::VectorXd& b) {
  return (.5 * (xk.transpose() * Sigma * xk)).sum() - (b.array() * (xk.array().log())).sum();
}
