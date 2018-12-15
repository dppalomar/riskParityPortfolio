#include "newton_nesterov.h"
// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace std;

// [[Rcpp::export]]
Eigen::VectorXd risk_parity_portfolio_nn(const Eigen::MatrixXd& Sigma,
                                         const Eigen::VectorXd& b,
                                         const double tol,
                                         const unsigned int maxiter) {
  const unsigned int N = b.size();
  Eigen::VectorXd xk = Eigen::VectorXd::Constant(N, 1);
  Eigen::VectorXd uk(N), d(N), rc(N);
  Eigen::MatrixXd Hk(N, N);
  double dx, lambdak, xk_sum, lambda_star = 0.3628676;

  // initial guess
  xk = std::sqrt(b.sum() / Sigma.sum()) * xk;
  // damped phase
  for (unsigned int i = 0; i < maxiter; ++i) {
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
  for(unsigned int i = 0; i < maxiter; ++i) {
    uk = gradient_log_formulation(Sigma, xk, b);
    Hk = hessian_log_formulation(Sigma, xk, b);
    d = Hk.llt().solve(uk);
    //lambdak = std::sqrt(uk.dot(d));
    xk = xk - d;
    xk_sum = xk.sum();
    rc = (xk.array() * (Sigma * xk).array() / (xk_sum * xk_sum)).matrix();
    if ((rc.array() / rc.sum() - b.array()).abs().maxCoeff() < tol)
      break;
    //if (lambdak < tol)
    //  break;
  }
  return xk / xk_sum;
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
