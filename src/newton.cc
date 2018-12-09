#include "newton.hh"
#include "newton_nesterov.hh"
// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace std;

// [[Rcpp::export]]
Eigen::VectorXd risk_parity_portfolio_newton(const Eigen::MatrixXd& Sigma,
                                             const Eigen::VectorXd& b,
                                             const double tol,
                                             const int maxiter) {
  int N = b.size();
  Eigen::VectorXd xk = Eigen::VectorXd::Constant(N, 1);
  Eigen::VectorXd uk(N), dx(N), Sigma_dx(N);
  Eigen::MatrixXd Hk(N, N);
  double lambda;

  // initial guess
  xk = std::sqrt(b.sum() / Sigma.sum()) * xk;
  // Newton method
  for (int i = 0; i < maxiter; ++i) {
    uk = gradient_log_formulation(Sigma, xk, b);
    Hk = hessian_log_formulation(Sigma, xk, b);
    dx = -Hk.llt().solve(uk);
    lambda = -uk.dot(dx);
    Sigma_dx = Sigma * dx;
    if (lambda < tol)
      break;
    xk = xk + dx;
  }
  return xk / xk.sum();
}
