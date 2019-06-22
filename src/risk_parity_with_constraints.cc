#include "risk_parity_with_constraints.h"
// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace std;

// [[Rcpp::export]]
Eigen::VectorXd rpp_equality_constraints_iteration(const Eigen::MatrixXd& Cmat,
                                                   const Eigen::VectorXd& cvec,
                                                   const Eigen::MatrixXd& Qk,
                                                   const Eigen::VectorXd& qk) {
  LLT<MatrixXd> lltOfQk(Qk);
  Eigen::MatrixXd Vk = Cmat * lltOfQk.solve(Cmat.transpose());
  ColPivHouseholderQR<MatrixXd> QROfVk(Vk);
  Eigen::VectorXd lambdak = -QROfVk.solve(Cmat * lltOfQk.solve(qk) + cvec);
  return -lltOfQk.solve(qk + Cmat.transpose() * lambdak);
}

// [[Rcpp::export]]
std::vector<Eigen::VectorXd>
rpp_eq_and_ineq_constraints_iteration(const Eigen::MatrixXd& Cmat, const Eigen::VectorXd& cvec,
                                      const Eigen::MatrixXd& Dmat, const Eigen::VectorXd& dvec,
                                      const Eigen::MatrixXd& Qk, const Eigen::VectorXd& qk,
                                      const Eigen::VectorXd& wk, Eigen::VectorXd& mu,
                                      Eigen::VectorXd& mu_prev, Eigen::VectorXd& mu_next,
                                      Eigen::VectorXd& lmd, Eigen::VectorXd& lmd_prev,
                                      Eigen::VectorXd& lmd_next, const unsigned int maxiter) {

  unsigned int n = Cmat.cols();
  Eigen::VectorXd w_tilde(n), w_prev(n), w_tilde_bar(n);
  LLT<MatrixXd> lltOfQk(Qk);
  Eigen::MatrixXd B(Cmat.cols() + Dmat.cols(), Cmat.rows());
  B << Cmat.transpose(), Dmat.transpose();
  double LC = (B * lltOfQk.solve(B.transpose())).norm();
  w_prev = wk;
  for (unsigned int i = 0; i < maxiter; ++i) {
    w_tilde = -lltOfQk.solve(qk + Cmat.transpose() * lmd + Dmat.transpose() * mu);
    w_tilde_bar = w_tilde + (i - 1) / (i + 2) * (w_tilde - w_prev);
    lmd_next = lmd + (i - 1) / (i + 2) * (lmd - lmd_prev) + (Cmat * w_tilde_bar - cvec) / LC;
    mu_next = (mu + (i - 1) / (i + 2) * (mu - mu_prev) + (Dmat * w_tilde_bar - dvec) / LC).array().max(0);
    mu_prev = mu;
    lmd_prev = lmd;
    mu = mu_next;
    lmd = lmd_next;
    if(((w_tilde - w_prev).array().abs() <= .5e-4 * (w_tilde.array().abs() + w_prev.array().abs())).all())
      break;
    w_prev = w_tilde;
  }
  std::vector<Eigen::VectorXd> params = {mu, mu_prev, mu_next, lmd, lmd_prev, lmd_next, w_tilde};
  return params;
}
