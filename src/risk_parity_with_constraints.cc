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
Eigen::VectorXd rpp_eq_and_ineq_constraints_iteration(const Eigen::MatrixXd& Cmat,
                                                      const Eigen::VectorXd& cvec,
                                                      const Eigen::MatrixXd& Dmat,
                                                      const Eigen::VectorXd& dvec,
                                                      const Eigen::MatrixXd& Qk,
                                                      const Eigen::VectorXd& qk,
                                                      const Eigen::VectorXd& wk,
                                                      const unsigned int maxiter) {

  unsigned int n = Cmat.cols();
  Eigen::VectorXd w_tilde(n), w_prev(n), w_tilde_bar(n);
  Eigen::VectorXd lmd = Eigen::VectorXd::Zeros(n);
  Eigen::VectorXd lmd_next = Eigen::VectorXd::Zeros(n);
  Eigen::VectorXd lmd_prev = Eigen::VectorXd::Zeros(n);
  Eigen::VectorXd mu = Eigen::VectorXd::Zeros(n);
  Eigen::VectorXd mu_next = Eigen::VectorXd::Zeros(n);
  Eigen::VectorXd mu_prev = Eigen::VectorXd::Zeros(n);
  LLT<MatrixXd> lltOfQk(Qk);
  Eigen::MatrixXd B(Cmat.cols() + Dmat.cols(), Cmat.rows());
  B << Cmat.transpose(), Dmat.transpose();
  double LC = (B * lltOfQk.solve(B.transpose())).norm();
  w_prev = wk;
  for (unsigned int i = 0; i < maxiter; ++i) {
    w_tilde = -lltOfQk.solve(qk + Cmat.transpose() * lmd + Dmat.transpose() * mu);
    w_tilde_bar = w_tilde + (i - 1) / (i + 2) * (w_tilde - w_prev);
    lmd_next = lmd + (i - 1) / (i + 2) * (lmd - lmd_prev) + (Cmat * w_tilde_bar - cvec) / LC;
    mu_next = (mu + (i - 1) / (i + 2) * (mu - mu_prev) + (Dmat * w_tilde_bar - dvec) / LC).maxElementWise;
    mu_prev = mu;
    mu = mu_next;
    lmd_prev = lmd;
    lmd = lmd_next;
  }
}
