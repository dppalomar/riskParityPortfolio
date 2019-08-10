#include "risk_parity_with_constraints.h"
// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace std;

// [[Rcpp::export]]
Eigen::VectorXd rpp_equality_constraints_iteration(const Eigen::MatrixXd& Cmat,
                                                   const Eigen::VectorXd& cvec,
                                                   const Eigen::MatrixXd& Qk,
                                                   const Eigen::VectorXd& qk) {
  LDLT<MatrixXd> ldltOfQk(Qk);
  Eigen::MatrixXd Vk = Cmat * ldltOfQk.solve(Cmat.transpose());
  ColPivHouseholderQR<MatrixXd> QRofVk(Vk);
  Eigen::VectorXd lambdak = -QRofVk.solve(Cmat * ldltOfQk.solve(qk) + cvec);
  return -ldltOfQk.solve(qk + Cmat.transpose() * lambdak);
}


// [[Rcpp::export]]
Eigen::VectorXd project_onto_equality_constraint_set(const Eigen::VectorXd& w,
                                                     const Eigen::MatrixXd& Cmat,
                                                     const Eigen::VectorXd& cvec) {
  ColPivHouseholderQR<MatrixXd> QRofCCt(Cmat * Cmat.transpose());
  return w - Cmat.transpose() * QRofCCt.solve(Cmat * w - cvec);
}

// [[Rcpp::export]]
std::vector<Eigen::VectorXd>
rpp_eq_and_ineq_constraints_iteration(const Eigen::MatrixXd& Cmat, const Eigen::VectorXd& cvec,
                                      const Eigen::MatrixXd& Dmat, const Eigen::VectorXd& dvec,
                                      const Eigen::MatrixXd& Qk, const Eigen::VectorXd& qk,
                                      const Eigen::VectorXd& wk,
                                      Eigen::VectorXd& dual_mu_0, Eigen::VectorXd& dual_mu_minus_1,
                                      Eigen::VectorXd& dual_lmd_0, Eigen::VectorXd& dual_lmd_minus_1,
                                      unsigned int maxiter, double tol) {

  std::vector<Eigen::VectorXd> params;
  unsigned int n = Cmat.cols();
  unsigned int dual_lmd_len = Cmat.rows();
  unsigned int dual_mu_len = Dmat.rows();
  Eigen::VectorXd w_tilde(n), w_prev(n), w_tilde_bar(n),
                  dual_mu_next(dual_mu_len), dual_lmd_next(dual_lmd_len);
  LLT<MatrixXd> lltOfQk(Qk);
  Eigen::MatrixXd B(dual_lmd_len + dual_mu_len, n);
  B << Cmat, Dmat;
  double LC = (B * lltOfQk.solve(B.transpose())).norm(), fac;
  w_prev = wk;
  for(unsigned int i = 0; i < maxiter; ++i) {
    fac = (i - 1.)/(i + 2.);
    w_tilde = -lltOfQk.solve(qk + Cmat.transpose() * dual_lmd_0 + Dmat.transpose() * dual_mu_0);
    w_tilde_bar = w_tilde +  fac * (w_tilde - w_prev);
    dual_lmd_next = dual_lmd_0 + fac * (dual_lmd_0 - dual_lmd_minus_1) + (Cmat * w_tilde_bar - cvec) / LC;
    dual_mu_next = (dual_mu_0 + fac * (dual_mu_0 - dual_mu_minus_1) + (Dmat * w_tilde_bar - dvec) / LC).array().max(0);
    dual_mu_minus_1 = dual_mu_0;
    dual_lmd_minus_1 = dual_lmd_0;
    dual_mu_0 = dual_mu_next;
    dual_lmd_0 = dual_lmd_next;
    if(((w_tilde - w_prev).array().abs() <=
        tol * (w_tilde.array().abs() + w_prev.array().abs())).all())
      break;
    w_prev = w_tilde;
  }
  params.push_back(dual_mu_0);
  params.push_back(dual_mu_next);
  params.push_back(dual_lmd_0);
  params.push_back(dual_lmd_next);
  params.push_back(w_tilde);
  return params;
}
