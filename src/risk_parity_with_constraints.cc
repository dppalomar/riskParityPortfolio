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
Eigen::VectorXd project_onto_equality_constraint_set(const Eigen::VectorXd& w,
                                                     const Eigen::MatrixXd& Cmat,
                                                     const Eigen::VectorXd& cvec) {
  LLT<MatrixXd> lltOfCCt(Cmat * Cmat.transpose());
  return w - Cmat.transpose() * lltOfCCt.solve(Cmat * w - cvec);
}

// [[Rcpp::export]]
std::vector<Eigen::VectorXd>
rpp_eq_and_ineq_constraints_iteration(const Eigen::MatrixXd& Cmat, const Eigen::VectorXd& cvec,
                                      const Eigen::MatrixXd& Dmat, const Eigen::VectorXd& dvec,
                                      const Eigen::MatrixXd& Qk, const Eigen::VectorXd& qk,
                                      const Eigen::VectorXd& wk, Eigen::VectorXd& chi,
                                      Eigen::VectorXd& chi_prev,
                                      Eigen::VectorXd& xi, Eigen::VectorXd& xi_prev) {

  std::vector<Eigen::VectorXd> params;
  unsigned int n = Cmat.cols();
  Eigen::VectorXd w_tilde(n), w_prev(n), w_tilde_bar(n),
                  chi_next(n), xi_next(n);
  LLT<MatrixXd> lltOfQk(Qk);
  Eigen::MatrixXd B(Cmat.rows() + Dmat.rows(), Cmat.cols());
  B << Cmat, Dmat;
  double LC = (B * lltOfQk.solve(B.transpose())).norm(), fac;
  w_prev = wk;
  unsigned int i = 1;
  while (true) {
    fac = (i - 1.)/(i + 2.);
    w_tilde = -lltOfQk.solve(qk + Cmat.transpose() * xi + Dmat.transpose() * chi);
    w_tilde_bar = w_tilde +  fac * (w_tilde - w_prev);
    xi_next = xi + fac * (xi - xi_prev) + (Cmat * w_tilde_bar - cvec) / LC;
    chi_next = (chi + fac * (chi - chi_prev) + (Dmat * w_tilde_bar - dvec) / LC).array().max(0);
    chi_prev = chi;
    xi_prev = xi;
    chi = chi_next;
    xi = xi_next;
    if(((w_tilde - w_prev).array().abs() <= .5e-5 * (w_tilde.array().abs() + w_prev.array().abs())).all())
      break;
    w_prev = w_tilde;
    ++i;
  }
  params.push_back(chi);
  params.push_back(chi_next);
  params.push_back(xi);
  params.push_back(xi_next);
  params.push_back(w_tilde);
  return params;
}
