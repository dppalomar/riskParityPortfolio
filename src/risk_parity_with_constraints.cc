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
