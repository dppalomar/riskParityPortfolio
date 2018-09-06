#include "utils.hh"
// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;

//' @export
// [[Rcpp::export]]
Eigen::MatrixXd computeACpp(const Eigen::VectorXd& w, const Eigen::MatrixXd& Sigma) {
  const int N = w.size();
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N * N, N);

  for (int i = 0; i < N; ++i) {
    Eigen::MatrixXd Mi = Eigen::MatrixXd::Zero(N, N);
    Mi.row(i) = Sigma.row(i);
    for (int j = i+1; j < N; ++j) {
      if (i == j)
        continue;
      Eigen::MatrixXd Mj = Eigen::MatrixXd::Zero(N, N);
      Mj.row(j) = Sigma.row(j);
      A.row(i + j * N) = (Mi + Mi.transpose() - Mj - Mj.transpose()) * w;
      A.row(j + i * N) = A.row(i + j * N);
    }
  }

  return A;
}
