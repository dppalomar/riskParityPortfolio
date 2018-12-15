#include "objfunctions.h"
// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace std;


double obj_function_spinu(const Eigen::MatrixXd& Sigma, const Eigen::VectorXd& x,
                          const Eigen::VectorXd& b) {
  return (.5 * (x.transpose() * Sigma * x)).sum() - (b.array() * (x.array().log())).sum();
}


double obj_function_roncalli(const Eigen::MatrixXd& Sigma, const Eigen::VectorXd& x,
                             const Eigen::VectorXd& b) {
  return std::sqrt((x.transpose() * Sigma * x).sum()) - (b.array() * (x.array().log())).sum();
}
