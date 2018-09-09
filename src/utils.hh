#ifndef UTILS_H
#define UTILS_H
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Eigen;

MatrixXd compute_A_double_index(const VectorXd&, const MatrixXd&, const int);
#endif
