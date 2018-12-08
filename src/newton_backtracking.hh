#ifndef NEWTON_BACKTRACKING_H
#define NEWTON_BACKTRACKING_H
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <cmath>

using namespace Eigen;
VectorXd risk_parity_portfolio_nn_bt(const MatrixXd&, const VectorXd&, const double, const int);
#endif
