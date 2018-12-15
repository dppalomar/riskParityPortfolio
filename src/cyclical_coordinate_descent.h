#ifndef CYCLICAL_COORDINATE_DESCENT_H
#define CYCLICAL_COORDINATE_DESCENT_H
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <cmath>

using namespace Eigen;
VectorXd risk_parity_portfolio_ccd(const MatrixXd&, const VectorXd&, const double, const int);
double obj_function_roncalli(const MatrixXd&, const VectorXd&, const VectorXd&);
#endif
