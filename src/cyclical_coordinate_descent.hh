#ifndef OPERATORS_H
#define OPERATORS_H
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <cmath>

using namespace Eigen;
VectorXd risk_parity_portfolio_ccd(const MatrixXd&, const VectorXd&, const double, const int);
double obj_function(const MatrixXd&, const VectorXd&, const VectorXd&);
#endif
