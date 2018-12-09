#ifndef NEWTON_H
#define NEWTON_H
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <cmath>

using namespace Eigen;
VectorXd risk_parity_portfolio_newton(const MatrixXd&, const VectorXd&, const double, const int);
#endif
