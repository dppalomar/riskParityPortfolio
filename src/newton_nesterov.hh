#ifndef OPERATORS_H
#define OPERATORS_H
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <cmath>

using namespace Eigen;
VectorXd risk_parity_portfolio_nn(const MatrixXd&, const VectorXd&, const double, const int);
VectorXd gradient_log_formulation(const MatrixXd&, const VectorXd&, const VectorXd&);
MatrixXd hessian_log_formulation(const MatrixXd&, const VectorXd&, const VectorXd&);
double obj_function_log_formulation(const MatrixXd&, const VectorXd&, const VectorXd&);
#endif
