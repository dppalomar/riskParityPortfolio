#ifndef RISK_PARITY_WITH_CONSTRAINTS_H
#define RISK_PARITY_WITH_CONSTRAINTS_H
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <cmath>

using namespace Eigen;
VectorXd rpp_equality_constraints_iteration(const MatrixXd&, const VectorXd&,
                                            const MatrixXd&, const VectorXd&);
#endif
