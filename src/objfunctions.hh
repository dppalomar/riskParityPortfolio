#ifndef OBJFUNCTIONS_H
#define OBJFUNCTIONS_H
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <cmath>

using namespace Eigen;

double obj_function_spinu(const Eigen::MatrixXd&, const Eigen::VectorXd&,
                          const Eigen::VectorXd&);
double obj_function_roncalli(const Eigen::MatrixXd&, const Eigen::VectorXd&,
                             const Eigen::VectorXd&);
#endif
