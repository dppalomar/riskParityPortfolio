#ifndef CYCLICAL_COORDINATE_DESCENT_H
#define CYCLICAL_COORDINATE_DESCENT_H
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <cmath>

using namespace Eigen;
VectorXd risk_parity_portfolio_ccd_spinu(const MatrixXd&, const VectorXd&, const double, const int);
VectorXd risk_parity_portfolio_ccd_roncalli(const MatrixXd&, const VectorXd&, const double, const int);
VectorXd active_risk_parity_portfolio_ccd(const MatrixXd&, const VectorXd&, const VectorXd&,
                                          const double, const double, const double, const unsigned int);
#endif
