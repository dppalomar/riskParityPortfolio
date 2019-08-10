#ifndef RISK_PARITY_WITH_CONSTRAINTS_H
#define RISK_PARITY_WITH_CONSTRAINTS_H
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <cmath>

using namespace Eigen;
VectorXd rpp_equality_constraints_iteration(const Eigen::MatrixXd&, const Eigen::VectorXd&,
                                            const Eigen::MatrixXd&, const Eigen::VectorXd&);
Eigen::VectorXd project_onto_equality_constraint_set(const Eigen::VectorXd&,
                                                     const Eigen::MatrixXd&,
                                                     const Eigen::VectorXd&);
std::vector<Eigen::VectorXd>
rpp_eq_and_ineq_constraints_iteration(const Eigen::MatrixXd&, const Eigen::VectorXd&,
                                      const Eigen::MatrixXd&, const Eigen::VectorXd&,
                                      const Eigen::MatrixXd&, const Eigen::VectorXd&,
                                      const Eigen::VectorXd&, Eigen::VectorXd&,
                                      Eigen::VectorXd&, Eigen::VectorXd&,
                                      Eigen::VectorXd&, unsigned int, double tol);
#endif
