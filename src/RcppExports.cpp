// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// risk_parity_portfolio_ccd
Eigen::VectorXd risk_parity_portfolio_ccd(const Eigen::MatrixXd& Sigma, const Eigen::VectorXd& b, const double tol, const unsigned int maxiter);
RcppExport SEXP _riskParityPortfolio_risk_parity_portfolio_ccd(SEXP SigmaSEXP, SEXP bSEXP, SEXP tolSEXP, SEXP maxiterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type maxiter(maxiterSEXP);
    rcpp_result_gen = Rcpp::wrap(risk_parity_portfolio_ccd(Sigma, b, tol, maxiter));
    return rcpp_result_gen;
END_RCPP
}
// risk_parity_portfolio_nn
Eigen::VectorXd risk_parity_portfolio_nn(const Eigen::MatrixXd& Sigma, const Eigen::VectorXd& b, const double tol, const unsigned int maxiter);
RcppExport SEXP _riskParityPortfolio_risk_parity_portfolio_nn(SEXP SigmaSEXP, SEXP bSEXP, SEXP tolSEXP, SEXP maxiterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type maxiter(maxiterSEXP);
    rcpp_result_gen = Rcpp::wrap(risk_parity_portfolio_nn(Sigma, b, tol, maxiter));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_riskParityPortfolio_risk_parity_portfolio_ccd", (DL_FUNC) &_riskParityPortfolio_risk_parity_portfolio_ccd, 4},
    {"_riskParityPortfolio_risk_parity_portfolio_nn", (DL_FUNC) &_riskParityPortfolio_risk_parity_portfolio_nn, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_riskParityPortfolio(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
