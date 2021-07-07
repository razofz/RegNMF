// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// CNmf
Rcpp::List CNmf(Eigen::Map<Eigen::MatrixXd> V, int K, int maxiter, Eigen::Map<Eigen::MatrixXd> W0, Eigen::Map<Eigen::MatrixXd> H0, int core);
RcppExport SEXP _packageTryV3_CNmf(SEXP VSEXP, SEXP KSEXP, SEXP maxiterSEXP, SEXP W0SEXP, SEXP H0SEXP, SEXP coreSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type V(VSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type W0(W0SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type H0(H0SEXP);
    Rcpp::traits::input_parameter< int >::type core(coreSEXP);
    rcpp_result_gen = Rcpp::wrap(CNmf(V, K, maxiter, W0, H0, core));
    return rcpp_result_gen;
END_RCPP
}
// CPPNMF_cluster_joint_cross_domain_try
Rcpp::List CPPNMF_cluster_joint_cross_domain_try(Eigen::Map<Eigen::MatrixXd> PeakO, Eigen::Map<Eigen::MatrixXd> X, Eigen::Map<Eigen::MatrixXd> Reg, int K, int maxiter, double lambda1, double lambda2, Eigen::Map<Eigen::MatrixXd> W10, Eigen::Map<Eigen::MatrixXd> W20, Eigen::Map<Eigen::MatrixXd> W30, Eigen::Map<Eigen::MatrixXd> H0, NumericVector c1, NumericVector c2, NumericVector Reg_w, int core);
RcppExport SEXP _packageTryV3_CPPNMF_cluster_joint_cross_domain_try(SEXP PeakOSEXP, SEXP XSEXP, SEXP RegSEXP, SEXP KSEXP, SEXP maxiterSEXP, SEXP lambda1SEXP, SEXP lambda2SEXP, SEXP W10SEXP, SEXP W20SEXP, SEXP W30SEXP, SEXP H0SEXP, SEXP c1SEXP, SEXP c2SEXP, SEXP Reg_wSEXP, SEXP coreSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type PeakO(PeakOSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type Reg(RegSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type lambda1(lambda1SEXP);
    Rcpp::traits::input_parameter< double >::type lambda2(lambda2SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type W10(W10SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type W20(W20SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type W30(W30SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type H0(H0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type c1(c1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type c2(c2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Reg_w(Reg_wSEXP);
    Rcpp::traits::input_parameter< int >::type core(coreSEXP);
    rcpp_result_gen = Rcpp::wrap(CPPNMF_cluster_joint_cross_domain_try(PeakO, X, Reg, K, maxiter, lambda1, lambda2, W10, W20, W30, H0, c1, c2, Reg_w, core));
    return rcpp_result_gen;
END_RCPP
}
// Fold_RE_TG_MultiAdjustCore
NumericMatrix Fold_RE_TG_MultiAdjustCore(NumericMatrix E2, NumericMatrix O2, NumericMatrix Symbol_location, NumericMatrix Peak_location);
RcppExport SEXP _packageTryV3_Fold_RE_TG_MultiAdjustCore(SEXP E2SEXP, SEXP O2SEXP, SEXP Symbol_locationSEXP, SEXP Peak_locationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type E2(E2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type O2(O2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Symbol_location(Symbol_locationSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Peak_location(Peak_locationSEXP);
    rcpp_result_gen = Rcpp::wrap(Fold_RE_TG_MultiAdjustCore(E2, O2, Symbol_location, Peak_location));
    return rcpp_result_gen;
END_RCPP
}
// eps
double eps(double a);
RcppExport SEXP _packageTryV3_eps(SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(eps(a));
    return rcpp_result_gen;
END_RCPP
}
// Cjaccard
NumericMatrix Cjaccard(NumericMatrix MM);
RcppExport SEXP _packageTryV3_Cjaccard(SEXP MMSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type MM(MMSEXP);
    rcpp_result_gen = Rcpp::wrap(Cjaccard(MM));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_packageTryV3_CNmf", (DL_FUNC) &_packageTryV3_CNmf, 6},
    {"_packageTryV3_CPPNMF_cluster_joint_cross_domain_try", (DL_FUNC) &_packageTryV3_CPPNMF_cluster_joint_cross_domain_try, 15},
    {"_packageTryV3_Fold_RE_TG_MultiAdjustCore", (DL_FUNC) &_packageTryV3_Fold_RE_TG_MultiAdjustCore, 4},
    {"_packageTryV3_eps", (DL_FUNC) &_packageTryV3_eps, 1},
    {"_packageTryV3_Cjaccard", (DL_FUNC) &_packageTryV3_Cjaccard, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_packageTryV3(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
