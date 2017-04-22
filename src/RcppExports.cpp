// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// rcpp_hello
List rcpp_hello();
RcppExport SEXP polenta_rcpp_hello() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello());
    return rcpp_result_gen;
END_RCPP
}
// Cmatrix
NumericMatrix Cmatrix(NumericMatrix msa);
RcppExport SEXP polenta_Cmatrix(SEXP msaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type msa(msaSEXP);
    rcpp_result_gen = Rcpp::wrap(Cmatrix(msa));
    return rcpp_result_gen;
END_RCPP
}
// nChoosek
int nChoosek(int n, int k);
RcppExport SEXP polenta_nChoosek(SEXP nSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(nChoosek(n, k));
    return rcpp_result_gen;
END_RCPP
}
// which_true
int which_true(LogicalVector x);
RcppExport SEXP polenta_which_true(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< LogicalVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(which_true(x));
    return rcpp_result_gen;
END_RCPP
}
// add_msa
NumericMatrix add_msa(NumericMatrix ref, NumericMatrix com);
RcppExport SEXP polenta_add_msa(SEXP refSEXP, SEXP comSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type ref(refSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type com(comSEXP);
    rcpp_result_gen = Rcpp::wrap(add_msa(ref, com));
    return rcpp_result_gen;
END_RCPP
}
