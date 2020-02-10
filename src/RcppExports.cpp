// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// solveBiCGSTAB
Eigen::VectorXd solveBiCGSTAB(Eigen::SparseMatrix<double>& lhs, Eigen::VectorXd& rhs, int maxit, double tol, bool print);
RcppExport SEXP _implicit_solveBiCGSTAB(SEXP lhsSEXP, SEXP rhsSEXP, SEXP maxitSEXP, SEXP tolSEXP, SEXP printSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double>& >::type lhs(lhsSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type rhs(rhsSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< bool >::type print(printSEXP);
    rcpp_result_gen = Rcpp::wrap(solveBiCGSTAB(lhs, rhs, maxit, tol, print));
    return rcpp_result_gen;
END_RCPP
}
// solveJacobi
Eigen::VectorXd solveJacobi(Eigen::SparseMatrix<double>& L, Eigen::VectorXd& rhs, Eigen::VectorXd guess, double tol);
RcppExport SEXP _implicit_solveJacobi(SEXP LSEXP, SEXP rhsSEXP, SEXP guessSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double>& >::type L(LSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type rhs(rhsSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type guess(guessSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(solveJacobi(L, rhs, guess, tol));
    return rcpp_result_gen;
END_RCPP
}
// makeLaplacian
void makeLaplacian(Eigen::SparseMatrix<double>& L, Eigen::VectorXd mu);
RcppExport SEXP _implicit_makeLaplacian(SEXP LSEXP, SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double>& >::type L(LSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type mu(muSEXP);
    makeLaplacian(L, mu);
    return R_NilValue;
END_RCPP
}
// make_I_Laplacian
void make_I_Laplacian(Eigen::SparseMatrix<double>& I_L, Eigen::VectorXd mu);
RcppExport SEXP _implicit_make_I_Laplacian(SEXP I_LSEXP, SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double>& >::type I_L(I_LSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type mu(muSEXP);
    make_I_Laplacian(I_L, mu);
    return R_NilValue;
END_RCPP
}
// lapTest
std::vector<Eigen::SparseMatrix<double> > lapTest(Eigen::VectorXd params, Eigen::MatrixXd X, Eigen::VectorXd init, Eigen::SparseMatrix<double> neighborMatrix, int nSpace, int nTime, bool KF, double tol, bool useExplicit);
RcppExport SEXP _implicit_lapTest(SEXP paramsSEXP, SEXP XSEXP, SEXP initSEXP, SEXP neighborMatrixSEXP, SEXP nSpaceSEXP, SEXP nTimeSEXP, SEXP KFSEXP, SEXP tolSEXP, SEXP useExplicitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type init(initSEXP);
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type neighborMatrix(neighborMatrixSEXP);
    Rcpp::traits::input_parameter< int >::type nSpace(nSpaceSEXP);
    Rcpp::traits::input_parameter< int >::type nTime(nTimeSEXP);
    Rcpp::traits::input_parameter< bool >::type KF(KFSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< bool >::type useExplicit(useExplicitSEXP);
    rcpp_result_gen = Rcpp::wrap(lapTest(params, X, init, neighborMatrix, nSpace, nTime, KF, tol, useExplicit));
    return rcpp_result_gen;
END_RCPP
}
// computeDiffusion
Eigen::MatrixXd computeDiffusion(Eigen::VectorXd params, Eigen::MatrixXd X, Eigen::VectorXd init, Eigen::SparseMatrix<double> neighborMatrix, int nSpace, int nTime, bool KF, double tol, bool useExplicit);
RcppExport SEXP _implicit_computeDiffusion(SEXP paramsSEXP, SEXP XSEXP, SEXP initSEXP, SEXP neighborMatrixSEXP, SEXP nSpaceSEXP, SEXP nTimeSEXP, SEXP KFSEXP, SEXP tolSEXP, SEXP useExplicitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type init(initSEXP);
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type neighborMatrix(neighborMatrixSEXP);
    Rcpp::traits::input_parameter< int >::type nSpace(nSpaceSEXP);
    Rcpp::traits::input_parameter< int >::type nTime(nTimeSEXP);
    Rcpp::traits::input_parameter< bool >::type KF(KFSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< bool >::type useExplicit(useExplicitSEXP);
    rcpp_result_gen = Rcpp::wrap(computeDiffusion(params, X, init, neighborMatrix, nSpace, nTime, KF, tol, useExplicit));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_implicit_solveBiCGSTAB", (DL_FUNC) &_implicit_solveBiCGSTAB, 5},
    {"_implicit_solveJacobi", (DL_FUNC) &_implicit_solveJacobi, 4},
    {"_implicit_makeLaplacian", (DL_FUNC) &_implicit_makeLaplacian, 2},
    {"_implicit_make_I_Laplacian", (DL_FUNC) &_implicit_make_I_Laplacian, 2},
    {"_implicit_lapTest", (DL_FUNC) &_implicit_lapTest, 9},
    {"_implicit_computeDiffusion", (DL_FUNC) &_implicit_computeDiffusion, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_implicit(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
