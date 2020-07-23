// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// dst1
Eigen::MatrixXd dst1(Eigen::MatrixXd U);
RcppExport SEXP _implicit_dst1(SEXP USEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type U(USEXP);
    rcpp_result_gen = Rcpp::wrap(dst1(U));
    return rcpp_result_gen;
END_RCPP
}
// homogeneousDiffusion_derivatives
Rcpp::List homogeneousDiffusion_derivatives(double mu_0, Eigen::VectorXd gamma, Eigen::VectorXd longLat, double sigma, double kappa, Eigen::MatrixXd coords, Eigen::MatrixXd X_reaction, Eigen::VectorXi positive, Eigen::VectorXi cell, Eigen::VectorXi time, int rows, int cols, int nTime, int diffusionType, double tol, int nIter, double lengthX, double lengthY, bool pad);
RcppExport SEXP _implicit_homogeneousDiffusion_derivatives(SEXP mu_0SEXP, SEXP gammaSEXP, SEXP longLatSEXP, SEXP sigmaSEXP, SEXP kappaSEXP, SEXP coordsSEXP, SEXP X_reactionSEXP, SEXP positiveSEXP, SEXP cellSEXP, SEXP timeSEXP, SEXP rowsSEXP, SEXP colsSEXP, SEXP nTimeSEXP, SEXP diffusionTypeSEXP, SEXP tolSEXP, SEXP nIterSEXP, SEXP lengthXSEXP, SEXP lengthYSEXP, SEXP padSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mu_0(mu_0SEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type longLat(longLatSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type coords(coordsSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X_reaction(X_reactionSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXi >::type positive(positiveSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXi >::type cell(cellSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXi >::type time(timeSEXP);
    Rcpp::traits::input_parameter< int >::type rows(rowsSEXP);
    Rcpp::traits::input_parameter< int >::type cols(colsSEXP);
    Rcpp::traits::input_parameter< int >::type nTime(nTimeSEXP);
    Rcpp::traits::input_parameter< int >::type diffusionType(diffusionTypeSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type nIter(nIterSEXP);
    Rcpp::traits::input_parameter< double >::type lengthX(lengthXSEXP);
    Rcpp::traits::input_parameter< double >::type lengthY(lengthYSEXP);
    Rcpp::traits::input_parameter< bool >::type pad(padSEXP);
    rcpp_result_gen = Rcpp::wrap(homogeneousDiffusion_derivatives(mu_0, gamma, longLat, sigma, kappa, coords, X_reaction, positive, cell, time, rows, cols, nTime, diffusionType, tol, nIter, lengthX, lengthY, pad));
    return rcpp_result_gen;
END_RCPP
}
// computeDiffusion2
Eigen::MatrixXd computeDiffusion2(Eigen::VectorXd params, Eigen::MatrixXd X, Eigen::VectorXd init, int rows, int cols, int nTime, int diffusionType, double tol, int nIter, double lengthX, double lengthY, bool pad);
RcppExport SEXP _implicit_computeDiffusion2(SEXP paramsSEXP, SEXP XSEXP, SEXP initSEXP, SEXP rowsSEXP, SEXP colsSEXP, SEXP nTimeSEXP, SEXP diffusionTypeSEXP, SEXP tolSEXP, SEXP nIterSEXP, SEXP lengthXSEXP, SEXP lengthYSEXP, SEXP padSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type init(initSEXP);
    Rcpp::traits::input_parameter< int >::type rows(rowsSEXP);
    Rcpp::traits::input_parameter< int >::type cols(colsSEXP);
    Rcpp::traits::input_parameter< int >::type nTime(nTimeSEXP);
    Rcpp::traits::input_parameter< int >::type diffusionType(diffusionTypeSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type nIter(nIterSEXP);
    Rcpp::traits::input_parameter< double >::type lengthX(lengthXSEXP);
    Rcpp::traits::input_parameter< double >::type lengthY(lengthYSEXP);
    Rcpp::traits::input_parameter< bool >::type pad(padSEXP);
    rcpp_result_gen = Rcpp::wrap(computeDiffusion2(params, X, init, rows, cols, nTime, diffusionType, tol, nIter, lengthX, lengthY, pad));
    return rcpp_result_gen;
END_RCPP
}
// computeDiffusion_general
Eigen::MatrixXd computeDiffusion_general(Eigen::VectorXd params, Eigen::MatrixXd X, Eigen::VectorXd init, int rows, int cols, int nTime, double tol, int nIter, Eigen::VectorXi internalPoints, bool dirichlet, double lengthX, double lengthY, bool pad, bool iterative, int diffusionType);
RcppExport SEXP _implicit_computeDiffusion_general(SEXP paramsSEXP, SEXP XSEXP, SEXP initSEXP, SEXP rowsSEXP, SEXP colsSEXP, SEXP nTimeSEXP, SEXP tolSEXP, SEXP nIterSEXP, SEXP internalPointsSEXP, SEXP dirichletSEXP, SEXP lengthXSEXP, SEXP lengthYSEXP, SEXP padSEXP, SEXP iterativeSEXP, SEXP diffusionTypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type init(initSEXP);
    Rcpp::traits::input_parameter< int >::type rows(rowsSEXP);
    Rcpp::traits::input_parameter< int >::type cols(colsSEXP);
    Rcpp::traits::input_parameter< int >::type nTime(nTimeSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type nIter(nIterSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXi >::type internalPoints(internalPointsSEXP);
    Rcpp::traits::input_parameter< bool >::type dirichlet(dirichletSEXP);
    Rcpp::traits::input_parameter< double >::type lengthX(lengthXSEXP);
    Rcpp::traits::input_parameter< double >::type lengthY(lengthYSEXP);
    Rcpp::traits::input_parameter< bool >::type pad(padSEXP);
    Rcpp::traits::input_parameter< bool >::type iterative(iterativeSEXP);
    Rcpp::traits::input_parameter< int >::type diffusionType(diffusionTypeSEXP);
    rcpp_result_gen = Rcpp::wrap(computeDiffusion_general(params, X, init, rows, cols, nTime, tol, nIter, internalPoints, dirichlet, lengthX, lengthY, pad, iterative, diffusionType));
    return rcpp_result_gen;
END_RCPP
}
// invert
Eigen::VectorXd invert(int rows, int cols, Eigen::VectorXd& mu, Eigen::VectorXd& x, Eigen::VectorXd& b, int diffusionType, int nIter, double lengthX, double lengthY, int preconditionerType, bool debug);
RcppExport SEXP _implicit_invert(SEXP rowsSEXP, SEXP colsSEXP, SEXP muSEXP, SEXP xSEXP, SEXP bSEXP, SEXP diffusionTypeSEXP, SEXP nIterSEXP, SEXP lengthXSEXP, SEXP lengthYSEXP, SEXP preconditionerTypeSEXP, SEXP debugSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type rows(rowsSEXP);
    Rcpp::traits::input_parameter< int >::type cols(colsSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type x(xSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type diffusionType(diffusionTypeSEXP);
    Rcpp::traits::input_parameter< int >::type nIter(nIterSEXP);
    Rcpp::traits::input_parameter< double >::type lengthX(lengthXSEXP);
    Rcpp::traits::input_parameter< double >::type lengthY(lengthYSEXP);
    Rcpp::traits::input_parameter< int >::type preconditionerType(preconditionerTypeSEXP);
    Rcpp::traits::input_parameter< bool >::type debug(debugSEXP);
    rcpp_result_gen = Rcpp::wrap(invert(rows, cols, mu, x, b, diffusionType, nIter, lengthX, lengthY, preconditionerType, debug));
    return rcpp_result_gen;
END_RCPP
}
// invert_Irregular
Eigen::VectorXd invert_Irregular(int rows, int cols, Eigen::VectorXd& mu, Eigen::VectorXd& x, Eigen::VectorXd& b, int diffusionType, double tol, int nIter, Eigen::VectorXi& internalPoints, bool dirichlet, double lengthX, double lengthY, int preconditionerType, bool debug);
RcppExport SEXP _implicit_invert_Irregular(SEXP rowsSEXP, SEXP colsSEXP, SEXP muSEXP, SEXP xSEXP, SEXP bSEXP, SEXP diffusionTypeSEXP, SEXP tolSEXP, SEXP nIterSEXP, SEXP internalPointsSEXP, SEXP dirichletSEXP, SEXP lengthXSEXP, SEXP lengthYSEXP, SEXP preconditionerTypeSEXP, SEXP debugSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type rows(rowsSEXP);
    Rcpp::traits::input_parameter< int >::type cols(colsSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type x(xSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type diffusionType(diffusionTypeSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type nIter(nIterSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXi& >::type internalPoints(internalPointsSEXP);
    Rcpp::traits::input_parameter< bool >::type dirichlet(dirichletSEXP);
    Rcpp::traits::input_parameter< double >::type lengthX(lengthXSEXP);
    Rcpp::traits::input_parameter< double >::type lengthY(lengthYSEXP);
    Rcpp::traits::input_parameter< int >::type preconditionerType(preconditionerTypeSEXP);
    Rcpp::traits::input_parameter< bool >::type debug(debugSEXP);
    rcpp_result_gen = Rcpp::wrap(invert_Irregular(rows, cols, mu, x, b, diffusionType, tol, nIter, internalPoints, dirichlet, lengthX, lengthY, preconditionerType, debug));
    return rcpp_result_gen;
END_RCPP
}
// logLikelihood_adjusted
double logLikelihood_adjusted(Eigen::VectorXi& positive, Eigen::VectorXi& cell, Eigen::VectorXi& time, Eigen::VectorXd& age, Eigen::VectorXd& sex, Eigen::MatrixXd& u, Eigen::VectorXd& eta);
RcppExport SEXP _implicit_logLikelihood_adjusted(SEXP positiveSEXP, SEXP cellSEXP, SEXP timeSEXP, SEXP ageSEXP, SEXP sexSEXP, SEXP uSEXP, SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXi& >::type positive(positiveSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXi& >::type cell(cellSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXi& >::type time(timeSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type age(ageSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type sex(sexSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type u(uSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(logLikelihood_adjusted(positive, cell, time, age, sex, u, eta));
    return rcpp_result_gen;
END_RCPP
}
// logLikelihood_raw
double logLikelihood_raw(Eigen::VectorXi& positive, Eigen::VectorXi& cell, Eigen::VectorXi& time, Eigen::MatrixXd& u);
RcppExport SEXP _implicit_logLikelihood_raw(SEXP positiveSEXP, SEXP cellSEXP, SEXP timeSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXi& >::type positive(positiveSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXi& >::type cell(cellSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXi& >::type time(timeSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(logLikelihood_raw(positive, cell, time, u));
    return rcpp_result_gen;
END_RCPP
}
// derivative_logLikelihood_adjusted
double derivative_logLikelihood_adjusted(Eigen::VectorXi& positive, Eigen::VectorXi& cell, Eigen::VectorXi& time, Eigen::VectorXd& age, Eigen::VectorXd& sex, Eigen::MatrixXd& u, Eigen::MatrixXd& du_dTheta, Eigen::VectorXd& eta);
RcppExport SEXP _implicit_derivative_logLikelihood_adjusted(SEXP positiveSEXP, SEXP cellSEXP, SEXP timeSEXP, SEXP ageSEXP, SEXP sexSEXP, SEXP uSEXP, SEXP du_dThetaSEXP, SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXi& >::type positive(positiveSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXi& >::type cell(cellSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXi& >::type time(timeSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type age(ageSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type sex(sexSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type u(uSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type du_dTheta(du_dThetaSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(derivative_logLikelihood_adjusted(positive, cell, time, age, sex, u, du_dTheta, eta));
    return rcpp_result_gen;
END_RCPP
}
// etaDerivative_logLikelihood_adjusted
Eigen::VectorXd etaDerivative_logLikelihood_adjusted(Eigen::VectorXi& positive, Eigen::VectorXi& cell, Eigen::VectorXi& time, Eigen::VectorXd& age, Eigen::VectorXd& sex, Eigen::MatrixXd& u, Eigen::VectorXd& eta);
RcppExport SEXP _implicit_etaDerivative_logLikelihood_adjusted(SEXP positiveSEXP, SEXP cellSEXP, SEXP timeSEXP, SEXP ageSEXP, SEXP sexSEXP, SEXP uSEXP, SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXi& >::type positive(positiveSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXi& >::type cell(cellSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXi& >::type time(timeSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type age(ageSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type sex(sexSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type u(uSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(etaDerivative_logLikelihood_adjusted(positive, cell, time, age, sex, u, eta));
    return rcpp_result_gen;
END_RCPP
}
// derivative_logLikelihood_raw
double derivative_logLikelihood_raw(Eigen::VectorXi& positive, Eigen::VectorXi& cell, Eigen::VectorXi& time, Eigen::MatrixXd& u, Eigen::MatrixXd& du_dTheta);
RcppExport SEXP _implicit_derivative_logLikelihood_raw(SEXP positiveSEXP, SEXP cellSEXP, SEXP timeSEXP, SEXP uSEXP, SEXP du_dThetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXi& >::type positive(positiveSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXi& >::type cell(cellSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXi& >::type time(timeSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type u(uSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type du_dTheta(du_dThetaSEXP);
    rcpp_result_gen = Rcpp::wrap(derivative_logLikelihood_raw(positive, cell, time, u, du_dTheta));
    return rcpp_result_gen;
END_RCPP
}
// logLikelihood_general
double logLikelihood_general(Eigen::MatrixXd& negatives, Eigen::MatrixXd& positives, Eigen::MatrixXd& probabilities, Eigen::VectorXi& internalPoints, int rows, int cols);
RcppExport SEXP _implicit_logLikelihood_general(SEXP negativesSEXP, SEXP positivesSEXP, SEXP probabilitiesSEXP, SEXP internalPointsSEXP, SEXP rowsSEXP, SEXP colsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type negatives(negativesSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type positives(positivesSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type probabilities(probabilitiesSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXi& >::type internalPoints(internalPointsSEXP);
    Rcpp::traits::input_parameter< int >::type rows(rowsSEXP);
    Rcpp::traits::input_parameter< int >::type cols(colsSEXP);
    rcpp_result_gen = Rcpp::wrap(logLikelihood_general(negatives, positives, probabilities, internalPoints, rows, cols));
    return rcpp_result_gen;
END_RCPP
}
// fick_L_f
Eigen::VectorXd fick_L_f(Eigen::VectorXd& mu_, Eigen::VectorXd& f_, int rows, int cols, double lengthX, double lengthY);
RcppExport SEXP _implicit_fick_L_f(SEXP mu_SEXP, SEXP f_SEXP, SEXP rowsSEXP, SEXP colsSEXP, SEXP lengthXSEXP, SEXP lengthYSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type mu_(mu_SEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type f_(f_SEXP);
    Rcpp::traits::input_parameter< int >::type rows(rowsSEXP);
    Rcpp::traits::input_parameter< int >::type cols(colsSEXP);
    Rcpp::traits::input_parameter< double >::type lengthX(lengthXSEXP);
    Rcpp::traits::input_parameter< double >::type lengthY(lengthYSEXP);
    rcpp_result_gen = Rcpp::wrap(fick_L_f(mu_, f_, rows, cols, lengthX, lengthY));
    return rcpp_result_gen;
END_RCPP
}
// homogeneous_L_f
Eigen::VectorXd homogeneous_L_f(Eigen::VectorXd& f_, int rows, int cols, double lengthX, double lengthY);
RcppExport SEXP _implicit_homogeneous_L_f(SEXP f_SEXP, SEXP rowsSEXP, SEXP colsSEXP, SEXP lengthXSEXP, SEXP lengthYSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type f_(f_SEXP);
    Rcpp::traits::input_parameter< int >::type rows(rowsSEXP);
    Rcpp::traits::input_parameter< int >::type cols(colsSEXP);
    Rcpp::traits::input_parameter< double >::type lengthX(lengthXSEXP);
    Rcpp::traits::input_parameter< double >::type lengthY(lengthYSEXP);
    rcpp_result_gen = Rcpp::wrap(homogeneous_L_f(f_, rows, cols, lengthX, lengthY));
    return rcpp_result_gen;
END_RCPP
}
// general_Homogeneous_Lf
Eigen::VectorXd general_Homogeneous_Lf(Eigen::VectorXd& f_, int rows, int cols, Eigen::VectorXi& boundaryPoints_, bool dirichlet, double lengthX, double lengthY);
RcppExport SEXP _implicit_general_Homogeneous_Lf(SEXP f_SEXP, SEXP rowsSEXP, SEXP colsSEXP, SEXP boundaryPoints_SEXP, SEXP dirichletSEXP, SEXP lengthXSEXP, SEXP lengthYSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type f_(f_SEXP);
    Rcpp::traits::input_parameter< int >::type rows(rowsSEXP);
    Rcpp::traits::input_parameter< int >::type cols(colsSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXi& >::type boundaryPoints_(boundaryPoints_SEXP);
    Rcpp::traits::input_parameter< bool >::type dirichlet(dirichletSEXP);
    Rcpp::traits::input_parameter< double >::type lengthX(lengthXSEXP);
    Rcpp::traits::input_parameter< double >::type lengthY(lengthYSEXP);
    rcpp_result_gen = Rcpp::wrap(general_Homogeneous_Lf(f_, rows, cols, boundaryPoints_, dirichlet, lengthX, lengthY));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_implicit_dst1", (DL_FUNC) &_implicit_dst1, 1},
    {"_implicit_homogeneousDiffusion_derivatives", (DL_FUNC) &_implicit_homogeneousDiffusion_derivatives, 19},
    {"_implicit_computeDiffusion2", (DL_FUNC) &_implicit_computeDiffusion2, 12},
    {"_implicit_computeDiffusion_general", (DL_FUNC) &_implicit_computeDiffusion_general, 15},
    {"_implicit_invert", (DL_FUNC) &_implicit_invert, 11},
    {"_implicit_invert_Irregular", (DL_FUNC) &_implicit_invert_Irregular, 14},
    {"_implicit_logLikelihood_adjusted", (DL_FUNC) &_implicit_logLikelihood_adjusted, 7},
    {"_implicit_logLikelihood_raw", (DL_FUNC) &_implicit_logLikelihood_raw, 4},
    {"_implicit_derivative_logLikelihood_adjusted", (DL_FUNC) &_implicit_derivative_logLikelihood_adjusted, 8},
    {"_implicit_etaDerivative_logLikelihood_adjusted", (DL_FUNC) &_implicit_etaDerivative_logLikelihood_adjusted, 7},
    {"_implicit_derivative_logLikelihood_raw", (DL_FUNC) &_implicit_derivative_logLikelihood_raw, 5},
    {"_implicit_logLikelihood_general", (DL_FUNC) &_implicit_logLikelihood_general, 6},
    {"_implicit_fick_L_f", (DL_FUNC) &_implicit_fick_L_f, 6},
    {"_implicit_homogeneous_L_f", (DL_FUNC) &_implicit_homogeneous_L_f, 5},
    {"_implicit_general_Homogeneous_Lf", (DL_FUNC) &_implicit_general_Homogeneous_Lf, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_implicit(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
