#include <Rcpp.h>
using namespace Rcpp;

#include <RcppEigen.h>
#include <math.h>
#include "convenienceFunctions.h"

#ifndef __LOGLIKE__
#define __LOGLIKE__


//[[Rcpp::export]]
double logLikelihood_raw(Eigen::VectorXi& positive,
                         Eigen::VectorXi& cell,
                         Eigen::VectorXi& time,
                         Eigen::MatrixXd& u){
  double logLike=0;
  for (int i=0;i<positive.size();i++){
    int time_i=time(i);
    int cell_i=cell(i);
    int pos_i=positive(i);
    double probRaw_i=u(cell_i,time_i);
    logLike+=pos_i*std::log(probRaw_i)+(1-pos_i)*std::log(1-probRaw_i);
  }
  return(logLike);
}

//[[Rcpp::export]]
double derivative_logLikelihood_adjusted(Eigen::VectorXi& positive,
                                         Eigen::VectorXi& cell,
                                         Eigen::VectorXi& time,
                                         Eigen::VectorXd& age,
                                         Eigen::VectorXd& sex,
                                         Eigen::MatrixXd& u,
                                         Eigen::MatrixXd& du_dTheta,
                                         Eigen::VectorXd& eta){
  double derivative=0;
  for (int i=0;i<positive.size();i++){
    int time_i=time(i);
    int cell_i=cell(i);
    int pos_i=positive(i);
    int sex_i=sex(i);
    int age_i=age(i);
    double probRaw_i=u(cell_i,time_i);
    double logit_i=logit(probRaw_i)+eta(0)*age(i)+eta(1)*sex(i);
    double probAdjusted_i=expit(logit_i);
    double rawVar=probRaw_i*(1-probRaw_i);
    double derivRaw_i=du_dTheta(cell_i,time_i);
    derivative+=(pos_i-probAdjusted_i)/rawVar*derivRaw_i;
  }
  return(derivative);
}

#endif
