#include <Rcpp.h>
using namespace Rcpp;

#include <RcppEigen.h>
#include <math.h>
#include "convenienceFunctions.h"

//[[Rcpp::export]]
double logLikelihood_adjusted(Eigen::VectorXi& positive,
                              Eigen::VectorXi& cell,
                              Eigen::VectorXi& time,
                              Eigen::VectorXd& age,
                              Eigen::VectorXd& sex,
                              Eigen::MatrixXd& u,
                              Eigen::VectorXd& eta){
  double logLike=0;
  for (int i=0;i<positive.size();i++){
    int time_i=time(i);
    int cell_i=cell(i);
    int pos_i=positive(i);
    int sex_i=sex(i);
    int age_i=age(i);
    double probRaw_i=u(cell_i,time_i);
    double logit_i=logit(probRaw_i)+eta(0)*age(i)+eta(1)*sex(i);
    double probAdjusted_i=expit(logit_i);
    logLike+=pos_i*std::log(probAdjusted_i)+(1-pos_i)*std::log(1-probAdjusted_i);
  }
  return(logLike);
}

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

//[[Rcpp::export]]
Eigen::VectorXd etaDerivative_logLikelihood_adjusted(Eigen::VectorXi& positive,
                                         Eigen::VectorXi& cell,
                                         Eigen::VectorXi& time,
                                         Eigen::VectorXd& age,
                                         Eigen::VectorXd& sex,
                                         Eigen::MatrixXd& u,
                                         Eigen::VectorXd& eta){
  Eigen::VectorXd derivative;
  derivative.setZero(2);
  for (int i=0;i<positive.size();i++){
    int time_i=time(i);
    int cell_i=cell(i);
    int pos_i=positive(i);
    int sex_i=sex(i);
    int age_i=age(i);
    double probRaw_i=u(cell_i,time_i);
    double logit_i=logit(probRaw_i)+eta(0)*age(i)+eta(1)*sex(i);
    double probAdjusted_i=expit(logit_i);
    derivative(0)+=(pos_i-probAdjusted_i)*age(i);
    derivative(1)+=(pos_i-probAdjusted_i)*sex(i);
  }
  return(derivative);
}


//[[Rcpp::export]]
double derivative_logLikelihood_raw(Eigen::VectorXi& positive,
                                    Eigen::VectorXi& cell,
                                    Eigen::VectorXi& time,
                                    Eigen::MatrixXd& u,
                                    Eigen::MatrixXd& du_dTheta){
  double derivative=0;
  for (int i=0;i<positive.size();i++){
    int time_i=time(i);
    int cell_i=cell(i);
    int pos_i=positive(i);
    double probRaw_i=u(cell_i,time_i);
    double rawVar=probRaw_i*(1-probRaw_i);
    double derivRaw_i=du_dTheta(cell_i,time_i);
    derivative+=(pos_i-probRaw_i)/rawVar*derivRaw_i;
  }
  return(derivative);
}



//[[Rcpp::export]]
double logLikelihood_general(Eigen::MatrixXd& negatives,
                             Eigen::MatrixXd& positives,
                             Eigen::MatrixXd& probabilities,
                             Eigen::VectorXi& internalPoints,
                             int rows,
                             int cols){
  int nTime=negatives.cols();
  Eigen::Map<Eigen::MatrixXi> internal(internalPoints.data(),rows,cols);
  
  double logLike=0;
  
  for (int i=0;i<nTime;i++){
    Eigen::Map<Eigen::MatrixXd> probs_i(probabilities.col(i).data(),rows+2,cols+2);
    Eigen::Map<Eigen::MatrixXd> pos_i(positives.col(i).data(),rows+2,cols+2);
    Eigen::Map<Eigen::MatrixXd> neg_i(negatives.col(i).data(),rows+2,cols+2);
    for (int colInd=0;colInd<cols;colInd++){
      for (int rowInd=0;rowInd<rows;rowInd++){
        if (internal(rowInd,colInd)!=0){
          logLike+=pos_i(rowInd+1,colInd+1)*std::log(probs_i(rowInd+1,colInd+1));
          logLike+=neg_i(rowInd+1,colInd+1)*std::log(1-probs_i(rowInd+1,colInd+1));
        }
      }
    }
  }
  return logLike;
}