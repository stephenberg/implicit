#include <RcppEigen.h>
#include <math.h>
#include "inversionFunctions.h"

//#include "C:/Users/saber/OneDrive/Documents/R/win-library/3.6/RcppEigen/include/unsupported/Eigen/src/AutoDiff/AutoDiffScalar.h"
//#include "C:/Users/saber/OneDrive/Documents/R/win-library/3.6/RcppEigen/include/unsupported/Eigen/src/AutoDiff/AutoDiffVector.h"
//#include "C:/Users/saber/OneDrive/Documents/R/win-library/3.6/RcppEigen/include/unsupported/Eigen/src/AutoDiff/AutoDiffJacobian.h"


using namespace Rcpp;
// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

//[[Rcpp::export]]
double logLikelihood(Eigen::MatrixXd& negatives,
                     Eigen::MatrixXd& positives,
                     Eigen::MatrixXd& probabilities){
  int nSpace=negatives.rows();
  int nTime=negatives.cols();
  double logLike=0;
  for (int i=0;i<nTime;i++){
    for (int j=0;j<nSpace;j++){
      if (negatives(j,i)>0){
        logLike+=(negatives(j,i)*std::log(1-probabilities(j,i)));
      }
      if (positives(j,i)>0){
        logLike+=(positives(j,i)*std::log(probabilities(j,i)));
      }
    }
  }
  return logLike;
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


//[[Rcpp::export]]
Eigen::MatrixXd computeDiffusion2(Eigen::VectorXd params,
                                  Eigen::MatrixXd X,
                                  Eigen::VectorXd init,
                                  int rows,
                                  int cols,
                                  int nTime,
                                  int diffusionType,
                                  double tol,
                                  int nIter,
                                  double lengthX=1,
                                  double lengthY=1,
                                  bool pad=true){
  
  int p=X.cols();
  
  Eigen::VectorXd alpha,gamma;
  alpha.setZero(p);
  gamma.setZero(p);
  
  alpha=params.segment(0,p);
  gamma=params.segment(p,p);
  
  
  //number of internal points, full grid points
  int nInternal=rows*cols;
  int nFull=(rows+2)*(cols+2);
  
  Eigen::MatrixXd U;
  
  U.setZero(nInternal,nTime);
  Eigen::Map<Eigen::MatrixXd> initMap(init.data(), rows + 2, cols + 2);
  Eigen::Map<Eigen::MatrixXd> U0_map(U.col(0).data(), rows, cols);
  U0_map=initMap.block(1,1,rows,cols);
  
  Eigen::VectorXd muVec, lambdaVec;
  muVec.setZero(X.rows());
  lambdaVec.setZero(X.rows());
  
  muVec=(X*alpha).array().exp();
  lambdaVec = (X*gamma).array().exp();
  
  
  Eigen::Map<Eigen::MatrixXd> mu(muVec.data(),nFull,nTime);
  Eigen::Map<Eigen::MatrixXd> lambda(lambdaVec.data(),nFull,nTime);
  
  Eigen::VectorXd rhs(nInternal);
  
  for (int i=1;i<nTime;i++){
    //Kolmogorov-Fisher
    //U2-U1=HU2+Lambda*U1(1-U1)
    //(I-H)U2=U1+Lambda*U1(1-U1)
    rhs=U.col(i-1).array()+lambda.col(i-1).array()*U.col(i-1).array()*(1-U.col(i-1).array());
    Eigen::VectorXd guess(U.col(i-1));
    Eigen::VectorXd mu_i(mu.col(i-1));
    U.col(i)=invert(rows,
          cols,
          mu_i,
          guess,
          rhs,
          diffusionType,
          nIter,
          lengthX,
          lengthY);
  }
  
  if (pad) {
    Eigen::MatrixXd U_full;
    U_full.setZero(nFull, nTime);
    for (int i = 0; i < nTime; i++) {
      Eigen::Map<Eigen::MatrixXd> U_i(U.col(i).data(), rows, cols);
      Eigen::Map<Eigen::MatrixXd> U_full_i(U_full.col(i).data(), rows + 2, cols + 2);
      U_full_i.block(1, 1, rows, cols) = U_i;
    }
    return U_full;
  }
  return U;
}

//[[Rcpp::export]]
Eigen::MatrixXd computeDiffusion_general(Eigen::VectorXd params,
                                         Eigen::MatrixXd X,
                                         Eigen::VectorXd init,
                                         int rows,
                                         int cols,
                                         int nTime,
                                         double tol,
                                         int nIter,
                                         Eigen::VectorXi internalPoints,
                                         bool dirichlet=true,
                                         double lengthX=1,
                                         double lengthY=1,
                                         bool pad=true,
                                         bool iterative=true,
                                         int diffusionType=1){
  int p=X.cols();
  
  Eigen::VectorXd alpha,gamma;
  alpha.setZero(p);
  gamma.setZero(p);
  
  alpha=params.segment(0,p);
  gamma=params.segment(p,p);
  
  
  //number of internal points, full grid points
  int nInternal=rows*cols;
  int nFull=(rows+2)*(cols+2);
  
  Eigen::MatrixXd U;
  
  U.setZero(nInternal,nTime);
  Eigen::Map<Eigen::MatrixXd> initMap(init.data(), rows + 2, cols + 2);
  Eigen::Map<Eigen::MatrixXd> U0_map(U.col(0).data(), rows, cols);
  
  U0_map=initMap.block(1,1,rows,cols);
  U.col(0)=fixBoundary(U.col(0),internalPoints);
  
  Eigen::VectorXd muVec, lambdaVec;
  muVec.setZero(X.rows());
  lambdaVec.setZero(X.rows());
  
  muVec=(X*alpha).array().exp();
  lambdaVec = (X*gamma).array().exp();
  
  
  Eigen::Map<Eigen::MatrixXd> mu(muVec.data(),nFull,nTime);
  Eigen::Map<Eigen::MatrixXd> lambda(lambdaVec.data(),nFull,nTime);
  
  Eigen::VectorXd rhs(nInternal);
  
  for (int i=1;i<nTime;i++){
    //Kolmogorov-Fisher
    //U2-U1=HU2+Lambda*U1(1-U1)
    //(I-H)U2=U1+Lambda*U1(1-U1)
    rhs=U.col(i-1).array()+lambda.col(i-1).array()*U.col(i-1).array()*(1-U.col(i-1).array());
    Eigen::VectorXd guess(U.col(i-1));
    Eigen::VectorXd mu_i(mu.col(i-1));
    U.col(i)=invert_Irregular(rows,
          cols,
          mu_i,
          guess,
          rhs,
          diffusionType,
          tol,
          nIter,
          internalPoints,
          dirichlet,
          lengthX,
          lengthY);
  }
  
  if (pad) {
    Eigen::MatrixXd U_full;
    U_full.setZero(nFull, nTime);
    for (int i = 0; i < nTime; i++) {
      Eigen::Map<Eigen::MatrixXd> U_i(U.col(i).data(), rows, cols);
      Eigen::Map<Eigen::MatrixXd> U_full_i(U_full.col(i).data(), rows + 2, cols + 2);
      U_full_i.block(1, 1, rows, cols) = U_i;
    }
    return U_full;
  }
  return U;
}