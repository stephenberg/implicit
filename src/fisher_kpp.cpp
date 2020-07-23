#include <RcppEigen.h>
#include <math.h>
#include "inversionFunctions.h"
#include "convenienceFunctions.h"
#include "loglikelihoodFunctions.h"

//#include "C:/Users/saber/OneDrive/Documents/R/win-library/3.6/RcppEigen/include/unsupported/Eigen/src/AutoDiff/AutoDiffScalar.h"
//#include "C:/Users/saber/OneDrive/Documents/R/win-library/3.6/RcppEigen/include/unsupported/Eigen/src/AutoDiff/AutoDiffVector.h"
//#include "C:/Users/saber/OneDrive/Documents/R/win-library/3.6/RcppEigen/include/unsupported/Eigen/src/AutoDiff/AutoDiffJacobian.h"


using namespace Rcpp;
// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]




//[[Rcpp::export]]
Rcpp::List homogeneousDiffusion_derivatives(double mu_0,
                                 Eigen::VectorXd gamma,
                                 Eigen::VectorXd longLat,
                                 double sigma,
                                 double kappa,
                                 Eigen::MatrixXd coords,
                                 Eigen::MatrixXd X_reaction,
                                 Eigen::VectorXi positive,
                                 Eigen::VectorXi cell,
                                 Eigen::VectorXi time,
                                 int rows,
                                 int cols,
                                 int nTime,
                                 int diffusionType,
                                 double tol,
                                 int nIter,
                                 double lengthX=1,
                                 double lengthY=1,
                                 bool differentiate = true,
                                 bool debug = false){
  
  //number of internal points, full grid points
  int nInternal=rows*cols;
  int nFull=(rows+2)*(cols+2);
  
  Eigen::MatrixXd U;
  
  U.setZero(nInternal,nTime);
  
  Eigen::MatrixXd init(rows+2,cols+2);
  int count=0;
  for (int cInd=0;cInd<cols+2;cInd++){
    for (int rInd=0;rInd<rows+2;rInd++){
      double d1=coords(count,0)-longLat(0);
      double d2=coords(count,1)-longLat(1);
      double dist2=std::pow(d1,2)+std::pow(d2,2);
      init(rInd,cInd)=expit(kappa)*std::exp(-dist2/(2*std::pow(sigma,2)));
      count=count+1;
    }
  }
  
  Eigen::Map<Eigen::MatrixXd> U0_map(U.col(0).data(), rows, cols);
  U0_map=init.block(1,1,rows,cols);
  
  Eigen::VectorXd muVec, lambdaVec;
  muVec.setZero(X_reaction.rows());
  muVec=muVec.array()+mu_0;
  
  lambdaVec.setZero(X_reaction.rows());
  lambdaVec = (X_reaction*gamma).array().exp();
  
  
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
  
  
  //////derivatives
  std::vector<Eigen::MatrixXd> du_dmu,du_dgamma,du_dlongLat,du_dsigma,du_dkappa;
  Eigen::VectorXd dl_dmu,dl_dgamma,dl_dlongLat,dl_dsigma,dl_dkappa;
  
  if (differentiate){
    
    //wrt mu_0
    Eigen::MatrixXd du_dTheta;
    du_dTheta.setZero(nInternal,nTime);
    dl_dmu.setZero(1);
    
    for (int i=0;i<nTime;i++){
      //(du_0/dTheta)=0
      if (i>0){
        rhs=U.col(i-1).array()+lambda.col(i-1).array()*U.col(i-1).array()*(1-U.col(i-1).array());
        Eigen::VectorXd guess,temp;
        guess.setZero(nInternal);
        temp.setZero(nInternal);
        Eigen::VectorXd mu_i(mu.col(i-1));
        du_dTheta.col(i)=invert(rows,
                      cols,
                      mu_i,
                      guess,
                      rhs,
                      diffusionType,
                      nIter,
                      lengthX,
                      lengthY);
        temp=du_dTheta.col(i);
        du_dTheta.col(i)=homogeneous_L_f(temp,
                      rows,
                      cols,
                      lengthX,
                      lengthY);
        temp=du_dTheta.col(i);
        du_dTheta.col(i)=invert(rows,
                      cols,
                      mu_i,
                      guess,
                      temp,
                      diffusionType,
                      nIter,
                      lengthX,
                      lengthY);
        
        Eigen::VectorXd rhs2;
        rhs2.setZero(nInternal);
        rhs2=du_dTheta.col(i-1).array()+lambda.col(i-1).array()*du_dTheta.col(i-1).array()+
          -2*lambda.col(i-1).array()*U.col(i-1).array()*du_dTheta.col(i-1).array();
          du_dTheta.col(i)=du_dTheta.col(i)+invert(rows,
                        cols,
                        mu_i,
                        guess,
                        rhs2,
                        diffusionType,
                        nIter,
                        lengthX,
                        lengthY);
      }
    }
    
    if (debug) du_dmu.push_back(du_dTheta);
    
    dl_dmu(0)=derivative_logLikelihood_raw(positive,
           cell,
           time,
           U,
           du_dTheta);
    
    //with respect to gamma
    int p=X_reaction.cols();
    dl_dgamma.setZero(p);
    for (int pInd=0;pInd<p;pInd++){
      du_dTheta.setZero(nInternal,nTime);
      for (int i=0;i<nTime;i++){
        //(du_0/dTheta)=0
        if (i>0){
          Eigen::VectorXd rhs2,guess;
          guess.setZero(nInternal);
          
          rhs2.setZero(nInternal);
          rhs2=du_dTheta.col(i-1).array()+
            X_reaction.col(pInd).array()*lambda.col(i-1).array()*U.col(i-1).array()*(1-U.col(i-1).array()).array()+
            lambda.col(i-1).array()*du_dTheta.col(i-1).array()+
            -2*lambda.col(i-1).array()*U.col(i-1).array()*du_dTheta.col(i-1).array();
            Eigen::VectorXd mu_i(mu.col(i-1));
            du_dTheta.col(i)=invert(rows,
                          cols,
                          mu_i,
                          guess,
                          rhs2,
                          diffusionType,
                          nIter,
                          lengthX,
                          lengthY);
        }
      }
      if (debug) du_dgamma.push_back(du_dTheta);
      dl_dgamma(pInd)=derivative_logLikelihood_raw(positive,
                cell,
                time,
                U,
                du_dTheta);
    }
    
    //with respect to longLat
    dl_dlongLat.setZero(2);
    Eigen::MatrixXd initFull(rows+2,cols+2);
    Eigen::Map<Eigen::MatrixXd> du_dTheta0_map(du_dTheta.col(0).data(),rows,cols);
    
    for (int llInd=0;llInd<2;llInd++){
      du_dTheta.setZero(nInternal,nTime);
      
      
      for (int i=0;i<nTime;i++){
        if (i==0){
          int count=0;
          for (int cInd=0;cInd<cols+2;cInd++){
            for (int rInd=0;rInd<rows+2;rInd++){
              Eigen::VectorXd d(2);
              d(0)=coords(count,0)-longLat(0);
              d(1)=coords(count,1)-longLat(1);
              double dist2=std::pow(d(0),2)+std::pow(d(1),2);
              initFull(rInd,cInd)=expit(kappa)*std::exp(-dist2/(2*std::pow(sigma,2)))*d(llInd)/std::pow(sigma,2);
              count=count+1;
            }
          }
          du_dTheta0_map=initFull.block(1,1,rows,cols);
        }
        else{
          Eigen::VectorXd rhs2,guess;
          guess.setZero(nInternal);
          
          rhs2.setZero(nInternal);
          rhs2=du_dTheta.col(i-1).array()+
            lambda.col(i-1).array()*du_dTheta.col(i-1).array()+
            -2*lambda.col(i-1).array()*U.col(i-1).array()*du_dTheta.col(i-1).array();
            Eigen::VectorXd mu_i(mu.col(i-1));
            du_dTheta.col(i)=invert(rows,
                          cols,
                          mu_i,
                          guess,
                          rhs2,
                          diffusionType,
                          nIter,
                          lengthX,
                          lengthY);
        }
      }
      if (debug) du_dlongLat.push_back(du_dTheta);
      dl_dlongLat(llInd)=derivative_logLikelihood_raw(positive,
                  cell,
                  time,
                  U,
                  du_dTheta);
    }
    
    //with respect to sigma
    du_dTheta.setZero(nInternal,nTime);
    dl_dsigma.setZero(1);
    for (int i=0;i<nTime;i++){
      if (i==0){
        int count=0;
        for (int cInd=0;cInd<cols+2;cInd++){
          for (int rInd=0;rInd<rows+2;rInd++){
            double d1=coords(count,0)-longLat(0);
            double d2=coords(count,1)-longLat(1);
            double dist2=std::pow(d1,2)+std::pow(d2,2);
            initFull(rInd,cInd)=expit(kappa)*std::exp(-dist2/(2*std::pow(sigma,2)))*dist2/std::pow(sigma,3);
            count=count+1;
          }
        }
        du_dTheta0_map=initFull.block(1,1,rows,cols);
      }
      else{
        Eigen::VectorXd rhs2,guess;
        guess.setZero(nInternal);
        
        rhs2.setZero(nInternal);
        rhs2=du_dTheta.col(i-1).array()+
          lambda.col(i-1).array()*du_dTheta.col(i-1).array()+
          -2*lambda.col(i-1).array()*U.col(i-1).array()*du_dTheta.col(i-1).array();
          Eigen::VectorXd mu_i(mu.col(i-1));
          du_dTheta.col(i)=invert(rows,
                        cols,
                        mu_i,
                        guess,
                        rhs2,
                        diffusionType,
                        nIter,
                        lengthX,
                        lengthY);
      }
    }
    if (debug) du_dsigma.push_back(du_dTheta);
    dl_dsigma(0)=derivative_logLikelihood_raw(positive,
              cell,
              time,
              U,
              du_dTheta);
    
    //with respect to kappa
    dl_dkappa.setZero(1);
    for (int i=0;i<nTime;i++){
      if (i==0){
        int count=0;
        for (int cInd=0;cInd<cols+2;cInd++){
          for (int rInd=0;rInd<rows+2;rInd++){
            double d1=coords(count,0)-longLat(0);
            double d2=coords(count,1)-longLat(1);
            double dist2=std::pow(d1,2)+std::pow(d2,2);
            initFull(rInd,cInd)=expit(kappa)*(1-expit(kappa))*std::exp(-dist2/(2*std::pow(sigma,2)));
            count=count+1;
          }
        }
        du_dTheta0_map=initFull.block(1,1,rows,cols);
      }
      else{
        Eigen::VectorXd rhs2,guess;
        guess.setZero(nInternal);
        
        rhs2.setZero(nInternal);
        rhs2=du_dTheta.col(i-1).array()+
          lambda.col(i-1).array()*du_dTheta.col(i-1).array()+
          -2*lambda.col(i-1).array()*U.col(i-1).array()*du_dTheta.col(i-1).array();
          Eigen::VectorXd mu_i(mu.col(i-1));
          du_dTheta.col(i)=invert(rows,
                        cols,
                        mu_i,
                        guess,
                        rhs2,
                        diffusionType,
                        nIter,
                        lengthX,
                        lengthY);
      }
    }
    if (debug) du_dkappa.push_back(du_dTheta);
    dl_dkappa(0)=derivative_logLikelihood_raw(positive,
              cell,
              time,
              U,
              du_dTheta);
  }
  

  Eigen::MatrixXd U_full;
  Eigen::MatrixXd du_dTheta_full;
  U_full.setZero(nFull, nTime);
  du_dTheta_full.setZero(nFull, nTime);
  for (int i = 0; i < nTime; i++) {
    Eigen::Map<Eigen::MatrixXd> U_i(U.col(i).data(), rows, cols);
    Eigen::Map<Eigen::MatrixXd> U_full_i(U_full.col(i).data(), rows + 2, cols + 2);
    U_full_i.block(1, 1, rows, cols) = U_i;
  }
  double logLike = logLikelihood_raw(positive,cell,time,U);
  
  
  if (differentiate==false){
    return Rcpp::List::create(Rcpp::Named("u") = U_full,
                              Rcpp::Named("logLike") = logLike);
  }
  
  if (debug){
    return Rcpp::List::create(Rcpp::Named("u") = U_full,
                              Rcpp::Named("logLike") = logLike,
                              Rcpp::Named("dl_dmu") = dl_dmu,
                              Rcpp::Named("dl_dgamma") = dl_dgamma,
                              Rcpp::Named("dl_dlongLat") = dl_dlongLat,
                              Rcpp::Named("dl_dsigma") = dl_dsigma,
                              Rcpp::Named("dl_dkappa") = dl_dkappa,
                              Rcpp::Named("du_dmu") = du_dmu,
                              Rcpp::Named("du_dgamma") = du_dgamma,
                              Rcpp::Named("du_dlongLat") = du_dlongLat,
                              Rcpp::Named("du_dsigma") = du_dsigma,
                              Rcpp::Named("du_dkappa") = du_dkappa);
  }
  return Rcpp::List::create(Rcpp::Named("u") = U_full,
                            Rcpp::Named("logLike") = logLike,
                            Rcpp::Named("dl_dmu") = dl_dmu,
                            Rcpp::Named("dl_dgamma") = dl_dgamma,
                            Rcpp::Named("dl_dlongLat") = dl_dlongLat,
                            Rcpp::Named("dl_dsigma") = dl_dsigma,
                            Rcpp::Named("dl_dkappa") = dl_dkappa);
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