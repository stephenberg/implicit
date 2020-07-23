#include <Rcpp.h>
#include <RcppEigen.h>
#include "fftPlan.h"
#include "multiplyByLaplacian.h"


Eigen::VectorXd fixBoundary(Eigen::VectorXd x,
                            Eigen::VectorXi& internalPoints) {
  for (int i = 0; i < x.size(); i++) {
    x(i) = x(i) * internalPoints(i);
  }
  return(x);
}

//solve (I-Laplacian)x=b for x, where the Laplacian corresponds to the diffusion
//type as described in L_f; the Laplacian corresponds to a rows by cols grid 
//of *internal points* surrounded implicitly by 
//0's (Dirichlet boundary)
//nIter is the number of iterations, the argument x is an initial guess
//[[Rcpp::export]]
Eigen::VectorXd invert(int rows,
                       int cols,
                       Eigen::VectorXd& mu,
                       Eigen::VectorXd& x,
                       Eigen::VectorXd& b,
                       int diffusionType,
                       int nIter,
                       double lengthX=1,
                       double lengthY=1,
                       int preconditionerType=1,
                       bool debug=false){
  
  int N=rows*cols;
  
  Eigen::VectorXd r(N),p(N),y(N),Ap(N);
  
  
  
  Eigen::Map<MatrixXd> muMatrix(mu.data(),rows+2,cols+2);
  Eigen::MatrixXd internalMuMatrix(muMatrix.block(1,1,rows,cols).array().pow(-1.0));
  Eigen::Map<VectorXd> internalMu_inverse(internalMuMatrix.data(),rows*cols);
  double delta=std::pow(internalMu_inverse.array().sum()/(rows*cols),-1);
  if (debug) Rcout<<"delta="<<delta<<"\n";
  
  fftPlan plan(rows,
               cols,
               r.data(),
               Ap.data(),
               y.data());
  if (diffusionType==1){
    r=internalMu_inverse.asDiagonal()*x
    -homogeneous_L_f(x,rows,cols,lengthX,lengthY)-b;
  }
  if (diffusionType==0){
    r=x-fick_L_f(mu,x,rows,cols,lengthX,lengthY)-b;
  }
  if (diffusionType==-1){
    r=internalMu_inverse.asDiagonal()*x-
      homogeneous_L_f(x,rows,cols,lengthX,lengthY)-
      internalMu_inverse.asDiagonal()*b;
  }
  
  
  if (preconditionerType==0){
    if (debug) Rcout<<"Preconditioning by L^{-1}"<<"\n";
    plan.multiplyBy_pow_neg_A(-1.0,lengthX,lengthY);
  }
  else{
    if (debug) Rcout<<"Preconditioning by (I-delta\\cdot L)^{-1}"<<"\n";
    plan.multiplyBy_I_A_inverse(lengthX,lengthY,delta);
  }
  
  p=-y;
  
  double r_norm=std::pow(r.array().pow(2).array().sum(),0.5)/r.size();
  double alphaNum=(r.array()*y.array()).array().sum();
  int count=0;
  while ((r_norm>std::pow(10,-10.0))& (count<nIter)){
    if (diffusionType==1){
      Ap=internalMu_inverse.asDiagonal()*p
      -homogeneous_L_f(p,rows,cols,lengthX,lengthY);
    }
    if (diffusionType==0){
      Ap=p-fick_L_f(mu,p,rows,cols,lengthX,lengthY);
    }
    if (diffusionType==-1){
      Ap=internalMu_inverse.asDiagonal()*p
      -homogeneous_L_f(p,rows,cols,lengthX,lengthY);
    }
    double alpha=alphaNum/((p.adjoint()*Ap)(0));
    x=x+alpha*p;
    r=r+alpha*Ap;
    
    if (preconditionerType==0){
      if (debug) Rcout<<"Preconditioning by L^{-1}"<<"\n";
      plan.multiplyBy_pow_neg_A(-1.0,lengthX,lengthY);
    }
    else{
      if (debug) Rcout<<"Preconditioning by (I-delta\\cdot L)^{-1}"<<"\n";
      plan.multiplyBy_I_A_inverse(lengthX,lengthY,delta);
    }
    
    double beta=(r.array()*y.array()).array().sum()/alphaNum;
    p=-y+beta*p;
    alphaNum=alphaNum*beta;
    r_norm=std::pow(r.array().pow(2).array().sum(),0.5)/r.size();
    count+=1;
    
    if (debug){
      Rcout<<"iter "<<count<<"\n";
      Rcout<<"norm "<<r_norm<<"\n\n";
    }
  }
  if (diffusionType==1){
    return internalMu_inverse.asDiagonal()*x;
  }
  return x;
}


//solve (I-Laplacian)x=b for x, where the Laplacian corresponds to the diffusion
//type as described in L_f; the Laplacian corresponds to a rows by cols grid 
//of *internal points* surrounded implicitly by 
//0's (Dirichlet boundary)
//nIter is the number of iterations, the argument x is an initial guess
//[[Rcpp::export]]
Eigen::VectorXd invert_Irregular(int rows,
                                 int cols,
                                 Eigen::VectorXd& mu,
                                 Eigen::VectorXd& x,
                                 Eigen::VectorXd& b,
                                 int diffusionType,
                                 double tol,
                                 int nIter,
                                 Eigen::VectorXi& internalPoints,
                                 bool dirichlet=true,
                                 double lengthX=1,
                                 double lengthY=1,
                                 int preconditionerType=1,
                                 bool debug=false){
  int N=rows*cols;
  
  Eigen::VectorXd r(N),p(N),y(N),Ap(N);
  
  
  
  Eigen::Map<MatrixXd> muMatrix(mu.data(),rows+2,cols+2);
  Eigen::MatrixXd internalMuMatrix(muMatrix.block(1,1,rows,cols).array().pow(-1.0));
  Eigen::Map<VectorXd> internalMu_inverse(internalMuMatrix.data(),rows*cols);
  
  //make sure we start out with properly zero padded mu, x, b
  //the remaining operations will not affect the zero padding?
  internalMu_inverse = fixBoundary(internalMu_inverse,internalPoints);
  double delta=std::pow(internalMu_inverse.array().sum()/(internalPoints.array().sum()),-1);
  if (debug) Rcout<<"delta="<<delta<<"\n";
  
  x = fixBoundary(x,internalPoints);
  b = fixBoundary(b,internalPoints);
  
  fftPlan plan(rows,
               cols,
               r.data(),
               Ap.data(),
               y.data());
  if (diffusionType==1){
    r=internalMu_inverse.asDiagonal()*x
    -general_Homogeneous_Lf(x,rows,cols,internalPoints,dirichlet,lengthX,lengthY)-b;
  }
  if (diffusionType==0){
    r=x-fixBoundary(fick_L_f(mu,x,rows,cols,lengthX,lengthY),internalPoints)-b;
  }
  if (diffusionType==-1){
    r=internalMu_inverse.asDiagonal()*x-
      general_Homogeneous_Lf(x,rows,cols,internalPoints,dirichlet,lengthX,lengthY)-
      internalMu_inverse.asDiagonal()*b;
  }
  
  
  if (preconditionerType==0){
    if (debug) Rcout<<"Preconditioning by L^{-1}"<<"\n";
    plan.multiplyBy_pow_neg_A(-1.0,lengthX,lengthY);
  }
  else{
    if (debug) Rcout<<"Preconditioning by (I-delta\\cdot L)^{-1}"<<"\n";
    plan.multiplyBy_I_A_inverse(lengthX,lengthY,delta);
  }
  
  p=-y;
  
  double r_norm=std::pow(r.array().pow(2).array().sum(),0.5)/r.size();
  double alphaNum=(r.array()*y.array()).array().sum();
  int count=0;
  while ((r_norm>tol)& (count<nIter)){
    if (diffusionType==1){
      Ap=internalMu_inverse.asDiagonal()*p
      -general_Homogeneous_Lf(p,rows,cols,internalPoints,dirichlet,lengthX,lengthY);
    }
    if (diffusionType==0){
      Ap=p-fixBoundary(fick_L_f(mu,p,rows,cols,lengthX,lengthY),internalPoints);
    }
    if (diffusionType==-1){
      Ap=internalMu_inverse.asDiagonal()*p
      - general_Homogeneous_Lf(p,rows,cols,internalPoints,dirichlet,lengthX,lengthY);
    }
    double alpha=alphaNum/((p.adjoint()*Ap)(0));
    x=x+alpha*p;
    r=r+alpha*Ap;
    
    if (preconditionerType==0){
      if (debug) Rcout<<"Preconditioning by L^{-1}"<<"\n";
      plan.multiplyBy_pow_neg_A(-1.0,lengthX,lengthY);
    }
    else{
      if (debug) Rcout<<"Preconditioning by (I-delta\\cdot L)^{-1}"<<"\n";
      plan.multiplyBy_I_A_inverse(lengthX,lengthY,delta);
    }
    
    double beta=(r.array()*y.array()).array().sum()/alphaNum;
    p=-y+beta*p;
    alphaNum=alphaNum*beta;
    r_norm=std::pow(r.array().pow(2).array().sum(),0.5)/r.size();
    
    count+=1;
    
    if (debug){
      Rcout<<"iter "<<count<<"\n";
      Rcout<<"norm "<<r_norm<<"\n\n";
    }
  }
  
  if (diffusionType==1){
    return internalMu_inverse.asDiagonal()*x;
  }
  return x;
}


