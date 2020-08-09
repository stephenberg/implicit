#ifndef __INVERSION__
#define __INVERSION__


#include <iostream>
//#include <Rcpp.h>
#include <Eigen>
#include "fftPlan.h"
#include "multiplyByLaplacian.h"
#include "grid.h"


//solve (I-Laplacian)x=b for x, where the Laplacian corresponds to the diffusion
//type as described in L_f;
//nIter is the number of iterations, the argument x is an initial guess
Eigen::VectorXd invert(Grid& grid,
                       const Ref<const VectorXd>& mu,
                       VectorXd& x,
                       VectorXd& b,
                       int diffusionType,
                       int nIter,
                       double tol=1.0e-10,
                       int preconditionerType=1,
                       bool dirichlet=true,
                       bool debug=false){
  
  int N = grid.nInternal;
  
  Eigen::VectorXd r(N),p(N),y(N),Ap(N);
  
  

  const Map<const MatrixXd> muMatrix(mu.data(),grid.rows,grid.cols);
  MatrixXd internalMuMatrix(muMatrix.block(1,1,grid.rows_internal,grid.cols_internal).array().pow(-1.0));
  Map<VectorXd> internalMu_inverse(internalMuMatrix.data(),grid.nInternal);
  double delta = std::pow(internalMu_inverse.array().sum() / (grid.internalPoints.array().sum()), -1);

  if (debug) std::cout<<"delta="<<delta<<"\n";
  int temp=0;
  fftPlan plan(grid.rows_internal,
               grid.cols_internal,
               r.data(),
               Ap.data(),
               y.data());
  if (diffusionType==1){
	  if (grid.regularGrid) {
		  r = internalMu_inverse.asDiagonal() * x
			  - homogeneous_Lf(x, grid, dirichlet) - b;
	  }
	  else {
		  r = internalMu_inverse.asDiagonal() * x
			  - general_Homogeneous_Lf(x, grid, dirichlet) - b;
	  }
  }
  if (diffusionType==0){
	  if (grid.regularGrid) {
		  r = x - fick_Lf(mu, x, grid, dirichlet) - b;
	  }
	  else {
		  r = x - general_Fick_Lf(mu, x, grid, dirichlet) - b;
	  }
  }
  if (diffusionType==-1){
	  if (grid.regularGrid){
		  r = internalMu_inverse.asDiagonal() * x -
			  homogeneous_Lf(x, grid, dirichlet) -
			  internalMu_inverse.asDiagonal() * b;
	  }
	  else {
		  r = internalMu_inverse.asDiagonal() * x -
			  general_Homogeneous_Lf(x, grid, dirichlet) -
			  internalMu_inverse.asDiagonal() * b;
	  }
  }
  
  
  if (preconditionerType==0){
    if (debug) std::cout<<"Preconditioning by L^{-1}"<<"\n";
    plan.multiplyBy_pow_neg_A(-1.0,grid.lengthX,grid.lengthY);
  }
  else{
    if (debug) std::cout<<"Preconditioning by (I-delta\\cdot L)^{-1}"<<"\n";
    plan.multiplyBy_I_A_inverse(grid.lengthX,grid.lengthY,delta);
  }
  
  p=-y;
  
  double r_norm=std::pow(r.array().pow(2).array().sum(),0.5)/r.size();
  double alphaNum=(r.array()*y.array()).array().sum();
  int count=0;
  while ((r_norm>tol)& (count<nIter)){
    if (diffusionType==1){
		if (grid.regularGrid) {
			Ap = internalMu_inverse.asDiagonal() * p
				- homogeneous_Lf(p, grid, dirichlet);
		}
		else {
			Ap = internalMu_inverse.asDiagonal() * p
				- general_Homogeneous_Lf(p, grid, dirichlet);
		}
    }
    if (diffusionType==0){
		if (grid.regularGrid) {
			Ap = p - fick_Lf(mu, p, grid, dirichlet);
		}
		else {
			Ap = p - general_Fick_Lf(mu, p, grid, dirichlet);
		}
    }
    if (diffusionType==-1){
		if (grid.regularGrid) {
			Ap = internalMu_inverse.asDiagonal() * p
				- homogeneous_Lf(p, grid, dirichlet);
		}
		else {
			Ap = internalMu_inverse.asDiagonal() * p
				- general_Homogeneous_Lf(p, grid, dirichlet);
		}
    }
    double alpha=alphaNum/((p.adjoint()*Ap)(0));
    x=x+alpha*p;
    r=r+alpha*Ap;
    
    if (preconditionerType==0){
      if (debug) std::cout<<"Preconditioning by L^{-1}"<<"\n";
      plan.multiplyBy_pow_neg_A(-1.0,grid.lengthX,grid.lengthY);
    }
    else{
      if (debug) std::cout<<"Preconditioning by (I-delta\\cdot L)^{-1}"<<"\n";
      plan.multiplyBy_I_A_inverse(grid.lengthX,grid.lengthY,delta);
    }
    
    double beta=(r.array()*y.array()).array().sum()/alphaNum;
    p=-y+beta*p;
    alphaNum=alphaNum*beta;
    r_norm=std::pow(r.array().pow(2).array().sum(),0.5)/r.size();
    count+=1;
    
    if (debug){
      std::cout<<"iter "<<count<<"\n";
      std::cout<<"norm "<<r_norm<<"\n\n";
    }
  }
  if (diffusionType==1){
    return internalMu_inverse.asDiagonal()*x;
  }
  return x;
}

#endif

