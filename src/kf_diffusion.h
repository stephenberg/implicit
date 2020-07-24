#include <Rcpp.h>
#include <RcppEigen.h>
#include "fftPlan.h"
#include "convenienceFunctions.h"
#include "inversionFunctions.h"
using namespace Rcpp;
using namespace Eigen;

class KF_diffusion
{
private:
  
  //diffusion parameters
  double mu_0; //diffusion coefficient
  Eigen::VectorXd gamma; //growth coefficients
  Eigen::VectorXd longLat; //center of gaussian kernel initial conditions
  Eigen::MatrixXd coords; //coordinates of grid cells
  double sigma,kappa; //spread and intensity of initial diffusion
  
  //covariate values
  Eigen::MatrixXd X_reaction;
  
  //diffusion values and diffusion/growth coefficient values (functions of diffusion parameters)
  Eigen::MatrixXd u;
  Eigen::Map<Eigen::MatrixXd> mu,lambda;
  
  //data and model properties
  Eigen::VectorXi positive,cell,time; //cwd status, location, and time sampled for each deer
  Eigen::VectorXd age; //age of each deer
  Eigen::VectorXi sex; //sex of each deer
  bool include_AgeSex; //include individual level covariates
  
  //properties of the grid
  bool rectangular; //indicator for whether the grid is rectangular
  Eigen::VectorXi internalPointIndicator; //indicator for whether each point is internal or on the boundary
  int rows,cols; //number of internal rows, internal columns
  int nInternal,nFull;  //number of locations in internal, full grid 
  int nTime; //number of time steps
  double lengthX,lengthY; //horizontal and vertical lengths of the grid, in units (m, km, miles, ...)
  
  //properties of the diffusion
  bool dirichlet; //false: use zero-flux boundary conditions, true: use Dirichlet boundary conditions
  int diffusionType; //-1,0,1
  
  
public:
  KF_diffusion(double mu_0_,
               Eigen::VectorXd gamma_,
               Eigen::VectorXd longLat_,
               double sigma_,
               double kappa_,
               Eigen::MatrixXd coords_,
               Eigen::MatrixXd X_reaction_,
               Eigen::VectorXi positive_,
               Eigen::VectorXi cell_,
               Eigen::VectorXi time_,
               int rows_,
               int cols_,
               int nTime_,
               int diffusionType_,
               double lengthX_=1,
               double lengthY_=1){
    
    //assignments
    mu_0=mu_0_;
    gamma=gamma_;
    longLat=longLat_;
    sigma=sigma_;
    kappa=kappa_;
    coords=coords_;
    X_reaction=X_reaction_;
    positive=positive_;
    cell=cell_;
    time=time_;
    rows=rows_;
    cols=cols_;
    nTime=nTime_;
    diffusionType=diffusionType_;
    lengthX=lengthX_;
    lengthY=lengthY_;
    rectangular=true;
    include_AgeSex=false;
    dirichlet=true;
    
    //number of internal points, full grid points
    nInternal=rows*cols;
    nFull=(rows+2)*(cols+2);
  
    u.setZero(nInternal,nTime); //initialize diffusion
    
    //initialize mu and lambea
    Eigen::VectorXd muVec, lambdaVec;
    muVec.setZero(X_reaction.rows());
    muVec=muVec.array()+mu_0;
    
    lambdaVec.setZero(X_reaction.rows());
    lambdaVec = (X_reaction*gamma).array().exp();
    
    mu(muVec.data(),nFull,nTime);
    lambda(lambdaVec.data(),nFull,nTime);
  }
  
  void setInitialConditions(){
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
    
    Eigen::Map<Eigen::MatrixXd> u0_map(u.col(0).data(), rows, cols);
    u0_map=init.block(1,1,rows,cols);
  }
  
  void computeDiffusion(int nIter,
                        double tol){
    Eigen::VectorXd rhs(nInternal);
    
    for (int i=1;i<nTime;i++){
      //Kolmogorov-Fisher
      //U2-U1=HU2+Lambda*U1(1-U1)
      //(I-H)U2=U1+Lambda*U1(1-U1)
      rhs=u.col(i-1).array()+lambda.col(i-1).array()*u.col(i-1).array()*(1-u.col(i-1).array());
      Eigen::VectorXd guess(u.col(i-1));
      Eigen::VectorXd mu_i(mu.col(i-1));
      u.col(i)=invert(rows,
            cols,
            mu_i,
            guess,
            rhs,
            diffusionType,
            nIter,
            lengthX,
            lengthY);
    }
  }
  
  double logLikelihood(){
    return 0;
  }
  
  Eigen::VectorXd dl_dmu(){
    Eigen::VectorXd x;
    return x;
  }
  
  Eigen::VectorXd dl_dgamma(){
    Eigen::VectorXd x;
    return x;
  }
  
  Eigen::VectorXd dl_dlongLat(){
    Eigen::VectorXd x;
    return x;
  }
  
  Eigen::VectorXd dl_dsigma(){
    Eigen::VectorXd x;
    return x;
  }
  
  Eigen::VectorXd dl_dkappa(){
    Eigen::VectorXd x;
    return x;
  }
  
  Eigen::VectorXd differentiateLogLikelihood(){
    Eigen::VectorXd x;
    return x;
  }
};





