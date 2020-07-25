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
  Eigen::VectorXd eta; //coefficients for individual/deer level covariates, eg age, sex
  
  //covariate values
  Eigen::MatrixXd X_diffusion;
  Eigen::MatrixXd X_reaction;
  Eigen::MatrixXd X_individual;
  
  //diffusion values and diffusion/growth coefficient values (functions of diffusion parameters)
  Eigen::MatrixXd u,mu,lambda;
  
  //data and model properties
  Eigen::VectorXi positive,cell,time; //cwd status, location, and time sampled for each deer
  bool individualCovariates; //include individual level covariates
  
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
  bool computed; //indicator for whether the diffusion has been computed
  
  
public:
  
  KF_diffusion(double mu_0_,
               Eigen::VectorXd gamma_,
               Eigen::VectorXd longLat_,
               double sigma_,
               double kappa_,
               Eigen::MatrixXd coords_,
               Eigen::MatrixXd X_reaction_,
               int rows_,
               int cols_,
               int nTime_,
               int diffusionType_,
               double lengthX_ = 1,
               double lengthY_= 1){
    
    //assignments
    mu_0=mu_0_;
    gamma=gamma_;
    longLat=longLat_;
    sigma=sigma_;
    kappa=kappa_;
    coords=coords_;
    X_reaction=X_reaction_;
    rows=rows_;
    cols=cols_;
    nTime=nTime_;
    diffusionType=diffusionType_;
    lengthX=lengthX_;
    lengthY=lengthY_;
    rectangular=true;
    individualCovariates=false;
    dirichlet=true;
    computed=false;
    
    //number of internal points, full grid points
    nInternal=rows*cols;
    nFull=(rows+2)*(cols+2);
    
    u.setZero(nInternal,nTime); //initialize diffusion
    
    //initialize mu and lambda
    Eigen::VectorXd muVec, lambdaVec;
    muVec.setZero(X_reaction.rows());
    muVec=muVec.array()+mu_0;
    
    lambdaVec.setZero(X_reaction.rows());
    lambdaVec = (X_reaction*gamma).array().exp();
    
    Eigen::Map<Eigen::MatrixXd> muMap(muVec.data(),nFull,nTime);
    Eigen::Map<Eigen::MatrixXd> lambdaMap(lambdaVec.data(),nFull,nTime);
    
    mu=muMap;
    lambda=lambdaMap;
  }
  
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
               double lengthX_ = 1,
               double lengthY_= 1){
    
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
    individualCovariates=false;
    dirichlet=true;
    computed=false;

    //number of internal points, full grid points
    nInternal=rows*cols;
    nFull=(rows+2)*(cols+2);

    u.setZero(nInternal,nTime); //initialize diffusion

    //initialize mu and lambda
    Eigen::VectorXd muVec, lambdaVec;
    muVec.setZero(X_reaction.rows());
    muVec=muVec.array()+mu_0;

    lambdaVec.setZero(X_reaction.rows());
    lambdaVec = (X_reaction*gamma).array().exp();
    
    Eigen::Map<Eigen::MatrixXd> muMap(muVec.data(),nFull,nTime);
    Eigen::Map<Eigen::MatrixXd> lambdaMap(lambdaVec.data(),nFull,nTime);
    
    mu=muMap;
    lambda=lambdaMap;
    
    computed=false;
    individualCovariates=false;
  }
  
  KF_diffusion(double mu_0_,
               Eigen::VectorXd gamma_,
               Eigen::VectorXd longLat_,
               double sigma_,
               double kappa_,
               Eigen::VectorXd eta_,
               Eigen::MatrixXd coords_,
               Eigen::MatrixXd X_reaction_,
               Eigen::MatrixXd X_individual_,
               Eigen::VectorXi positive_,
               Eigen::VectorXi cell_,
               Eigen::VectorXi time_,
               int rows_,
               int cols_,
               int nTime_,
               int diffusionType_,
               double lengthX_ = 1,
               double lengthY_= 1){
    
    //assignments
    mu_0=mu_0_;
    gamma=gamma_;
    longLat=longLat_;
    sigma=sigma_;
    kappa=kappa_;
    eta=eta_;
    coords=coords_;
    X_reaction=X_reaction_;
    X_individual=X_individual_;
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
    individualCovariates=true;
    dirichlet=true;
    computed=false;
    
    //number of internal points, full grid points
    nInternal=rows*cols;
    nFull=(rows+2)*(cols+2);
    
    u.setZero(nInternal,nTime); //initialize diffusion
    
    //initialize mu and lambda
    Eigen::VectorXd muVec, lambdaVec;
    muVec.setZero(X_reaction.rows());
    muVec=muVec.array()+mu_0;
    
    lambdaVec.setZero(X_reaction.rows());
    lambdaVec = (X_reaction*gamma).array().exp();
    
    Eigen::Map<Eigen::MatrixXd> muMap(muVec.data(),nFull,nTime);
    Eigen::Map<Eigen::MatrixXd> lambdaMap(lambdaVec.data(),nFull,nTime);
    
    mu=muMap;
    lambda=lambdaMap;
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
  
  void computeDiffusion(int nIter=100,
                        double tol=1.0e-10){
    if (computed) return;
    
    Eigen::VectorXd rhs(nInternal);
    setInitialConditions();
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
            tol,
            lengthX,
            lengthY);
    }
    computed=true;
  }
  
  double logLikelihood(){
    double logLike=0;
    for (int i=0;i<positive.size();i++){
      int time_i=time(i);
      int cell_i=cell(i);
      int pos_i=positive(i);
      
      double probRaw_i=u(cell_i,time_i);
      double probAdjusted_i=probRaw_i;
      
      if (individualCovariates){  
        double logit_i=logit(probRaw_i);
        logit_i+=(X_individual.row(i).array()*eta.array()).array().sum();
        probAdjusted_i=expit(logit_i);
      }
      
      logLike+=pos_i*std::log(probAdjusted_i)+(1-pos_i)*std::log(1-probAdjusted_i);
    }
    return(logLike);
  }
  
  
  
  //wrt mu_0
  Eigen::VectorXd dlogLike_dmu(int nIter=100,
                               double tol=1.0e-10){
    computeDiffusion();
    
    Eigen::MatrixXd du_dTheta;
    du_dTheta.setZero(nInternal,nTime);
    Eigen::VectorXd dl_dmu;
    dl_dmu.setZero(1);
    
    for (int i=0;i<nTime;i++){
      //(du_0/dTheta)=0
      if (i>0){
        Eigen::VectorXd rhs(nInternal);
        rhs=u.col(i-1).array()+lambda.col(i-1).array()*u.col(i-1).array()*(1-u.col(i-1).array());
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
                      tol,
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
          -2*lambda.col(i-1).array()*u.col(i-1).array()*du_dTheta.col(i-1).array();
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
    
    dl_dmu(0)=dl_du(du_dTheta);
    return dl_dmu;
  }
  
  //with respect to gamma
  Eigen::VectorXd dlogLike_dgamma(int nIter=100,
                                  double tol=1.0e-10){
    computeDiffusion();
    
    Eigen::MatrixXd du_dTheta;
    du_dTheta.setZero(nInternal,nTime);
    
    int p=X_reaction.cols();
    Eigen::VectorXd dl_dgamma(p);
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
            X_reaction.col(pInd).array()*lambda.col(i-1).array()*u.col(i-1).array()*(1-u.col(i-1).array()).array()+
            lambda.col(i-1).array()*du_dTheta.col(i-1).array()+
            -2*lambda.col(i-1).array()*u.col(i-1).array()*du_dTheta.col(i-1).array();
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

      dl_dgamma(pInd)=dl_du(du_dTheta);
    }
    return dl_dgamma;
  }
  
  Eigen::VectorXd dlogLike_dlongLat(int nIter=100,
                                    double tol=1.0e-10){
    computeDiffusion();
    
    Eigen::MatrixXd du_dTheta;
    du_dTheta.setZero(nInternal,nTime);
    
    Eigen::VectorXd dl_dlongLat(2);
    dl_dlongLat.setZero();
    
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
            -2*lambda.col(i-1).array()*u.col(i-1).array()*du_dTheta.col(i-1).array();
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

      dl_dlongLat(llInd)=dl_du(du_dTheta);
    }
    return dl_dlongLat;
  }
  
  //with respect to sigma  
  Eigen::VectorXd dlogLike_dsigma(int nIter=100,
                                  double tol=1.0e-10){
    computeDiffusion();
    
    Eigen::MatrixXd du_dTheta;
    du_dTheta.setZero(nInternal,nTime);
    
    Eigen::VectorXd dl_dsigma(1);
    dl_dsigma.setZero();
    
    Eigen::MatrixXd initFull(rows+2,cols+2);
    Eigen::Map<Eigen::MatrixXd> du_dTheta0_map(du_dTheta.col(0).data(),rows,cols);
    
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
          -2*lambda.col(i-1).array()*u.col(i-1).array()*du_dTheta.col(i-1).array();
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
    
    dl_dsigma(0)=dl_du(du_dTheta);
    
    return dl_dsigma;
  }

  //with respect to kappa  
  Eigen::VectorXd dlogLike_dkappa(int nIter=100,
                                  double tol=1.0e-10){
    computeDiffusion();
    
    Eigen::MatrixXd du_dTheta;
    du_dTheta.setZero(nInternal,nTime);
    
    Eigen::VectorXd dl_dkappa(1);
    dl_dkappa.setZero();
    
    Eigen::MatrixXd initFull(rows+2,cols+2);
    Eigen::Map<Eigen::MatrixXd> du_dTheta0_map(du_dTheta.col(0).data(),rows,cols);
    
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
          -2*lambda.col(i-1).array()*u.col(i-1).array()*du_dTheta.col(i-1).array();
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
    
    dl_dkappa(0)=dl_du(du_dTheta);
    
    return dl_dkappa;
  }
  
  Eigen::VectorXd dlogLike_dEta(int nIter=100,
                                double tol=1.0e-10){
    computeDiffusion();
    
    Eigen::VectorXd derivative;
    derivative.setZero(X_individual.cols());
    for (int i=0;i<positive.size();i++){
      int time_i=time(i);
      int cell_i=cell(i);
      int pos_i=positive(i);
      
      double probRaw_i=u(cell_i,time_i);
      double logit_i=logit(probRaw_i);
      
      logit_i+=(X_individual.row(i).array()*eta.array()).array().sum();
      double probAdjusted_i=expit(logit_i);
      
      derivative=derivative+(pos_i-probAdjusted_i)*X_individual.row(i);
    }
    return(derivative);
  }
  
  double dl_du(Eigen::MatrixXd& du_dTheta){
    double derivative=0;
    for (int i=0;i<positive.size();i++){
      int time_i=time(i);
      int cell_i=cell(i);
      int pos_i=positive(i);
      
      
      double probRaw_i=u(cell_i,time_i);
      double probAdjusted_i=probRaw_i;
      
      if (individualCovariates){  
        double logit_i=logit(probRaw_i);
        logit_i+=(X_individual.row(i).array()*eta.array()).array().sum();
        probAdjusted_i=expit(logit_i);
      }
      
      double rawVar=probRaw_i*(1-probRaw_i);
      double derivRaw_i=du_dTheta(cell_i,time_i);
      derivative+=(pos_i-probAdjusted_i)/rawVar*derivRaw_i;
      
    }
    return(derivative);
  }
  
  Eigen::VectorXd differentiateLogLikelihood(){
    Eigen::VectorXd x;
    return x;
  }
};





