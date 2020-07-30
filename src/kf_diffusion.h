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
  VectorXd gamma; //growth coefficients
  VectorXd longLat; //center of gaussian kernel initial conditions
  MatrixXd coords; //coordinates of grid cells
  double sigma,kappa; //spread and intensity of initial diffusion
  VectorXd eta; //coefficients for individual/deer level covariates, eg age, sex
  
  //covariate values
  MatrixXd X_diffusion;
  MatrixXd X_reaction;
  MatrixXd X_individual;
  
  //diffusion values and diffusion/growth coefficient values (functions of diffusion parameters)
  MatrixXd u,mu,lambda;
  
  //data and model properties
  VectorXi positive,cell,time; //cwd status, location, and time sampled for each deer
  bool individualCovariates; //include individual level covariates
  
  //properties of the grid
  bool rectangular; //indicator for whether the grid is rectangular
  VectorXi internalPointIndicator; //indicator for whether each point is internal or on the boundary
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
               VectorXd gamma_,
               VectorXd longLat_,
               double sigma_,
               double kappa_,
               MatrixXd coords_,
               MatrixXd X_reaction_,
               int rows_,
               int cols_,
               int nTime_,
               int diffusionType_,
               bool dirichlet_ = true,
               double lengthX_ = 1,
               double lengthY_ = 1){
    
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
    dirichlet=dirichlet_;
    lengthX=lengthX_;
    lengthY=lengthY_;
    rectangular=true;
    individualCovariates=false;
    computed=false;
    
    //number of internal points, full grid points
    nInternal=rows*cols;
    nFull=(rows+2)*(cols+2);
    
    u.setZero(nInternal,nTime); //initialize diffusion
    
    //initialize mu and lambda
    VectorXd muVec, lambdaVec;
    muVec.setZero(X_reaction.rows());
    muVec=muVec.array()+mu_0;
    
    lambdaVec.setZero(X_reaction.rows());
    lambdaVec = (X_reaction*gamma).array().exp();
    
    Map<MatrixXd> muMap(muVec.data(),nFull,nTime);
    Map<MatrixXd> lambdaMap(lambdaVec.data(),nFull,nTime);
    
    mu=muMap;
    lambda=lambdaMap;
  }
  
  KF_diffusion(double mu_0_,
               VectorXd gamma_,
               VectorXd longLat_,
               double sigma_,
               double kappa_,
               MatrixXd coords_,
               MatrixXd X_reaction_,
               VectorXi cell_,
               VectorXi positive_,
               VectorXi time_,
               int rows_,
               int cols_,
               int nTime_,
               int diffusionType_,
               bool dirichlet_ = true,
               double lengthX_ = 1,
               double lengthY_ = 1){
    
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
    dirichlet=dirichlet_;
    lengthX=lengthX_;
    lengthY=lengthY_;
    rectangular=true;
    individualCovariates=false;
    computed=false;
    cell=cell_;
    positive=positive_;
    time=time_;
    
    //number of internal points, full grid points
    nInternal=rows*cols;
    nFull=(rows+2)*(cols+2);
    
    u.setZero(nInternal,nTime); //initialize diffusion
    
    //initialize mu and lambda
    VectorXd muVec, lambdaVec;
    muVec.setZero(X_reaction.rows());
    muVec=muVec.array()+mu_0;
    
    lambdaVec.setZero(X_reaction.rows());
    lambdaVec = (X_reaction*gamma).array().exp();
    
    Map<MatrixXd> muMap(muVec.data(),nFull,nTime);
    Map<MatrixXd> lambdaMap(lambdaVec.data(),nFull,nTime);
    
    mu=muMap;
    lambda=lambdaMap;
  }
  
  KF_diffusion(double mu_0_,
               VectorXd gamma_,
               VectorXd longLat_,
               double sigma_,
               double kappa_,
               VectorXd eta_,
               MatrixXd coords_,
               MatrixXd X_reaction_,
               MatrixXd X_individual_,
               VectorXi cell_,
               VectorXi positive_,
               VectorXi time_,
               int rows_,
               int cols_,
               int nTime_,
               int diffusionType_,
               bool dirichlet_ = true,
               double lengthX_ = 1,
               double lengthY_ = 1){
    
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
    rows=rows_;
    cols=cols_;
    nTime=nTime_;
    diffusionType=diffusionType_;
    dirichlet=dirichlet_;
    lengthX=lengthX_;
    lengthY=lengthY_;
    rectangular=true;
    individualCovariates=true;
    computed=false;
    cell=cell_;
    positive=positive_;
    time=time_;
    
    //number of internal points, full grid points
    nInternal=rows*cols;
    nFull=(rows+2)*(cols+2);
    
    u.setZero(nInternal,nTime); //initialize diffusion
    
    //initialize mu and lambda
    VectorXd muVec, lambdaVec;
    muVec.setZero(X_reaction.rows());
    muVec=muVec.array()+mu_0;
    
    lambdaVec.setZero(X_reaction.rows());
    lambdaVec = (X_reaction*gamma).array().exp();
    
    Map<MatrixXd> muMap(muVec.data(),nFull,nTime);
    Map<MatrixXd> lambdaMap(lambdaVec.data(),nFull,nTime);
    
    mu=muMap;
    lambda=lambdaMap;
  }
  
  
  void setInitialConditions(){
    MatrixXd init(rows+2,cols+2);
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
    
    Map<MatrixXd> u0_map(u.col(0).data(), rows, cols);
    u0_map=init.block(1,1,rows,cols);
  }
  
  void computeDiffusion(int nIter=100,
                        double tol=1.0e-10){
    if (computed) return;
    
    VectorXd rhs(nInternal);
    setInitialConditions();
    for (int i=1;i<nTime;i++){
      
      //get growth coefficients for previous time point
      Map<MatrixXd> lambdaMatrix(lambda.col(i-1).data(),rows+2,cols+2);
      MatrixXd internalLambdaMatrix(lambdaMatrix.block(1,1,rows,cols));
      Map<VectorXd> lambda_i_1(internalLambdaMatrix.data(),rows*cols);
      
      //Kolmogorov-Fisher
      //U2-U1=HU2+Lambda*U1(1-U1)
      //(I-H)U2=U1+Lambda*U1(1-U1)
      rhs=u.col(i-1).array()+lambda_i_1.array()*u.col(i-1).array()*(1-u.col(i-1).array());
      VectorXd guess(u.col(i-1));
      VectorXd mu_i(mu.col(i-1));
      u.col(i)=solver(mu_i,guess,rhs,nIter,tol);
    }
    computed=true;
  }
  
  double logLikelihood(){
    computeDiffusion();
    double logLike=0;
    VectorXd shifts;
    if (individualCovariates) shifts=X_individual*eta;
    for (int i=0;i<positive.size();i++){
      int time_i=time(i);
      int cell_i=cell(i);
      int pos_i=positive(i);
      
      double probRaw_i=u(cell_i,time_i);
      double probAdjusted_i=probRaw_i;
      
      if (individualCovariates){  
        double logit_i=logit(probRaw_i);
        logit_i+=shifts(i);
        probAdjusted_i=expit(logit_i);
      }
      logLike+=pos_i*std::log(probAdjusted_i)+(1-pos_i)*std::log(1-probAdjusted_i);
    }
    return(logLike);
  }
  
  
  //derivative of diffusion wrt mu_0
  MatrixXd du_dmu(int nIter=100,
                  double tol=1.0e-10){
    computeDiffusion();
    
    MatrixXd du_dTheta;
    du_dTheta.setZero(nInternal,nTime);
    
    for (int i=0;i<nTime;i++){
      //(du_0/dTheta)=0
      if (i>0){
        
        //get growth coefficients for previous time point
        Map<MatrixXd> lambdaMatrix(lambda.col(i-1).data(),rows+2,cols+2);
        MatrixXd internalLambdaMatrix(lambdaMatrix.block(1,1,rows,cols));
        Map<VectorXd> lambda_i_1(internalLambdaMatrix.data(),rows*cols);
        
        
        VectorXd rhs(nInternal);
        rhs=u.col(i-1).array()+lambda_i_1.array()*u.col(i-1).array()*(1-u.col(i-1).array());
        VectorXd guess,temp;
        guess.setZero(nInternal);
        temp.setZero(nInternal);
        VectorXd mu_i(mu.col(i-1));
        du_dTheta.col(i)=solver(mu_i,guess,rhs,nIter,tol);
        
        temp=du_dTheta.col(i);
        du_dTheta.col(i)=homogeneous_L_f(temp,
                      rows,
                      cols,
                      dirichlet,
                      lengthX,
                      lengthY);
        
        temp=du_dTheta.col(i);
        du_dTheta.col(i)=solver(mu_i,guess,temp,nIter,tol);
        VectorXd rhs2;
        rhs2.setZero(nInternal);
        rhs2=du_dTheta.col(i-1).array()+lambda_i_1.array()*du_dTheta.col(i-1).array()+
          -2*lambda_i_1.array()*u.col(i-1).array()*du_dTheta.col(i-1).array();
          du_dTheta.col(i)=du_dTheta.col(i)+solver(mu_i,guess,rhs2,nIter,tol);
      }
    }
    
    return du_dTheta;
  }
  
  //derivative of diffusion with respect to gamma
  std::vector<MatrixXd> du_dgamma(int nIter=100,
                                  double tol=1.0e-10){
    computeDiffusion();
    std::vector<MatrixXd> du_dTheta_list;
    MatrixXd du_dTheta;
    du_dTheta.setZero(nInternal,nTime);
    
    int p=X_reaction.cols();
    
    for (int pInd=0;pInd<p;pInd++){
      du_dTheta.setZero(nInternal,nTime);
      for (int i=0;i<nTime;i++){
        //(du_0/dTheta)=0
        if (i>0){
          
          //get growth coefficients for previous time point
          Map<MatrixXd> lambdaMatrix(lambda.col(i-1).data(),rows+2,cols+2);
          MatrixXd internalLambdaMatrix(lambdaMatrix.block(1,1,rows,cols));
          Map<VectorXd> lambda_i_1(internalLambdaMatrix.data(),rows*cols);
          
          //get covariate values for previous time point
          Map<MatrixXd> XpMatrix=time_i_covariate(X_reaction,i-1,pInd);
          MatrixXd internalXpMatrix(XpMatrix.block(1,1,rows,cols));
          Map<VectorXd> Xp_i_1(internalXpMatrix.data(),rows*cols);
          
          VectorXd rhs2,guess;
          guess.setZero(nInternal);
          
          rhs2.setZero(nInternal);
          rhs2=du_dTheta.col(i-1).array()+
            Xp_i_1.array()*lambda_i_1.array()*u.col(i-1).array()*(1-u.col(i-1).array()).array()+
            lambda_i_1.array()*du_dTheta.col(i-1).array()+
            -2*lambda_i_1.array()*u.col(i-1).array()*du_dTheta.col(i-1).array();
            VectorXd mu_i(mu.col(i-1));
          du_dTheta.col(i)=solver(mu_i,guess,rhs2,nIter,tol);
        }
      }
      du_dTheta_list.push_back(du_dTheta);
    }
    return du_dTheta_list;
  }
  
  std::vector<MatrixXd> du_dlongLat(int nIter=100,
                                    double tol=1.0e-10){
    computeDiffusion();
    
    MatrixXd du_dTheta;
    du_dTheta.setZero(nInternal,nTime);
    
    
    MatrixXd initFull(rows+2,cols+2);
    Map<MatrixXd> du_dTheta0_map(du_dTheta.col(0).data(),rows,cols);
    
    
    std::vector<MatrixXd> du_dTheta_list;
    
    for (int llInd=0;llInd<2;llInd++){
      du_dTheta.setZero(nInternal,nTime);
      
      
      for (int i=0;i<nTime;i++){
        if (i==0){
          int count=0;
          for (int cInd=0;cInd<cols+2;cInd++){
            for (int rInd=0;rInd<rows+2;rInd++){
              VectorXd d(2);
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
          Map<MatrixXd> lambdaMatrix(lambda.col(i-1).data(),rows+2,cols+2);
          MatrixXd internalLambdaMatrix(lambdaMatrix.block(1,1,rows,cols));
          Map<VectorXd> lambda_i_1(internalLambdaMatrix.data(),rows*cols);
          
          
          VectorXd rhs2,guess;
          guess.setZero(nInternal);
          
          rhs2.setZero(nInternal);
          rhs2=du_dTheta.col(i-1).array()+
            lambda_i_1.array()*du_dTheta.col(i-1).array()+
            -2*lambda_i_1.array()*u.col(i-1).array()*du_dTheta.col(i-1).array();
            VectorXd mu_i(mu.col(i-1));
            du_dTheta.col(i)=solver(mu_i,guess,rhs2,nIter,tol);
        }
      }
      
      du_dTheta_list.push_back(du_dTheta);
    }
    return du_dTheta_list;
  }
  
  MatrixXd du_dsigma(int nIter=100,
                            double tol=1.0e-10){
    computeDiffusion();
    
    MatrixXd du_dTheta;
    du_dTheta.setZero(nInternal,nTime);
    
    
    MatrixXd initFull(rows+2,cols+2);
    Map<MatrixXd> du_dTheta0_map(du_dTheta.col(0).data(),rows,cols);
    
    
    du_dTheta.setZero(nInternal,nTime);
    
    
    for (int i=0;i<nTime;i++){
      if (i==0){
        int count=0;
        for (int cInd=0;cInd<cols+2;cInd++){
          for (int rInd=0;rInd<rows+2;rInd++){
            VectorXd d(2);
            d(0)=coords(count,0)-longLat(0);
            d(1)=coords(count,1)-longLat(1);
            double dist2=std::pow(d(0),2)+std::pow(d(1),2);
            initFull(rInd,cInd)=expit(kappa)*std::exp(-dist2/(2*std::pow(sigma,2)))*dist2/std::pow(sigma,3);
            count=count+1;
          }
        }
        du_dTheta0_map=initFull.block(1,1,rows,cols);
      }
      else{
        Map<MatrixXd> lambdaMatrix(lambda.col(i-1).data(),rows+2,cols+2);
        MatrixXd internalLambdaMatrix(lambdaMatrix.block(1,1,rows,cols));
        Map<VectorXd> lambda_i_1(internalLambdaMatrix.data(),rows*cols);
        
        
        VectorXd rhs2,guess;
        guess.setZero(nInternal);
        
        rhs2.setZero(nInternal);
        rhs2=du_dTheta.col(i-1).array()+
          lambda_i_1.array()*du_dTheta.col(i-1).array()+
          -2*lambda_i_1.array()*u.col(i-1).array()*du_dTheta.col(i-1).array();
          VectorXd mu_i(mu.col(i-1));
          du_dTheta.col(i)=solver(mu_i,guess,rhs2,nIter,tol);
      }
    }
    
    return du_dTheta;
  }
  MatrixXd du_dkappa(int nIter=100,
                            double tol=1.0e-10){
    computeDiffusion();
    
    MatrixXd du_dTheta;
    du_dTheta.setZero(nInternal,nTime);
    
    
    MatrixXd initFull(rows+2,cols+2);
    Map<MatrixXd> du_dTheta0_map(du_dTheta.col(0).data(),rows,cols);
    
    
    du_dTheta.setZero(nInternal,nTime);
    
    
    for (int i=0;i<nTime;i++){
      if (i==0){
        int count=0;
        for (int cInd=0;cInd<cols+2;cInd++){
          for (int rInd=0;rInd<rows+2;rInd++){
            VectorXd d(2);
            d(0)=coords(count,0)-longLat(0);
            d(1)=coords(count,1)-longLat(1);
            double dist2=std::pow(d(0),2)+std::pow(d(1),2);
            initFull(rInd,cInd)=expit(kappa)*(1-expit(kappa))*std::exp(-dist2/(2*std::pow(sigma,2)));
            count=count+1;
          }
        }
        du_dTheta0_map=initFull.block(1,1,rows,cols);
      }
      else{
        Map<MatrixXd> lambdaMatrix(lambda.col(i-1).data(),rows+2,cols+2);
        MatrixXd internalLambdaMatrix(lambdaMatrix.block(1,1,rows,cols));
        Map<VectorXd> lambda_i_1(internalLambdaMatrix.data(),rows*cols);
        
        
        VectorXd rhs2,guess;
        guess.setZero(nInternal);
        
        rhs2.setZero(nInternal);
        rhs2=du_dTheta.col(i-1).array()+
          lambda_i_1.array()*du_dTheta.col(i-1).array()+
          -2*lambda_i_1.array()*u.col(i-1).array()*du_dTheta.col(i-1).array();
          VectorXd mu_i(mu.col(i-1));
          du_dTheta.col(i)=solver(mu_i,guess,rhs2,nIter,tol);
      }
    }
    
    return du_dTheta;
  }
  
  VectorXd dl_dEta(){
    computeDiffusion();
    
    VectorXd derivative;
    derivative.setZero(X_individual.cols());
    
    VectorXd shifts;
    if (individualCovariates) shifts=X_individual*eta;
    for (int i=0;i<positive.size();i++){
      int time_i=time(i);
      int cell_i=cell(i);
      int pos_i=positive(i);
      
      double probRaw_i=u(cell_i,time_i);
      double logit_i=logit(probRaw_i);
      
      logit_i+=shifts(i);
      double probAdjusted_i=expit(logit_i);
      for (int pInd=0;pInd<X_individual.cols();pInd++){
        derivative(pInd)=derivative(pInd)+(pos_i-probAdjusted_i)*X_individual(i,pInd);
      }
    }
    return(derivative);
  }
  
  double dl_du(MatrixXd& du_dTheta){
    double derivative=0;

    // std::cout<<"positives: "<<positive.size()<<"\n";
    // std::cout<<"min "<<positive.minCoeff()<<"\n";
    // std::cout<<"max "<<positive.maxCoeff()<<"\n";
    // std::cout<<"cell: "<<cell.size()<<"\n";
    // std::cout<<"min "<<cell.minCoeff()<<"\n";
    // std::cout<<"max "<<cell.maxCoeff()<<"\n";
    // std::cout<<"time: "<<time.size()<<"\n";
    // std::cout<<"min "<<time.minCoeff()<<"\n";
    // std::cout<<"max "<<time.maxCoeff()<<"\n";
    // Rcout<<"positives: "<<positive.size()<<"\n";
    // Rcout<<"min "<<positive.minCoeff()<<"\n";
    // Rcout<<"max "<<positive.maxCoeff()<<"\n";
    // Rcout<<"cell: "<<cell.size()<<"\n";
    // Rcout<<"min "<<cell.minCoeff()<<"\n";
    // Rcout<<"max "<<cell.maxCoeff()<<"\n";
    // Rcout<<"time: "<<time.size()<<"\n";
    // Rcout<<"min "<<time.minCoeff()<<"\n";
    // Rcout<<"max "<<time.maxCoeff()<<"\n";

    // 
    VectorXd shifts;
    if (individualCovariates) shifts=X_individual*eta;
    for (int i=0;i<positive.size();i++){
      int time_i=time(i);
      int cell_i=cell(i);
      int pos_i=positive(i);


      double probRaw_i=u(cell_i,time_i);
      double probAdjusted_i=probRaw_i;

      if (individualCovariates){
        double logit_i=logit(probRaw_i);
        logit_i+=shifts(i);
        probAdjusted_i=expit(logit_i);
      }

      double rawVar=probRaw_i*(1-probRaw_i);
      double derivRaw_i=du_dTheta(cell_i,time_i);
      derivative+=(pos_i-probAdjusted_i)/rawVar*derivRaw_i;

    }
    return(derivative);
  }
  
  MatrixXd padDiffusion(MatrixXd& u_){
    MatrixXd u_Full(nFull,nTime);
    u_Full.setZero();
    
    for (int i=0;i<nTime;i++){
      Map<MatrixXd> u_i(u_.col(i).data(),rows,cols);
      Map<MatrixXd> u_Full_i(u_Full.col(i).data(),rows+2,cols+2);
      u_Full_i.block(1,1,rows,cols)=u_i;
      
      if (!dirichlet){
        for (int rowInd=1;rowInd<rows+1;rowInd++){
          u_Full_i(rowInd,0)=u_Full_i(rowInd,1);
          u_Full_i(rowInd,cols+1)=u_Full_i(rowInd,cols);
        }
        for (int colInd=1;colInd<cols+1;colInd++){
          u_Full_i(0,colInd)=u_Full_i(1,colInd);
          u_Full_i(rows+1,colInd)=u_Full_i(rows,colInd);
        }
      }
    }
    return(u_Full);
  }
  
  VectorXd differentiateLogLikelihood(){
    VectorXd x;
    if (individualCovariates){
      //1 diffusion coefficient, p1 reaction coefficients, p2 individual coefficients
      //2 long lat coefficients, sigma, kappa
      x.setZero(1+X_reaction.cols()+X_individual.cols()+4);
    }
    else{
      x.setZero(1+X_reaction.cols()+4);
    }
    
    //diffusion
    int count=0;
    Eigen::MatrixXd du_dTheta(nInternal,nTime);
    du_dTheta=du_dmu();
    x(count)=dl_du(du_dTheta);
    count+=1;
    // 
    //gamma
    std::vector<MatrixXd> du_dTheta_list;
    du_dTheta_list=du_dgamma();
    for (int i=0;i<du_dTheta_list.size();i++){
      x(count)=dl_du(du_dTheta_list[i]);
      count+=1;
    }
    // 
    //longLat
    du_dTheta_list=du_dlongLat();
    for (int i=0;i<du_dTheta_list.size();i++){
      x(count)=dl_du(du_dTheta_list[i]);
      count+=1;
    }
    //
    // //sigma
    du_dTheta=du_dsigma();
    x(count)=dl_du(du_dTheta);
    count+=1;

    // //kappa
    du_dTheta=du_dkappa();
    x(count)=dl_du(du_dTheta);
    count+=1;
    // 
    // //eta
    if (individualCovariates){
      x.segment(count,X_individual.cols())=dl_dEta();
    }
    // 
    return x;
  }
  
  MatrixXd get_u(){
    return u;
  }
  
  Map<MatrixXd> time_i_covariate(MatrixXd& X_, int i=0,int p=0){
    Map<MatrixXd> X_p(X_.col(p).data(),nFull,nTime);
    Map<MatrixXd> X_p_i(X_p.col(i).data(),rows+2,cols+2);
    return X_p_i;
  }
  
  VectorXd solver(VectorXd& mu_i,
                         VectorXd& guess,
                         VectorXd& rhs,
                         int nIter,
                         double tol){
    return invert(rows,
         cols,
         mu_i,
         guess,
         rhs,
         diffusionType,
         nIter,
         tol,
         lengthX,
         lengthY,
         1,//preconditioner type
         dirichlet,
         false);//debug
  }
};





