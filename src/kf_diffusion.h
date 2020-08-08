#include <Rcpp.h>
#include <RcppEigen.h>
#include "fftPlan.h"
#include "convenienceFunctions.h"
#include "inversionFunctions.h"
#include "grid.h"

using namespace Rcpp;
using namespace Eigen;

class KF_diffusion
{
private:
  
  //diffusion parameters
  double mu_0; //diffusion coefficient
  VectorXd alpha; //diffusion coefficients
  VectorXd gamma; //growth coefficients
  VectorXd longLat; //center of gaussian kernel initial conditions
  MatrixXd coords; //coordinates of grid cells
  double sigma,kappa; //spread and intensity of initial diffusion
  VectorXd eta; //coefficients for individual/deer level covariates, eg age, sex
  
  //covariate values
  MatrixXd X_diffusion;
  MatrixXd X_reaction;
  MatrixXd X_individual;
  
  //data and model properties
  int nTime; //number of time steps
  VectorXi positive,cell,time; //cwd status, location, and time sampled for each deer
  bool individualCovariates; //include individual level covariates
  bool diffusionCovariates; //include covariates with the diffusion parameter
  
  //properties of the grid
  Grid grid;

  //properties of the diffusion
  bool dirichlet; //false: use zero-flux boundary conditions, true: use Dirichlet boundary conditions
  int diffusionType; //-1,0,1; use 0 for Fickian, 1 for ecological
  bool computed; //indicator for whether the diffusion has been computed

  
  //diffusion values and diffusion/growth coefficient values (functions of diffusion parameters)
  MatrixXd u, mu, lambda;

public:
  
  //for computing diffusion given parameter values, without passing in observed data or internal points indicator
  KF_diffusion(double mu_0_,
	           VectorXd alpha_,
               VectorXd gamma_,
               VectorXd longLat_,
               double sigma_,
               double kappa_,
               MatrixXd coords_,
	           MatrixXd X_diffusion_,
               MatrixXd X_reaction_,
               int rows_,
               int cols_,
               int nTime_,
               int diffusionType_,
               bool dirichlet_ = true,
               double lengthX_ = 1,
               double lengthY_ = 1):grid(rows_,cols_,lengthX_,lengthY_){
    
    //assignments
    mu_0=mu_0_;
	alpha = alpha_;
    gamma=gamma_;
    longLat=longLat_;
    sigma=sigma_;
    kappa=kappa_;
    coords=coords_;
	X_diffusion = X_diffusion_;
    X_reaction=X_reaction_;
    nTime=nTime_;
    diffusionType=diffusionType_;
    dirichlet=dirichlet_;
	diffusionCovariates = true;
    individualCovariates=false;
    computed=false;
    
    u.setZero(grid.nInternal,nTime); //initialize diffusion
    
    //initialize mu and lambda
    VectorXd muVec, lambdaVec;
    muVec.setZero(X_reaction.rows());
	muVec = (X_diffusion * alpha).array().exp().array() * mu_0;
    
    lambdaVec.setZero(X_reaction.rows());
    lambdaVec = (X_reaction*gamma).array().exp();
    
    Map<MatrixXd> muMap(muVec.data(),grid.nFull,nTime);
    Map<MatrixXd> lambdaMap(lambdaVec.data(),grid.nFull,nTime);
    
    mu=muMap;
    lambda=lambdaMap;
  }

  //for computing diffusion given parameter values, without passing in observed data but passing in internal points indicator
  KF_diffusion(double mu_0_,
	  VectorXd alpha_,
	  VectorXd gamma_,
	  VectorXd longLat_,
	  double sigma_,
	  double kappa_,
	  MatrixXd coords_,
	  MatrixXd X_diffusion_,
	  MatrixXd X_reaction_,
	  int rows_,
	  int cols_,
	  VectorXi internalPoints_,
	  int nTime_,
	  int diffusionType_,
	  bool dirichlet_ = true,
	  double lengthX_ = 1,
	  double lengthY_ = 1) :grid(rows_, cols_, internalPoints_, lengthX_, lengthY_) {

	  //assignments
	  mu_0 = mu_0_;
	  alpha = alpha_;
	  gamma = gamma_;
	  longLat = longLat_;
	  sigma = sigma_;
	  kappa = kappa_;
	  coords = coords_;
	  X_diffusion = X_diffusion_;
	  X_reaction = X_reaction_;
	  nTime = nTime_;
	  diffusionType = diffusionType_;
	  dirichlet = dirichlet_;
	  diffusionCovariates = true;
	  individualCovariates = false;
	  computed = false;

	  u.setZero(grid.nInternal, nTime); //initialize diffusion

	  //initialize mu and lambda
	  VectorXd muVec, lambdaVec;
	  muVec.setZero(X_reaction.rows());
	  muVec = (X_diffusion * alpha).array().exp().array() * mu_0;

	  lambdaVec.setZero(X_reaction.rows());
	  lambdaVec = (X_reaction * gamma).array().exp();

	  Map<MatrixXd> muMap(muVec.data(), grid.nFull, nTime);
	  Map<MatrixXd> lambdaMap(lambdaVec.data(), grid.nFull, nTime);

	  mu = muMap;
	  lambda = lambdaMap;
  }
  
  //for computing diffusion and loglikelihood given parameter values, without passing in internal points indicator
  KF_diffusion(double mu_0_,
	           VectorXd alpha_,
               VectorXd gamma_,
               VectorXd longLat_,
               double sigma_,
               double kappa_,
               VectorXd eta_,
               MatrixXd coords_,
	           MatrixXd X_diffusion_,
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
               double lengthY_ = 1):grid(rows_, cols_, lengthX_, lengthY_) {
    
    //assignments
    mu_0=mu_0_;
	alpha = alpha_;
    gamma=gamma_;
    longLat=longLat_;
    sigma=sigma_;
    kappa=kappa_;
    eta=eta_;
    coords=coords_;
	X_diffusion = X_diffusion_;
    X_reaction=X_reaction_;
    X_individual=X_individual_;
    nTime=nTime_;
    diffusionType=diffusionType_;
    dirichlet=dirichlet_;
	diffusionCovariates = true;
    individualCovariates=true;
    computed=false;
    cell=cell_;
    positive=positive_;
    time=time_;
    
   
    u.setZero(grid.nInternal,nTime); //initialize diffusion
    
    //initialize mu and lambda
    VectorXd muVec, lambdaVec;
    muVec.setZero(X_reaction.rows());
    muVec=muVec.array()+mu_0;
	muVec = muVec.array() * (X_diffusion * alpha).array().exp().array();
    
    lambdaVec.setZero(X_reaction.rows());
    lambdaVec = (X_reaction*gamma).array().exp();
    
    Map<MatrixXd> muMap(muVec.data(),grid.nFull,nTime);
    Map<MatrixXd> lambdaMap(lambdaVec.data(),grid.nFull,nTime);
    
    mu=muMap;
    lambda=lambdaMap;
  }

  //for computing diffusion and loglikelihood given parameter values, passing in internal points indicator
  KF_diffusion(double mu_0_,
	  VectorXd alpha_,
	  VectorXd gamma_,
	  VectorXd longLat_,
	  double sigma_,
	  double kappa_,
	  VectorXd eta_,
	  MatrixXd coords_,
	  MatrixXd X_diffusion_,
	  MatrixXd X_reaction_,
	  MatrixXd X_individual_,
	  VectorXi cell_,
	  VectorXi positive_,
	  VectorXi time_,
	  int rows_,
	  int cols_,
	  VectorXi internalPoints_,
	  int nTime_,
	  int diffusionType_,
	  bool dirichlet_ = true,
	  double lengthX_ = 1,
	  double lengthY_ = 1) :grid(rows_, cols_, internalPoints_, lengthX_, lengthY_) {

	  //assignments
	  mu_0 = mu_0_;
	  alpha = alpha_;
	  gamma = gamma_;
	  longLat = longLat_;
	  sigma = sigma_;
	  kappa = kappa_;
	  eta = eta_;
	  coords = coords_;
	  X_diffusion = X_diffusion_;
	  X_reaction = X_reaction_;
	  X_individual = X_individual_;
	  nTime = nTime_;
	  diffusionType = diffusionType_;
	  dirichlet = dirichlet_;
	  diffusionCovariates = true;
	  individualCovariates = true;
	  computed = false;
	  cell = cell_;
	  positive = positive_;
	  time = time_;


	  u.setZero(grid.nInternal, nTime); //initialize diffusion

	  //initialize mu and lambda
	  VectorXd muVec, lambdaVec;
	  muVec.setZero(X_reaction.rows());
	  muVec = muVec.array() + mu_0;
	  muVec = muVec.array() * (X_diffusion * alpha).array().exp().array();

	  lambdaVec.setZero(X_reaction.rows());
	  lambdaVec = (X_reaction * gamma).array().exp();

	  Map<MatrixXd> muMap(muVec.data(), grid.nFull, nTime);
	  Map<MatrixXd> lambdaMap(lambdaVec.data(), grid.nFull, nTime);

	  mu = muMap;
	  lambda = lambdaMap;
  }
  
  
  void setInitialConditions(){
    MatrixXd init(grid.rows,grid.cols);
    int count=0;
    for (int cInd=0;cInd<grid.cols;cInd++){
      for (int rInd=0;rInd<grid.rows;rInd++){
        double d1=coords(count,0)-longLat(0);
        double d2=coords(count,1)-longLat(1);
        double dist2=std::pow(d1,2)+std::pow(d2,2);
        init(rInd,cInd)=expit(kappa)*std::exp(-dist2/(2*std::pow(sigma,2)));
        count=count+1;
      }
    }
    
    Map<MatrixXd> u0_map(u.col(0).data(), grid.rows_internal, grid.cols_internal);
    u0_map=init.block(1,1,grid.rows_internal,grid.cols_internal);
  }
  
  void computeDiffusion(int nIter=100,
                        double tol=1.0e-16){
    if (computed) return;
    
    VectorXd rhs(grid.nInternal);
    setInitialConditions();
    for (int i=1;i<nTime;i++){
      
      //get growth coefficients for previous time point
      Map<MatrixXd> lambdaMatrix(lambda.col(i-1).data(),grid.rows,grid.cols);
      MatrixXd internalLambdaMatrix(lambdaMatrix.block(1,1,grid.rows_internal,grid.cols_internal));
      Map<VectorXd> lambda_i_1(internalLambdaMatrix.data(),grid.nInternal);
      
      //Kolmogorov-Fisher
      //U2-U1=HU2+Lambda*U1(1-U1)
      //(I-H)U2=U1+Lambda*U1(1-U1)
      rhs=u.col(i-1).array()+lambda_i_1.array()*u.col(i-1).array()*(1-u.col(i-1).array());
      VectorXd guess(u.col(i-1));

	  //debugging
      u.col(i)=solver(mu.col(i-1),guess,rhs,nIter,tol);
	  //

	  //u.col(i) = u.col(i - 1);
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
  
  
  //derivative of diffusion wrt mu_0/alpha
  MatrixXd du_dmu(int covariateIndex = -1,
	              int nIter=100,
                  double tol=1.0e-16){
    computeDiffusion();
    
    MatrixXd du_dTheta;
    du_dTheta.setZero(grid.nInternal,nTime);
    
    for (int i=0;i<nTime;i++){
      //(du_0/dTheta)=0
      if (i>0){
        
        //get growth coefficients for previous time point
        Map<MatrixXd> lambdaMatrix(lambda.col(i-1).data(),grid.rows,grid.cols);
        MatrixXd internalLambdaMatrix(lambdaMatrix.block(1,1,grid.rows_internal,grid.cols_internal));
        Map<VectorXd> lambda_i_1(internalLambdaMatrix.data(),grid.nInternal);
        
        
        VectorXd rhs(grid.nInternal);
        rhs=u.col(i-1).array()+lambda_i_1.array()*u.col(i-1).array()*(1-u.col(i-1).array());
        VectorXd guess,temp;
        guess.setZero(grid.nInternal);
        temp.setZero(grid.nInternal);
        Ref<VectorXd> mu_i=mu.col(i-1);
        du_dTheta.col(i)=solver(mu_i,guess,rhs,nIter,tol);
        
        temp=du_dTheta.col(i);

		VectorXd dmu_i_dalpha(grid.nFull);
		if (covariateIndex == -1) {
			dmu_i_dalpha.setZero(grid.nFull);
			dmu_i_dalpha = mu_i / mu_0;
		}
		else {
			//get covariate values for previous time point
			Map<MatrixXd> XpMap = time_i_covariate(X_diffusion, i - 1, covariateIndex);
			Map<VectorXd> Xp_i_1(XpMap.data(), grid.nFull);
			dmu_i_dalpha = Xp_i_1.array() * mu_i.array();
		}
		du_dTheta.col(i) = laplacianMultiply(dmu_i_dalpha, temp);
        temp=du_dTheta.col(i);
        du_dTheta.col(i)=solver(mu_i,guess,temp,nIter,tol);

        VectorXd rhs2;
        rhs2.setZero(grid.nInternal);
        rhs2=drhs_dut_1(du_dTheta, lambda_i_1, i);
        du_dTheta.col(i)=du_dTheta.col(i)+solver(mu_i,guess,rhs2,nIter,tol);
      }
    }
    return du_dTheta;
  }
  
  //derivative of diffusion with respect to gamma
  MatrixXd du_dgamma(int covariateIndex=-1,int nIter=100,double tol=1.0e-16){
    computeDiffusion();
    MatrixXd du_dTheta;
    du_dTheta.setZero(grid.nInternal,nTime);
    
    int p=X_reaction.cols();
    
      for (int i=0;i<nTime;i++){
        //(du_0/dTheta)=0
        if (i>0){
          
          //get growth coefficients for previous time point
          Map<MatrixXd> lambdaMatrix(lambda.col(i-1).data(),grid.rows,grid.cols);
          MatrixXd internalLambdaMatrix(lambdaMatrix.block(1,1,grid.rows_internal,grid.cols_internal));
          Map<VectorXd> lambda_i_1(internalLambdaMatrix.data(),grid.nInternal);
          
          //get covariate values for previous time point
          Map<MatrixXd> XpMatrix=time_i_covariate(X_reaction,i-1,covariateIndex);
		  MatrixXd internalXpMatrix(XpMatrix.block(1, 1, grid.rows_internal, grid.cols_internal));
          Map<VectorXd,0,OuterStride<> > Xp_i_1(internalXpMatrix.data(),grid.nInternal, OuterStride<>(internalXpMatrix.outerStride()));
          
          VectorXd rhs2,guess;
          guess.setZero(grid.nInternal);
          
          rhs2.setZero(grid.nInternal);
		  rhs2=drhs_dut_1(du_dTheta, lambda_i_1, i);
		  rhs2 = rhs2.array() + Xp_i_1.array() * lambda_i_1.array() * u.col(i - 1).array() * (1 - u.col(i - 1).array()).array();
          du_dTheta.col(i)=solver(mu.col(i-1),guess,rhs2,nIter,tol);
        }
      }
    return du_dTheta;
  }
  
  std::vector<MatrixXd> du_dlongLat(int nIter=100,
                                    double tol=1.0e-16){
    computeDiffusion();
    
    MatrixXd du_dTheta;
    du_dTheta.setZero(grid.nInternal,nTime);
    
    
    MatrixXd initFull(grid.rows,grid.cols);
    Map<MatrixXd> du_dTheta0_map(du_dTheta.col(0).data(),grid.rows_internal,grid.cols_internal);
    
    
    std::vector<MatrixXd> du_dTheta_list;
    
    for (int llInd=0;llInd<2;llInd++){
      du_dTheta.setZero(grid.nInternal,nTime);
      
      
      for (int i=0;i<nTime;i++){
        if (i==0){
          int count=0;
          for (int cInd=0;cInd<grid.cols;cInd++){
            for (int rInd=0;rInd<grid.rows;rInd++){
              VectorXd d(2);
              d(0)=coords(count,0)-longLat(0);
              d(1)=coords(count,1)-longLat(1);
              double dist2=std::pow(d(0),2)+std::pow(d(1),2);
              initFull(rInd,cInd)=expit(kappa)*std::exp(-dist2/(2*std::pow(sigma,2)))*d(llInd)/std::pow(sigma,2);
              count=count+1;
            }
          }
          du_dTheta0_map=initFull.block(1,1,grid.rows_internal,grid.cols_internal);
        }
        else{
          Map<MatrixXd> lambdaMatrix(lambda.col(i-1).data(),grid.rows,grid.cols);
          MatrixXd internalLambdaMatrix(lambdaMatrix.block(1,1,grid.rows_internal,grid.cols_internal));
          Map<VectorXd> lambda_i_1(internalLambdaMatrix.data(),grid.nInternal);

          
          
          VectorXd rhs2,guess;
          guess.setZero(grid.nInternal);
          
          rhs2.setZero(grid.nInternal);
          rhs2=drhs_dut_1(du_dTheta, lambda_i_1, i);
          du_dTheta.col(i)=solver(mu.col(i-1),guess,rhs2,nIter,tol);
        }
      }
      
      du_dTheta_list.push_back(du_dTheta);
    }
    return du_dTheta_list;
  }
  
  MatrixXd du_dsigma(int nIter=100,
                            double tol=1.0e-16){
    computeDiffusion();
    
    MatrixXd du_dTheta;
    du_dTheta.setZero(grid.nInternal,nTime);
    
    
    MatrixXd initFull(grid.rows,grid.cols);
    Map<MatrixXd> du_dTheta0_map(du_dTheta.col(0).data(),grid.rows_internal,grid.cols_internal);
    
    
    du_dTheta.setZero(grid.nInternal,nTime);
    
    
    for (int i=0;i<nTime;i++){
      if (i==0){
        int count=0;
        for (int cInd=0;cInd<grid.cols;cInd++){
          for (int rInd=0;rInd<grid.rows;rInd++){
            VectorXd d(2);
            d(0)=coords(count,0)-longLat(0);
            d(1)=coords(count,1)-longLat(1);
            double dist2=std::pow(d(0),2)+std::pow(d(1),2);
            initFull(rInd,cInd)=expit(kappa)*std::exp(-dist2/(2*std::pow(sigma,2)))*dist2/std::pow(sigma,3);
            count=count+1;
          }
        }
        du_dTheta0_map=initFull.block(1,1,grid.rows_internal,grid.cols_internal);
      }
      else{
        Map<MatrixXd> lambdaMatrix(lambda.col(i-1).data(),grid.rows,grid.cols);
        MatrixXd internalLambdaMatrix(lambdaMatrix.block(1,1,grid.rows_internal,grid.cols_internal));
        Map<VectorXd> lambda_i_1(internalLambdaMatrix.data(),grid.nInternal);
        
        
        VectorXd rhs2,guess;
        guess.setZero(grid.nInternal);
        
        rhs2.setZero(grid.nInternal);
        rhs2= drhs_dut_1(du_dTheta, lambda_i_1, i);
        du_dTheta.col(i)=solver(mu.col(i-1),guess,rhs2,nIter,tol);
      }
    }
    
    return du_dTheta;
  }

  MatrixXd du_dkappa(int nIter=100,
                            double tol=1.0e-16){
    computeDiffusion();
    
    MatrixXd du_dTheta;
    du_dTheta.setZero(grid.nInternal,nTime);
    
    
    MatrixXd initFull(grid.rows,grid.cols);
    Map<MatrixXd> du_dTheta0_map(du_dTheta.col(0).data(),grid.rows_internal,grid.cols_internal);
    
    
    du_dTheta.setZero(grid.nInternal,nTime);
    
    
    for (int i=0;i<nTime;i++){
      if (i==0){
        int count=0;
        for (int cInd=0;cInd<grid.cols;cInd++){
          for (int rInd=0;rInd<grid.rows;rInd++){
            VectorXd d(2);
            d(0)=coords(count,0)-longLat(0);
            d(1)=coords(count,1)-longLat(1);
            double dist2=std::pow(d(0),2)+std::pow(d(1),2);
            initFull(rInd,cInd)=expit(kappa)*(1-expit(kappa))*std::exp(-dist2/(2*std::pow(sigma,2)));
            count=count+1;
          }
        }
        du_dTheta0_map=initFull.block(1,1,grid.rows_internal,grid.cols_internal);
      }
      else{
        Map<MatrixXd> lambdaMatrix(lambda.col(i-1).data(),grid.rows,grid.cols);
        MatrixXd internalLambdaMatrix(lambdaMatrix.block(1,1,grid.rows_internal,grid.cols_internal));
        Map<VectorXd> lambda_i_1(internalLambdaMatrix.data(),grid.nInternal);
        
        VectorXd rhs2,guess;
        guess.setZero(grid.nInternal);
		rhs2.setZero(grid.nInternal);
		rhs2 = drhs_dut_1(du_dTheta, lambda_i_1, i);
        
        du_dTheta.col(i)=solver(mu.col(i-1),guess,rhs2,nIter,tol);
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
    MatrixXd u_Full(grid.nFull,nTime);
    u_Full.setZero();
    
    for (int i=0;i<nTime;i++){
      Map<MatrixXd> u_i(u_.col(i).data(),grid.rows_internal,grid.cols_internal);
      Map<MatrixXd> u_Full_i(u_Full.col(i).data(),grid.rows,grid.cols);
      u_Full_i.block(1,1,grid.rows_internal,grid.cols_internal)=u_i;
      
      if (!dirichlet){
        for (int rowInd=1;rowInd<grid.rows_internal+1;rowInd++){
          u_Full_i(rowInd,0)=u_Full_i(rowInd,1);
          u_Full_i(rowInd,grid.cols_internal+1)=u_Full_i(rowInd,grid.cols_internal);
        }
        for (int colInd=1;colInd<grid.cols_internal+1;colInd++){
          u_Full_i(0,colInd)=u_Full_i(1,colInd);
          u_Full_i(grid.rows_internal+1,colInd)=u_Full_i(grid.rows_internal,colInd);
        }
      }
    }
    return(u_Full);
  }
  
  VectorXd differentiateLogLikelihood(){
    VectorXd x;
    if (individualCovariates & diffusionCovariates){
      //1 diffusion coefficient, p1 reaction coefficients, p2 individual coefficients
      //2 long lat coefficients, sigma, kappa
		x.setZero(1 + X_diffusion.cols() + X_reaction.cols() + X_individual.cols() + 4);
    }
    else{
      x.setZero(1+X_reaction.cols()+4);
    }
    
    //diffusion
    int count=0;
    Eigen::MatrixXd du_dTheta(grid.nInternal,nTime);
    du_dTheta=du_dmu();
    x(count)=dl_du(du_dTheta);
    count+=1;

	for (int i = 0; i < X_diffusion.cols(); i++) {
		Eigen::MatrixXd du_dTheta(grid.nInternal, nTime);
		du_dTheta = du_dmu(i);
		x(count) = dl_du(du_dTheta);
		count += 1;
	}

    //growth 
    for (int i=0;i<X_reaction.cols();i++){
		Eigen::MatrixXd du_dTheta(grid.nInternal, nTime);
		du_dTheta = du_dgamma(i);
		x(count) = dl_du(du_dTheta);
		count += 1;
    }

    //longLat
	std::vector<MatrixXd> du_dTheta_list;
    du_dTheta_list=du_dlongLat();
    for (int i=0;i<du_dTheta_list.size();i++){
      x(count)=dl_du(du_dTheta_list[i]);
      count+=1;
    }

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
    Map<MatrixXd> X_p(X_.col(p).data(),grid.nFull,nTime);
    Map<MatrixXd> X_p_i(X_p.col(i).data(),grid.rows,grid.cols);
    return X_p_i;
  }

  VectorXd drhs_dut_1(MatrixXd& dut, Map<VectorXd> lambda_i_1, int i) {
	  VectorXd drhs(grid.nInternal);

	  drhs = dut.col(i - 1).array() +
		  lambda_i_1.array() * dut.col(i - 1).array() +
		  -2 * lambda_i_1.array() * u.col(i - 1).array() * dut.col(i - 1).array();
	  return drhs;
  }

  VectorXd laplacianMultiply(const Ref<const Eigen::VectorXd>& mu,Eigen::VectorXd& rhs) {

	  const Map<const MatrixXd> muMatrix(mu.data(), grid.rows,grid.cols);
	  MatrixXd muInternal(muMatrix.block(1, 1, grid.rows_internal, grid.cols_internal));
	  Map<VectorXd> muVec(muInternal.data(), grid.rows_internal * grid.cols_internal);

	  if (diffusionType == 1) {
		  return homogeneous_Lf(muVec.asDiagonal()*rhs,grid,dirichlet);
	  }
	  if (diffusionType == -1) {
		  return muVec.asDiagonal() * homogeneous_Lf(rhs, grid, dirichlet);
	  }
	  return fick_Lf(mu, rhs, grid, dirichlet);
  }
  
  VectorXd solver(const Ref<const VectorXd>& mu_i, VectorXd& guess, VectorXd& rhs, int nIter, double tol){
    return invert(grid, mu_i, guess, rhs, diffusionType, nIter, tol, 1,dirichlet,false);
  }
};


