#include "kf_diffusion.h"

//[[Rcpp::export]]
Eigen::MatrixXd computeDiffusion(double mu_0_,
	                             Eigen::VectorXd alpha_,
                                 Eigen::VectorXd gamma_,
                                 Eigen::VectorXd longLat_,
                                 double sigma_,
                                 double kappa_,
                                 Eigen::MatrixXd coords_,
	                             Eigen::MatrixXd X_diffusion_,
                                 Eigen::MatrixXd X_reaction_,
                                 int rows_,
                                 int cols_,
                                 int nTime_,
                                 int diffusionType_,
                                 bool dirichlet_=true,
	                             double lengthX_=1,
                                 double lengthY_=1,
                                 bool pad=true){
  KF_diffusion diffusion(mu_0_,
	                     alpha_,
                         gamma_,
                         longLat_,
                         sigma_,
                         kappa_,
                         coords_,
	                     X_diffusion_,
                         X_reaction_,
                         rows_,
                         cols_,
                         nTime_,
                         diffusionType_,
                         dirichlet_,
                         lengthX_,
                         lengthY_);
  diffusion.computeDiffusion();
  if (pad){
    Eigen::MatrixXd u((rows_)*(cols_),nTime_);
    u=diffusion.get_u();
    
    Eigen::MatrixXd u_padded((rows_+2)*(cols_+2),nTime_);
    u_padded=diffusion.padDiffusion(u);
    return u_padded;
  }
  return diffusion.get_u();
}

//[[Rcpp::export]]
Eigen::MatrixXd computeDiffusion_NR(double mu_0_,
	Eigen::VectorXd alpha_,
	Eigen::VectorXd gamma_,
	Eigen::VectorXd longLat_,
	double sigma_,
	double kappa_,
	Eigen::MatrixXd coords_,
	Eigen::MatrixXd X_diffusion_,
	Eigen::MatrixXd X_reaction_,
	int rows_,
	int cols_,
	Eigen::VectorXi internalPoints_,
	int nTime_,
	int diffusionType_,
	bool dirichlet_ = true,
	double lengthX_ = 1,
	double lengthY_ = 1,
	bool pad = true) {

	KF_diffusion diffusion(mu_0_,
		alpha_,
		gamma_,
		longLat_,
		sigma_,
		kappa_,
		coords_,
		X_diffusion_,
		X_reaction_,
		rows_,
		cols_,
		internalPoints_,
		nTime_,
		diffusionType_,
		dirichlet_,
		lengthX_,
		lengthY_);
	diffusion.computeDiffusion();
	if (pad) {
		Eigen::MatrixXd u((rows_) * (cols_), nTime_);
		u = diffusion.get_u();

		Eigen::MatrixXd u_padded((rows_ + 2) * (cols_ + 2), nTime_);
		u_padded = diffusion.padDiffusion(u);
		return u_padded;
	}
	return diffusion.get_u();
}


//[[Rcpp::export]]
std::vector<Eigen::MatrixXd> du_dmu(double mu_0_,
	                   Eigen::VectorXd alpha_,
                       Eigen::VectorXd gamma_,
                       Eigen::VectorXd longLat_,
                       double sigma_,
                       double kappa_,
                       Eigen::MatrixXd coords_,
	                   Eigen::MatrixXd X_diffusion_,
                       Eigen::MatrixXd X_reaction_,
                       int rows_,
                       int cols_,
                       int nTime_,
                       int diffusionType_,
                       bool dirichlet_=true,
                       double lengthX_=1,
                       double lengthY_=1,
                       bool pad=true){
  KF_diffusion diffusion(mu_0_,
	                     alpha_,
                         gamma_,
                         longLat_,
                         sigma_,
                         kappa_,
                         coords_,
	                     X_diffusion_,
                         X_reaction_,
                         rows_,
                         cols_,
                         nTime_,
                         diffusionType_,
                         dirichlet_,
                         lengthX_,
                         lengthY_);
  
  std::vector<Eigen::MatrixXd> du_dThetaL;

  Eigen::MatrixXd du_dTheta((rows_) * (cols_), nTime_);
  du_dTheta = diffusion.du_dmu();
  if (pad) {
	  Eigen::MatrixXd u_padded((rows_ + 2) * (cols_ + 2), nTime_);
	  u_padded = diffusion.padDiffusion(du_dTheta);
	  du_dThetaL.push_back(u_padded);
  }
  else {
	  du_dThetaL.push_back(du_dTheta);
  }
  
  for (int i = 0; i < X_diffusion_.cols(); i++) {
	  du_dTheta = diffusion.du_dmu(i);
	  if (pad) {
		  Eigen::MatrixXd u_padded((rows_ + 2) * (cols_ + 2), nTime_);
		  u_padded = diffusion.padDiffusion(du_dTheta);
		  du_dThetaL.push_back(u_padded);
	  }
	  else {
		  du_dThetaL.push_back(du_dTheta);
	  }
  }
  return du_dThetaL;
}

//[[Rcpp::export]]
Eigen::MatrixXd du_dkappa(double mu_0_,
	                      Eigen::VectorXd alpha_,
                          Eigen::VectorXd gamma_,
                          Eigen::VectorXd longLat_,
                          double sigma_,
                          double kappa_,
                          Eigen::MatrixXd coords_,
	                      Eigen::MatrixXd X_diffusion_,
                          Eigen::MatrixXd X_reaction_,
                          int rows_,
                          int cols_,
                          int nTime_,
                          int diffusionType_,
                          bool dirichlet_=true,
                          double lengthX_=1,
                          double lengthY_=1,
                          bool pad=true){
  KF_diffusion diffusion(mu_0_,
	                     alpha_,
                         gamma_,
                         longLat_,
                         sigma_,
                         kappa_,
                         coords_,
	                     X_diffusion_,
                         X_reaction_,
                         rows_,
                         cols_,
                         nTime_,
                         diffusionType_,
                         dirichlet_,
                         lengthX_,
                         lengthY_);
  
  Eigen::MatrixXd du_dTheta((rows_)*(cols_),nTime_);
  du_dTheta=diffusion.du_dkappa();
  
  if (pad){
    
    Eigen::MatrixXd u_padded((rows_+2)*(cols_+2),nTime_);
    u_padded=diffusion.padDiffusion(du_dTheta);
    return u_padded;
  }
  return du_dTheta;
}

//[[Rcpp::export]]
Eigen::MatrixXd du_dsigma(double mu_0_,
	                      Eigen::VectorXd alpha_,
                          Eigen::VectorXd gamma_,
                          Eigen::VectorXd longLat_,
                          double sigma_,
                          double kappa_,
                          Eigen::MatrixXd coords_,
	                      Eigen::MatrixXd X_diffusion_,
                          Eigen::MatrixXd X_reaction_,
                          int rows_,
                          int cols_,
                          int nTime_,
                          int diffusionType_,
                          bool dirichlet_=true,
                          double lengthX_=1,
                          double lengthY_=1,
                          bool pad=true){
  KF_diffusion diffusion(mu_0_,
	                     alpha_,
                         gamma_,
                         longLat_,
                         sigma_,
                         kappa_,
                         coords_,
	                     X_diffusion_,
                         X_reaction_,
                         rows_,
                         cols_,
                         nTime_,
                         diffusionType_,
                         dirichlet_,
                         lengthX_,
                         lengthY_);
  
  Eigen::MatrixXd du_dTheta((rows_)*(cols_),nTime_);
  du_dTheta=diffusion.du_dsigma();
  
  if (pad){
    
    Eigen::MatrixXd u_padded((rows_+2)*(cols_+2),nTime_);
    u_padded=diffusion.padDiffusion(du_dTheta);
    return u_padded;
  }
  return du_dTheta;
}

//[[Rcpp::export]]
std::vector<Eigen::MatrixXd> du_dgamma(double mu_0_,
	                                   Eigen::VectorXd alpha_,
                                       Eigen::VectorXd gamma_,
                                       Eigen::VectorXd longLat_,
                                       double sigma_,
                                       double kappa_,
                                       Eigen::MatrixXd coords_,
	                                   Eigen::MatrixXd X_diffusion_,
                                       Eigen::MatrixXd X_reaction_,
                                       int rows_,
                                       int cols_,
                                       int nTime_,
                                       int diffusionType_,
                                       bool dirichlet_=true,
                                       double lengthX_=1,
                                       double lengthY_=1,
                                       bool pad=true){
  KF_diffusion diffusion(mu_0_,
	                     alpha_,
                         gamma_,
                         longLat_,
                         sigma_,
                         kappa_,
                         coords_,
	                     X_diffusion_,
                         X_reaction_,
                         rows_,
                         cols_,
                         nTime_,
                         diffusionType_,
                         dirichlet_,
                         lengthX_,
                         lengthY_);
  
  std::vector<Eigen::MatrixXd> du_dThetaL;

  Eigen::MatrixXd du_dTheta((rows_) * (cols_), nTime_);
  for (int i = 0; i < X_reaction_.cols(); i++) {
	  du_dTheta = diffusion.du_dgamma(i);
	  if (pad) {
		  Eigen::MatrixXd u_padded((rows_ + 2) * (cols_ + 2), nTime_);
		  u_padded = diffusion.padDiffusion(du_dTheta);
		  du_dThetaL.push_back(u_padded);
	  }
	  else {
		  du_dThetaL.push_back(du_dTheta);
	  }
  }
  return du_dThetaL;
}

//[[Rcpp::export]]
std::vector<Eigen::MatrixXd> du_dlongLat(double mu_0_,
	                                     Eigen::VectorXd alpha_,
                                         Eigen::VectorXd gamma_,
                                         Eigen::VectorXd longLat_,
                                         double sigma_,
                                         double kappa_,
                                         Eigen::MatrixXd coords_,
	                                     Eigen::MatrixXd X_diffusion_,
                                         Eigen::MatrixXd X_reaction_,
                                         int rows_,
                                         int cols_,
                                         int nTime_,
                                         int diffusionType_,
                                         bool dirichlet_=true,
                                         double lengthX_=1,
                                         double lengthY_=1,
                                         bool pad=true){
  KF_diffusion diffusion(mu_0_,
	                     alpha_,
                         gamma_,
                         longLat_,
                         sigma_,
                         kappa_,
                         coords_,
	                     X_diffusion_,
                         X_reaction_,
                         rows_,
                         cols_,
                         nTime_,
                         diffusionType_,
                         dirichlet_,
                         lengthX_,
                         lengthY_);
  
  std::vector<Eigen::MatrixXd> du_dTheta,du_dTheta_padded;
  du_dTheta=diffusion.du_dlongLat();
  
  
  if (pad){
    for (int i=0;i<du_dTheta.size();i++){
      du_dTheta_padded.push_back(diffusion.padDiffusion(du_dTheta[i]));
    }
    return(du_dTheta_padded);
  }
  return du_dTheta;
}

//[[Rcpp::export]]
Eigen::VectorXd dl_dtheta(double mu_0_,
	                      Eigen::VectorXd alpha_,
                          Eigen::VectorXd gamma_,
                          Eigen::VectorXd longLat_,
                          double sigma_,
                          double kappa_,
                          Eigen::VectorXd eta_,
                          Eigen::MatrixXd coords_,
	                      Eigen::MatrixXd X_diffusion_,
                          Eigen::MatrixXd X_reaction_,
                          Eigen::MatrixXd X_individual_,
                          Eigen::VectorXi cell_,
                          Eigen::VectorXi positive_,
                          Eigen::VectorXi time_,
                          int rows_,
                          int cols_,
                          int nTime_,
                          int diffusionType_,
                          bool dirichlet_=true,
                          double lengthX_=1,
                          double lengthY_=1,
                          bool pad=true){
  KF_diffusion diffusion(mu_0_,
	                     alpha_,
                         gamma_,
                         longLat_,
                         sigma_,
                         kappa_,
                         eta_,
                         coords_,
	                     X_diffusion_,
                         X_reaction_,
                         X_individual_,
                         cell_,
                         positive_,
                         time_,
                         rows_,
                         cols_,
                         nTime_,
                         diffusionType_,
                         dirichlet_,
                         lengthX_,
                         lengthY_);
  return diffusion.differentiateLogLikelihood();
}

//[[Rcpp::export]]
double loglikelihood(double mu_0_,
                     Eigen::VectorXd alpha_,
                     Eigen::VectorXd gamma_,
                     Eigen::VectorXd longLat_,
                     double sigma_,
                     double kappa_,
                     Eigen::VectorXd eta_,
                     Eigen::MatrixXd coords_,
                     Eigen::MatrixXd X_diffusion_,
                     Eigen::MatrixXd X_reaction_,
                     Eigen::MatrixXd X_individual_,
                     Eigen::VectorXi cell_,
                     Eigen::VectorXi positive_,
                     Eigen::VectorXi time_,
                     int rows_,
                     int cols_,
                     int nTime_,
                     int diffusionType_,
                     bool dirichlet_=true,
                     double lengthX_=1,
                     double lengthY_=1,
                     bool pad=true){
  KF_diffusion diffusion(mu_0_,
                         alpha_,
                         gamma_,
                         longLat_,
                         sigma_,
                         kappa_,
                         eta_,
                         coords_,
                         X_diffusion_,
                         X_reaction_,
                         X_individual_,
                         cell_,
                         positive_,
                         time_,
                         rows_,
                         cols_,
                         nTime_,
                         diffusionType_,
                         dirichlet_,
                         lengthX_,
                         lengthY_);

  
  return diffusion.logLikelihood();
}

//[[Rcpp::export]]
Rcpp::List getPredictions(double mu_0_,
	Eigen::VectorXd alpha_,
	Eigen::VectorXd gamma_,
	Eigen::VectorXd longLat_,
	double sigma_,
	double kappa_,
	Eigen::VectorXd eta_,
	Eigen::MatrixXd coords_,
	Eigen::MatrixXd X_diffusion_,
	Eigen::MatrixXd X_reaction_,
	Eigen::MatrixXd X_individual_,
	Eigen::VectorXi cell_,
	Eigen::VectorXi positive_,
	Eigen::VectorXi time_,
	int rows_,
	int cols_,
	int nTime_,
	int diffusionType_,
	bool dirichlet_ = true,
	double lengthX_ = 1,
	double lengthY_ = 1,
	bool pad = true) {
	KF_diffusion diffusion(mu_0_,
		alpha_,
		gamma_,
		longLat_,
		sigma_,
		kappa_,
		eta_,
		coords_,
		X_diffusion_,
		X_reaction_,
		X_individual_,
		cell_,
		positive_,
		time_,
		rows_,
		cols_,
		nTime_,
		diffusionType_,
		dirichlet_,
		lengthX_,
		lengthY_);


	diffusion.computeDiffusion();
	VectorXd rawProbabilities(cell_.size());
	VectorXd adjustedProbabilities(cell_.size());
	MatrixXd u(diffusion.get_u());
	VectorXd shifts(X_individual_ * eta_);
	for (int i = 0; i < cell_.size(); i++) {
		double rawProbability_i = u(cell_(i), time_(i));
		double logit_i = logit(rawProbability_i);
		double logitAdjusted_i = logit_i + shifts(i);
		rawProbabilities(i) = rawProbability_i;
		adjustedProbabilities(i) = expit(logitAdjusted_i);
	}
	return Rcpp::List::create(Rcpp::Named("rawProbabilities")=rawProbabilities,
		                      Rcpp::Named("adjustedProbabilities")=adjustedProbabilities);
}
