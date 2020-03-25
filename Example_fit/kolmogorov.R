#params: vector of model parameters
#X: covariate matrix (nSpace*nTime) by nCovariates matrix
#initCov: x and y coordinates of each grid cell, in order from grid cell 1 to grid cell nSpace: nSpace by 2 matrix
#neighborMatrix: adjacency matrix of grid A_{ij}=1 if grid cell i,j neighbors otherwise A_{ij}=0
#nTime is number of time points in simulated diffusion
#KF is just whether or not to use KF
#y length nDeer vector of binary outcomes y_i=1 if ith deer was positive, 0 if negative
#yIndices references ith deer to grid cell s_i and time point t_i, each entry between 1 and (nSpace*nTime)
#tol how precise should the solution to (I-H)^{-1} be, default 10^{-12}, ie U_{t+1}=(I-H)^{-1}U_{t}
#
logLikelihood<-function(params,
                        X,
                        neighborMatrix,
                        geom,
                        nTime,
                        negatives,
                        positives,
                        tol){
  p=dim(X)[2]
  alpha=params[1:p]
  gamma=params[(p+1):(2*p)]

  nSpace=length(geom)
  coords=st_coordinates(st_centroid(geom))
  initParams=params[(2*p+1):length(params)]

  init=NULL
  long=initParams[1]
  lat=initParams[2]
  distances=sqrt((coords[,1]-long)^2+(coords[,2]-lat)^2)
  sigma=exp(initParams[3])
  kappa=1/(1+exp(-initParams[4]))
  init=kappa*exp(-distances^2/(2*sigma^2))
  
  # Eigen::MatrixXd computeDiffusion(Eigen::VectorXd params,
  #                                  Eigen::MatrixXd X,
  #                                  Eigen::VectorXd init,
  #                                  Eigen::SparseMatrix<double> neighborMatrix,
  #                                  int nSpace,
  #                                  int nTime,
  #                                  bool KF,
  #                                  double tol,
  #                                  bool useExplicit=false
  
  KF=TRUE
  U=implicit::computeDiffusion(c(alpha,gamma),
                               X,
                               init,
                               neighborMatrix,
                               nSpace,
                               nTime,
                               KF,
                               tol)
  return(implicit::logLikelihood(negatives,positives,U))
}

#params: vector of model parameters
#X: covariate matrix (nSpace*nTime) by nCovariates matrix
#initCov: x and y coordinates of each grid cell, in order from grid cell 1 to grid cell nSpace: nSpace by 2 matrix
#neighborMatrix: adjacency matrix of grid A_{ij}=1 if grid cell i,j neighbors otherwise A_{ij}=0
#nTime is number of time points in simulated diffusion
#KF is just whether or not to use KF
#y length nDeer vector of binary outcomes y_i=1 if ith deer was positive, 0 if negative
#yIndices references ith deer to grid cell s_i and time point t_i, each entry between 1 and (nSpace*nTime)
#tol how precise should the solution to (I-H)^{-1} be, default 10^{-12}, ie U_{t+1}=(I-H)^{-1}U_{t}
#
computeDiffusion_wrapper<-function(params,
                        X,
                        neighborMatrix,
                        geom,
                        nTime,
                        y_sf,
                        tol){
  p=dim(X)[2]
  alpha=params[1:p]
  gamma=params[(p+1):(2*p)]
  
  nSpace=length(geom)
  coords=st_coordinates(st_centroid(geom))
  initParams=params[(2*p+1):length(params)]
  
  init=NULL
  long=initParams[1]
  lat=initParams[2]
  distances=sqrt((coords[,1]-long)^2+(coords[,2]-lat)^2)
  sigma=exp(initParams[3])
  kappa=1/(1+exp(-initParams[4]))
  init=kappa*exp(-distances^2/(2*sigma^2))
  
  # Eigen::MatrixXd computeDiffusion(Eigen::VectorXd params,
  #                                  Eigen::MatrixXd X,
  #                                  Eigen::VectorXd init,
  #                                  Eigen::SparseMatrix<double> neighborMatrix,
  #                                  int nSpace,
  #                                  int nTime,
  #                                  bool KF,
  #                                  double tol,
  #                                  bool useExplicit=false
  
  KF=TRUE
  U=implicit::computeDiffusion(c(alpha,gamma),
                               X,
                               init,
                               neighborMatrix,
                               nSpace,
                               nTime,
                               KF,
                               tol)
  return(st_as_sf(as.data.frame(U),geom=geom))
}