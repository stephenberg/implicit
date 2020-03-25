library(implicit)
library(Matrix)
library(sf)
setwd("C:/Dropbox/CWD/")

source("wiscScripts/makeAdjacency.R")
source("wiscScripts/makeDesignMatrix.R")
source("implicit/Example_fit/kolmogorov.R")

load("Data/gridCovariates.RData")

#positive and negatives observations
negatives=as.matrix(st_drop_geometry(y_sf)[,1:20])
positives=as.matrix(st_drop_geometry(y_sf)[,21:40])

#intercept only design matrix
X_intercept=makeDesignMatrix(nTime=20,geom=grid)

#design matrix with covariates
timeIndependentCovariates=list(nlcd_sf,soil_sf,river_sf)
timeDependentCovariates=list(culling=culling_sf)
X=makeDesignMatrix(timeIndependentCovariates = timeIndependentCovariates,
                        timeDependentCovariates = timeDependentCovariates)
X=X[,colnames(X)%in%c("intercept",
                 "Pasture_Hay",
                 "Cultivated_Crops",
                 "clayTot",
                 "river",
                 "culling")]

#get adjacency matrix based on the grid geometry
A=makeAdjacency(grid)
boundary=which(colSums(A)<4)
negatives[boundary,]=0
positives[boundary,]=0

#fit model
par0=c(0,-5,500000,000,log(100000),-2)

#argument list of loglikelihood function:
# logLikelihood<-function(params,
#                         X,
#                         neighborMatrix,
#                         geom,
#                         timeSteps,
#                         negatives,
#                         positives,
#                         tol){
logLikeGradient<-function(params,
                          X,
                          neighborMatrix,
                          geom,
                          timeSteps,
                          negatives,
                          positives,
                          tol){
  return(numDeriv::grad(logLikelihood,params,method="Richardson",side=NULL,method.args=list(eps=0.01),
              X,
              neighborMatrix,
              geom,
              timeSteps,
              negatives,
              positives,
              tol))
}

res_IO=optim(par0,
          fn=logLikelihood,
          gr=NULL,
          X_intercept,
          A,
          grid,
          20,
          negatives,
          positives,
          10^-8,
          method="BFGS",
          control=list(fnscale=-1,
                       trace=1,
                       REPORT=1,
                       parscale=c(1,1,10000,10000,1,1),
                       maxit=300)
          )

U_IO_sf=computeDiffusion_wrapper(res_IO$par,
                              X_intercept,
                              A,
                              grid,
                              20,
                              y_sf,
                              10^-8)
plot(U_IO_sf,max.plot=20,border=NA)
parIO=res_IO$par
params=c(parIO[1],0,0,0,0,0,parIO[2],0,0,0,0,0,parIO[3:6])

res_covariates=optim(params,
             fn=logLikelihood,
             gr=NULL,
             X,
             A,
             grid,
             20,
             negatives,
             positives,
             10^-8,
             method="BFGS",
             control=list(fnscale=-1,
                          trace=1,
                          REPORT=1,
                          parscale=c(1,1,1,1,1,0.01,1,1,1,1,1,0.01,10000,10000,1,1),
                          maxit=300)
)

U_covariates_I=numDeriv::jacobian(func=logLikeGradient,
                                  x=res_covariates$par,
                                  method="Richardson",
                                  side=NULL,
                                  method.args=list(),
                                  X,
                                  A,
                                  grid,
                                  20,
                                  negatives,
                                  positives,
                                  10^-8)

U_covariates_sf=computeDiffusion_wrapper(res_covariates$par,
                                 X,
                                 A,
                                 grid,
                                 20,
                                 y_sf,
                                 10^-8)
plot(U_covariates_sf,max.plot=20,border=NA)
par(mfrow=c(1,2))
plot(U_covariates_sf[,17],border=NA,main="Predicted, 2018")
plot(y_sf[,37],max.plot=20,border=NA,main="Actual, 2018")

toMatrix<-function(par,X,covariance=NULL,grid){
  cellSize=st_area(grid[1])
  names=colnames(X)
  p=dim(X)[2]
  alphaInds=1:p
  gammaInds=(p+1):(2*p)
  initInds=(2*p+1):length(par)
  alpha=par[alphaInds]
  gamma=par[gammaInds]
  init=par[initInds]
  init[3]=exp(init[3])
  init[4]=1/(1+exp(-init[4]))
  
  SE=sqrt(diag(covariance))
  SE_alpha=SE[alphaInds]
  SE_gamma=SE[gammaInds]
  SE_init=SE[initInds]
  SE_init[3]=SE[initInds[3]]*init[3]
  SE_init[4]=SE[initInds[4]]*(init[4]*(1-init[4]))
  alpha=cbind(alpha=alpha,SE=SE_alpha)
  gamma=cbind(gamma=gamma,SE=SE_gamma)
  init=cbind(init,SE=SE_init)
  
  
  
  rownames(alpha)=names
  rownames(gamma)=names
  rownames(init)=c("long","lat","sigma","kappa")
  return(list(alpha,gamma,init))
}

speed<-function(X,alpha,gamma,grid){
  cell_area=st_area(grid[1])
  diffusion=matrix(exp(X%*%alpha),ncol=20)*cell_area/10^6
  growth=matrix(exp(X%*%gamma),ncol=20)
  C=matrix(2*sqrt(diffusion*growth),ncol=20)
  C=as.data.frame(C)
  colnames(C)=1:20
  
  return(list(diffusion=st_as_sf(as.data.frame(diffusion),geom=grid),
              growth=st_as_sf(as.data.frame(growth),geom=grid),
              st_as_sf(C,geom=grid)))
}

coefs=toMatrix(res_covariates$par,X,-solve(U_covariates_I),grid)
alpha=coefs[[1]][,1]
gamma=coefs[[2]][,1]
C=speed(X,alpha,gamma,grid)

X_without_culling=X
X[,6]=0
U_withoutCulling_sf=computeDiffusion_wrapper(res_covariates$par,
                                           X_without_culling,
                                           A,
                                           grid,
                                           20,
                                           y_sf,
                                           10^-8)
for (i in 1:10){
  par(mfrow=c(1,2))
  plot(U_covariates_sf[i])
  plot(U_withoutCulling_sf[i])
}
