setwd("C:/Dropbox/CWD/cwd_covariates/")
load("Data/analysisCovariates.RData")
source("../cwdCode/R/kolmogorov.R")
library(sp)
library(ggplot2)
library(implicit)

params=c(-5,-5,c(-90,43.,log(2),log(10)))

X_intercept=X[,1,drop=FALSE]

KF=TRUE
steps=20
tol=10^(-12)
coords=coordinates(longLats)
neighborMatrix=A


gradFunc<-function(params,X,X_spline,neighborMatrix,steps,KF,y,yIndices,tol,method="simple",kernel){
  val=numDeriv::grad(func=logLikelihood,params,method=method,side=NULL,method.args=list(eps=10^-4),X,X_spline,neighborMatrix,steps,KF,y,yIndices,tol,method,kernel)
  return(val)
}


logLikelihood(params,
              X_intercept,
              coords,
              neighborMatrix,
              steps,
              KF,
              y,
              yIndices,
              tol,
              "simple",
              TRUE)

res_io_KF=optim(par=params,fn=logLikelihood,gr=gradFunc,
             X_intercept,
             coords,
             neighborMatrix,
             steps,
             KF,
             y,
             yIndices,
             tol,
             "simple",
             TRUE,
          method="BFGS",control=list(fnscale=-1,
                                     trace=1,
                                     REPORT=1,
                                     reltol=10^-16,
                                     maxit=200))

U_IO_KF=computeDiffusion_wrapper(res_io_KF$par,
  X_intercept,
  coords,
  neighborMatrix,
  steps,
  KF,
  y,
  yIndices,
  tol,
  "simple",
  TRUE)

H_IO_KF=numDeriv::jacobian(gradFunc,res_io_KF$par,method="simple",side=NULL,list(),
                      X_intercept,
                      coords,
                      neighborMatrix,
                      steps,
                      KF,
                      y,
                      yIndices,
                      tol,
                      "simple",
                      TRUE)

for (i in seq(1,20,by=1)){
  df=data.frame(x=coordinates(harvest_red)[,1],y=coordinates(harvest_red)[,2],probs=U_IO_KF[,i])
  p1=ggplot(df,aes(x=x,y=y,col=probs))+geom_point(size=4,shape="square")+ggtitle(label=round(i))#+scale_color_gradient(limits=c(0,1))
  print(p1)
}



