load("analysisCovariates.RData")
source("kolmogorov.R")
library(sp)
library(ggplot2)

params=c(-4,-4,c(-90,43.,log(0.5),log(10)))

X_intercept=X[,1,drop=FALSE]

KF=TRUE
steps=20
tol=10^-6
coords=coordinates(longLats)
neighborMatrix=A


gradFunc<-function(params,X,X_spline,neighborMatrix,steps,KF,y,yIndices,tol,method="simple",kernel){
  val=numDeriv::grad(func=logLikelihood,params,method=method,side=NULL,method.args=list(eps=10^-4),X,X_spline,neighborMatrix,steps,KF,y,yIndices,tol,method,kernel)
  return(val)
}

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

params=c(-5,-5,c(-90,43.2,log(2),log(0.1)))
res_io_exp=optim(par=params,fn=logLikelihood,gr=gradFunc,
                 X_intercept,
                 coords,
                 neighborMatrix,
                 steps,
                 !KF,
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


params=res_io_KF$par
params=c(params[1],rep(0,dim(X)[2]-1),params[2],rep(0,dim(X)[2]-1),params[3:length(params)])
res_KF=optim(par=params,fn=logLikelihood,gr=gradFunc,X,coords,neighborMatrix,steps,KF,y,yIndices,tol,"simple",TRUE,
             method="BFGS",control=list(fnscale=-1,
                                        trace=1,
                                        REPORT=1,
                                        maxit=200,
                                        reltol=10^-16))

params=res_io_exp$par
params=c(params[1],rep(0,dim(X)[2]-1),params[2],rep(0,dim(X)[2]-1),params[3:length(params)])

res_exp=optim(par=params,fn=logLikelihood,gr=gradFunc,X,coords,neighborMatrix,steps,!KF,y,yIndices,tol,"simple",TRUE,
              method="BFGS",control=list(fnscale=-1,
                                         trace=1,
                                         REPORT=1,
                                         maxit=200,
                                         reltol=10^-16))
 
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
# 
U_IO_exp=computeDiffusion_wrapper(res_io_exp$par,
                                  X_intercept,
                                  coords,
                                  neighborMatrix,
                                  steps,
                                  !KF,
                                  y,
                                  yIndices,
                                  tol,
                                  "simple",
                                  TRUE)
# 
U_cov_KF=computeDiffusion_wrapper(res_KF$par,
                                  X,
                                  coords,
                                  neighborMatrix,
                                  steps,
                                  KF,
                                  y,
                                  yIndices,
                                  tol,
                                  "simple",
                                  TRUE)
# 
U_cov_exp=computeDiffusion_wrapper(res_exp$par,
                                   X,
                                   coords,
                                   neighborMatrix,
                                   steps,
                                   !KF,
                                   y,
                                   yIndices,
                                   tol,
                                   "simple",
                                   TRUE)
# 
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
# 
H_IO_exp=numDeriv::jacobian(gradFunc,res_io_exp$par,method="simple",side=NULL,list(),
                            X_intercept,
                            coords,
                            neighborMatrix,
                            steps,
                            !KF,
                            y,
                            yIndices,
                            tol,
                            "simple",
                            TRUE)
# 
H_cov_exp=numDeriv::jacobian(gradFunc,res_exp$par,method="simple",side=NULL,list(),
                             X,
                             coords,
                             neighborMatrix,
                             steps,
                             !KF,
                             y,
                             yIndices,
                             tol,
                             "simple",
                             TRUE)
# 
H_cov_KF=numDeriv::jacobian(gradFunc,res_KF$par,method="simple",side=NULL,list(),
                            X,
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
for (i in seq(1,20,by=1)){
  df=data.frame(x=coordinates(harvest_red)[,1],y=coordinates(harvest_red)[,2],probs=U_cov_KF[,i])
  p1=ggplot(df,aes(x=x,y=y,col=probs))+geom_point(size=4,shape="square")+ggtitle(label=round(i))#+scale_color_gradient(limits=c(0,1))
  print(p1)
}

for (i in seq(1,20,by=1)){
  df=data.frame(x=coordinates(harvest_red)[,1],y=coordinates(harvest_red)[,2],probs=U_IO_exp[,i])
  p1=ggplot(df,aes(x=x,y=y,col=probs))+geom_point(size=4,shape="square")+ggtitle(label=round(i))#+scale_color_gradient(limits=c(0,1))
  print(p1)
}
for (i in seq(1,20,by=1)){
  df=data.frame(x=coordinates(harvest_red)[,1],y=coordinates(harvest_red)[,2],probs=U_cov_exp[,i])
  p1=ggplot(df,aes(x=x,y=y,col=probs))+geom_point(size=4,shape="square")+ggtitle(label=round(i))#+scale_color_gradient(limits=c(0,1))
  print(p1)
}


save(res_io_KF,res_io_exp,res_KF,res_exp,H_cov_KF,H_cov_exp,H_IO_exp,H_IO_KF,U_cov_exp,U_cov_KF,U_IO_exp,U_IO_KF,X,file="cwdFits.RData")


# load("Data/harvestProportions.RData")
# load("cwdFits.RData")
# awHarvest=matrix(X[,2],ncol=20)
# propHarvest=matrix(X[,3],ncol=20)
# diffusion=X%*%res_KF$par[1:4]
#
# diffusion=matrix(X%*%res_KF$par[1:4],ncol=20)
# growth=matrix(X%*%res_KF$par[5:8],ncol=20)
# mu=exp(diffusion)
# lambda=exp(growth)
# for (i in 1:20){
#   df=data.frame(x=coordinates(coreGrid)[,1],y=coordinates(coreGrid)[,2],probs=U_cov_KF[,i],awHarvest=awHarvest[,i],propHarvest=propHarvest[,i],lambda=lambda[,i],mu=mu[,i],growth=growth[,i],diffusion=diffusion[,i],hay=matrix(X[,4],ncol=20)[,i])
#   p1=ggplot(df,aes(x=x,y=y,col=growth))+geom_point(size=4,shape="square")+ggtitle(label=round(i))#+scale_color_gradient(limits=c(0,1))
#   print(p1)
# }

