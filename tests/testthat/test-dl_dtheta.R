# 
# 
# test_that("dl_dtheta works", {
# <<<<<<< HEAD
#   load("../../misc/setup.RData")
#   diffusionType=0
# =======
#   source("../../misc/diffusion.R")
#   load("../../misc/setup.RData")
#   
# >>>>>>> 263e79fd2f86d4863e510a160792143bb747d30a
#   unwrap<-function(params,X_reaction,X_individual){
#     mu_0=params[1]
#     p0=ncol(X_reaction)
#     p1=ncol(X_individual)
#     gamma=params[2:(p0+1)]
#     longLat=params[(p0+2):(p0+3)]
#     sigma=params[(p0+4)]
#     kappa=params[p0+5]
#     individual=params[(p0+6):length(params)]
#     names(gamma)=names(X_reaction)
#     names(individual)=names(X_individual)
#     return(list(mu_0=mu_0,
#                 gamma=gamma,
#                 longLat=longLat,
#                 sigma=sigma,
#                 kappa=kappa,
#                 individual=individual))
#   }
#   llwrapper<-function(params){
#     
#     mu_0=params[1]
#     pr=ncol(X_reaction)
#     gamma=params[2:(pr+1)]
#     longLat=params[(pr+2):(pr+3)]
#     sigma=params[pr+4]
#     kappa=params[pr+5]
#     eta=params[(pr+6):length(params)]
#     return(loglikelihood(mu_0,
#       gamma,
#       longLat,
#       sigma,
#       kappa,
#       eta,
#       coords,
#       X_reaction,
#       X_individual,
#       cell,
#       positive,
#       time,
#       rows-2,
#       cols-2,
#       nTime,
#       diffusionType,
#       TRUE,
#       lengthX,
#       lengthY,
#       TRUE))
#   }
#   diffusionwrapper<-function(params,X_reaction,X_individual){
#     par=unwrap(params,X_reaction,X_individual)
#     computeDiffusion(par$mu_0,
#                      par$gamma,
#                      par$longLat,
#                      par$sigma,
#                      par$kappa,
#                      coords,
#                      X_reaction,
#                      rows-2,
#                      cols-2,
#                      nTime,
#                      diffusionType,
#                      TRUE,
#                      lengthX,
#                      lengthY,
#                      TRUE)
#   }
#   derivwrapper<-function(params){
#     mu_0=params[1]
#     pr=ncol(X_reaction)
#     gamma=params[2:(pr+1)]
#     longLat=params[(pr+2):(pr+3)]
#     sigma=params[pr+4]
#     kappa=params[pr+5]
#     eta=params[(pr+6):length(params)]
#     dl_dtheta0=dl_dtheta(mu_0,
#                          gamma,
#                          longLat,
#                          sigma,
#                          kappa,
#                          eta,
#                          coords,
#                          X_reaction,
#                          X_individual,
#                          cell,
#                          positive,
#                          time,
#                          rows-2,
#                          cols-2,
#                          nTime,
#                          diffusionType,
#                          TRUE,
#                          lengthX,
#                          lengthY,
#                          TRUE)
#     return(dl_dtheta0)
#   }
# <<<<<<< HEAD
#   longLat=c(450000,245000)/100000
#   coords=coords/100000
#   sigma=sigma/100000
#   params=c(mu_0,gamma,longLat,sigma,kappa,eta)
#   t1=derivwrapper(params)
#   #t2=numDeriv::grad(llwrapper,params)
#   #H=numDeriv::jacobian(derivwrapper,params,method="Richardson")
# =======
#   
#   H=numDeriv::jacobian(derivwrapper,params,method="Richardson")
# >>>>>>> 263e79fd2f86d4863e510a160792143bb747d30a
#   
#   scales=diag(abs(H))
#   res=optim(par=params,
#             fn=llwrapper,
#             gr=derivwrapper,
#             method="BFGS",
#             control=list(fnscale=-1,
#                          trace=2,
# <<<<<<< HEAD
#                          maxit=1000,
#                          REPORT=5,
# =======
#                          maxit=500,
#                          REPORT=1,
# >>>>>>> 263e79fd2f86d4863e510a160792143bb747d30a
#                          reltol=1e-12))
# 
# 
#   dl_dtheta0=derivWrapper(res$par)/nrow(X_individual)
#   dl_dtheta1=numDeriv::grad(llwrapper,res$par)/nrow(X_individual)
# <<<<<<< HEAD
#   H1=numDeriv::jacobian(derivwrapper,res$par,method="Richardson",method.args=list(eps=10^-6))
#   H2=numDeriv::hessian(llwrapper,res$par)
#   u_init=diffusionwrapper(params,X_reaction,X_individual)
#   u=diffusionwrapper(res$par,X_reaction,X_individual)
#   sds1=sqrt(diag(solve(abs(H1))))
#   sds2=sqrt(diag(solve(abs(H2))))
#   sf_init=data.frame(u_init,geom=grid)%>%st_as_sf()
#   sf1=data.frame(u,geom=grid)%>%st_as_sf()
#   plot(sf_init,max.plot=20,border=NA)
# =======
#   H=numDeriv::jacobian(derivwrapper,res$par,method="Richardson")/nrow(X_individual)
#   u=diffusionwrapper(res$par,X_reaction,X_individual)
#   
#   sf1=data.frame(u,geom=grid)%>%st_as_sf()
# >>>>>>> 263e79fd2f86d4863e510a160792143bb747d30a
#   plot(sf1,max.plot=20,border=NA)
# })