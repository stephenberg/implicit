# setwd("C:/Dropbox/CWD/cwd_covariates/")
# load("Data/analysisCovariates.RData")
# source("../cwdCode/R/kolmogorov.R")
# library(sp)
# library(ggplot2)
# library(implicit)
# library(spatstat)
# 
# d1=cbind(y_df[y_df$year==1,-(2:3)],longLats[y_df$point[y_df$year==1]]@coords)
# colnames(d1)=c("outcome","x","y")
# d1=as.data.frame(d1)
# #ggplot(data=d1,aes(x=jitter(x,factor=50),y=jitter(y,factor=50),color=outcome,alpha=0.1))+
# #  geom_point(size=1)
# 
# help("ppp")
# 
# fitList=list()
# for (i in 1:13){
#   d1=cbind(y_df[y_df$year==i,-(2:3)],longLats[y_df$point[y_df$year==i]]@coords)
#   colnames(d1)=c("outcome","x","y")
#   d1=as.data.frame(d1)
#   markTest=ppp(x=d1$x,y=d1$y,window=owin(xrange=range(longLats@coords[,1]),yrange=range(longLats@coords[,2])),marks=d1$outcome)
#   fit=Smooth.ppp(markTest,sigma=0.1)
#   fitList[[i]]=fit
#   image(fit$v,zlim=c(0,0.2))
# }
