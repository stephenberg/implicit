rm(list=ls())
setwd("E:/Dropbox/CWD/wiscScripts/modelFitting/")
load("E:/implicit/misc/setup.RData")
load("cwdFits.RData")
library(dplyr)
library(implicit)
library(reshape2)
library(ggplot2)
getPredictions_wrapper<-function(params){
  par=unwrap(params)
  return(implicit::getPredictions(par$mu_0,
                   par$alpha,
                   par$gamma,
                   par$longLat,
                   par$sigma,
                   par$kappa,
                   par$individual,
                   coords,
                   X_diffusion,
                   X_reaction,
                   X_individual,
                   cell,
                   positive,
                   time,
                   rows,
                   cols,
                   nTime,
                   diffusionType,
                   TRUE,
                   lengthX,
                   lengthY,
                   TRUE))
}

preds=getPredictions_wrapper(modelList[[3]]$res$par)

df=data.frame(raw=preds$rawProbabilities,adjusted=preds$adjustedProbabilities)
df=cbind(df,positive=positive,sex=1*(y_internal$Sex=="M"),age=y_internal$Age,time=time+2002)


#plots of prevalence by age and sex
df_sa=df%>%group_by(sex,age)%>%summarise(positive=mean(positive),prediction_adj=mean(adjusted),prediction_raw=mean(raw))%>%ungroup()
df_sa$sex[df_sa$sex==1]="M"
df_sa$sex[df_sa$sex=="0"]="F"
df_sa$sex=as.factor(df_sa$sex)
ggplot(data=df_sa,aes(x=age,y=positive,col=sex))+geom_point()
ggplot(data=df_sa,aes(x=age,y=prediction_adj,col=sex))+geom_point()
ggplot(data=df_sa,aes(x=age,y=prediction_raw,col=sex))+geom_point()

#pattern over time

#spatial plots
internal=c(matrix(1:(rows*cols),ncol=cols)[2:(rows-1),][,2:(cols-1)])

grid_internal=grid[internal]
df=cbind(df,cell=cell+1) #go from 0 to 1 based indexing

df_temp=df%>%group_by(cell)%>%summarise(sex=mean(sex),age=mean(age),positive=mean(positive),prediction_adj=mean(adjusted),prediction_raw=mean(raw))%>%ungroup()
df_spatial=cbind(df_temp,geom=grid_internal[df_temp$cell])%>%st_as_sf()
ggplot(data=df_spatial,aes(fill=sex))+geom_sf()+scale_fill_gradient(limits=c(0,1))
ggplot(data=df_spatial,aes(fill=age))+geom_sf()
ggplot(data=df_spatial,aes(fill=positive))+geom_sf()+scale_fill_gradient(limits=c(0,0.3))
ggplot(data=df_spatial,aes(fill=prediction_adj))+geom_sf()+scale_fill_gradient(limits=c(0,0.3))

df_time_int=df%>%group_by(time)%>%summarise(observed=mean(positive),fitted=mean(adjusted))%>%ungroup()
df_time=df_time_int%>%melt(measure.vars=c("observed","fitted"),value.name = "prevalence")
ggplot(data=df_time,aes(x=time,y=prevalence,col=variable))+geom_point()

       