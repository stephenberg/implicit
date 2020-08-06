library(sf)
setwd("E:/implicit")
source("misc/makeDesignMatrix.R")

load("misc/covariates_withHarvest.RData")

nTime=20
coords=grid%>%st_centroid()%>%st_coordinates()


#deal with deer on boundary
internal=c(matrix(1:(rows*cols),ncol=cols)[2:(rows-1),][,2:(cols-1)])
y_internal=y_sf[y_sf$Cell%in%internal,]

#0 based indexing for internal points
indices=matrix(1:(rows*cols),ncol=cols)
indices[2:(rows-1),][,2:(cols-1)]=matrix(1:((rows-2)*(cols-2)),ncol=cols-2)
cell=c(indices)[y_internal$Cell]-1


#intercept only design matrix
X_intercept=makeDesignMatrix(nTime=20,geom=grid)

#design matrix with covariates
timeIndependentCovariates=list(nlcd_sf,river_sf)

culling_sf=cbind(culling_sf%>%st_drop_geometry(),matrix(0,1600,14),geom=grid)%>%st_as_sf()
colnames(culling_sf)[1:20]=paste0("year",2002:2021)

landowner_sf=cbind(landowner_sf%>%st_drop_geometry(),matrix(0,1600,14),geom=grid)%>%st_as_sf()
colnames(landowner_sf)[1:20]=paste0("year",2002:2021)


harvest18=harvest_sf[,18]%>%st_drop_geometry()%>%as.matrix()
harvest_sf=cbind(harvest_sf%>%st_drop_geometry(),matrix(rep(harvest18,2),1600,2),geom=grid)%>%st_as_sf()
colnames(harvest_sf)[1:20]=paste0("year",2002:2021)

timeDependentCovariates=list(culling=culling_sf,landowner=landowner_sf,harvest=harvest_sf)
X_reaction=makeDesignMatrix(timeIndependentCovariates = timeIndependentCovariates,
                   timeDependentCovariates = timeDependentCovariates)
X_reaction=X_reaction[,colnames(X_reaction)%in%c("intercept",
                      "Pasture_Hay",
                      "Cultivated_Crops",
                      "clayTot",
                      "river",
                      "culling",
                      "landowner",
                      "harvest")]

X_diffusion=X_reaction[,2:ncol(X_reaction)]
alpha=runif(ncol(X_diffusion))

X_individual=cbind(y_sf$Age,as.integer(y_sf$Sex=="M"))
colnames(X_individual)=c("Age","Sex")


mu_0=0.001875
gamma=c(-2,runif(ncol(X_reaction)-1))
longLat=c(535000,295000)
sigma=10000
kappa=-2.5
eta=rep(0,ncol(X_individual))


time=y_internal$Year-2002 #0 based indexing
positive=as.integer(y_internal$Positive=="Y")




diffusionType=0
lengthX=1
lengthY=1

save(mu_0,
     alpha,
     gamma,
     longLat,
     sigma,
     kappa,
     eta,
     rows,
     cols,
     grid,
     nTime,
     coords,
     X_diffusion,
     X_reaction,
     X_individual,
     cell,
     time,
     positive,
     diffusionType,
     lengthX,
     lengthY,
     file="misc/setup.RData")
