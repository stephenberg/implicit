makeDesignMatrix<-function(timeIndependentCovariates=list(),
                           timeDependentCovariates=list(),
                           nTime=1L,
                           addIntercept=TRUE,
                           geom=NULL){
  
  #type checking of arguments
  if (!(class(timeIndependentCovariates)=="list")) timeIndependentCovariates=list(timeIndependentCovariates)
  if (!(class(timeDependentCovariates)=="list")) timeDependentCovariates=list(timeDependentCovariates)
  
  independentClasses=lapply(timeIndependentCovariates,class)
  dependentClasses=lapply(timeIndependentCovariates,class)
  classCheck<-function(x,class) x!=class
  if (any(unlist(lapply(independentClasses,classCheck,class=c("sf","data.frame"))))){
    stop("Time-independent covariates not sf dataframes.\n Each provided covariate dataframe should have class 'sf','data.frame'")
  }
  if (any(unlist(lapply(independentClasses,classCheck,class=c("sf","data.frame"))))){
    stop("Time-dependent covariates not sf dataframes.\n Each provided covariate dataframe should have class 'sf','data.frame'")
  }
  
  #subtract 1 from each sf length since each sf has a geometry column
  p=sum(unlist(lapply(timeIndependentCovariates,length)))-
    length(timeIndependentCovariates)+
    length(timeDependentCovariates)
  
  
  if (p==0 & (!addIntercept)) stop("No covariates or intercept")
  
  #compare geometries
  geomListInd=lapply(timeIndependentCovariates,st_geometry)
  geomListDep=lapply(timeDependentCovariates,st_geometry)
  equalGeometries=c(geomListInd,geomListDep)%>%
    lapply(function(x) identical(x,c(geomListInd,geomListDep)[[1]]))%>%
    unlist()%>%
    all()
  if (!equalGeometries) stop("Non-identical covariate geometries")
  if ((length(c(geomListInd,geomListDep))>0)) {
    if (!is.null(geom)) warning("geom argument is ignored when covariates provide a geometry")
    geom=st_geometry(c(geomListInd,geomListDep)[[1]])
  }
  if ((length(c(geomListInd,geomListDep))==0)&(is.null(geom))) stop("no geom provided.")
  if (any(class(geom)!=c("sfc_POLYGON","sfc"))) stop("for intercept only models, geom should be an object of class 'sfc_POLYGON' 'sfc' ")
  if (any(unlist(lapply(c(timeIndependentCovariates,timeDependentCovariates),nrow))!=length(geom))) stop("geom length must equal number of covariate rows")
  nSpace=length(geom)
  
  #deal with number of time steps
  nTime=as.integer(nTime)
  if (length(timeDependentCovariates)>0){
    if (!all(unlist(lapply(timeDependentCovariates,ncol))==ncol(timeDependentCovariates[[1]]))) stop("Time dependent covariates must have equal numbers of columns")
    if (nTime!=1){
      if (!all(unlist(lapply(timeDependentCovariates,ncol))==nTime)) warning("nTime argument ignored when time dependent covariate matrix provided")
    }
    nTime=ncol(st_drop_geometry(timeDependentCovariates[[1]]))
  }
  
  if (nTime<1)stop("nTime must be a positive integer")
  
  #make X matrix
  X=NULL
  if (addIntercept) X=cbind(intercept=rep(1,nSpace*nTime),X)
  for (cov_sf in timeIndependentCovariates){
    dataMat=as.matrix(st_drop_geometry(cov_sf))
    tempMat=NULL
    for (j in 1:nTime){
      tempMat=rbind(tempMat,dataMat)
    }
    X=cbind(X,tempMat)
  }
  
  
  tdNames=names(timeDependentCovariates)
  if (length(unique(names(timeDependentCovariates)))<length(timeDependentCovariates)){
    warning("Not all time dependent covariates have unique names: renaming")
    tdNames=paste0("timeDep",1:length(timeDependentCovariates))
  }
  if (length(timeDependentCovariates)>0){
    tempMat=NULL
    for (cov_sf in timeDependentCovariates){
      tempMat=cbind(tempMat,c(as.matrix(st_drop_geometry(cov_sf))))
    }
    colnames(tempMat)=tdNames
    X=cbind(X,as.matrix(tempMat))
  }
  return(X)
}
