logLikelihood<-function(params,X,initCov,neighborMatrix,timeSteps,KF=TRUE,y,yIndices,tol,method=NULL,kernel=TRUE){
  p=dim(X)[2]
  alpha=params[1:p]
  gamma=params[(p+1):(2*p)]

  nSpace=dim(initCov)[1]
  initParams=params[(2*p+1):length(params)]

  init=NULL
  if (!kernel){
    if (KF){
      init=1/(1+exp(-as.numeric(initCov%*%initParams)))
    }
    if (!KF){
      init=exp(as.numeric(initCov%*%initParams))
    }
  }

  if (kernel){
    long=initParams[1]
    lat=initParams[2]
    distances=sqrt((initCov[,1]-long)^2+(initCov[,2]-lat)^2)
    sigma=exp(initParams[3])
    if (KF){
      kappa=exp(initParams[4])
      init=1/(1+kappa*exp(distances^2/sigma^2))
    }
    if (!KF){
      kappa=exp(initParams[4])
      init=exp(-distances^2/sigma^2)
      init=init*kappa
    }
  }

  
  # Eigen::MatrixXd computeDiffusion(Eigen::VectorXd params,
  #                                  Eigen::MatrixXd X,
  #                                  Eigen::VectorXd init,
  #                                  Eigen::SparseMatrix<double> neighborMatrix,
  #                                  int nSpace,
  #                                  int nTime,
  #                                  bool KF,
  #                                  double tol,
  #                                  bool useExplicit=false
  U=implicit::computeDiffusion(c(alpha,gamma),
                               X,
                               init,
                               neighborMatrix,
                               nSpace,
                               timeSteps,
                               KF,
                               tol)
  if (KF){
    probs=U[yIndices]
  }
  if (!KF){
    probs=2*pnorm(U[yIndices])-1
  }
  return(sum(dbinom(x=y,size=1,prob=probs,log=TRUE)))
}

laplaceTest<-function(params,X,initCov,neighborMatrix,timeSteps,KF=TRUE,y,yIndices,tol,method=NULL,kernel=TRUE){
  p=dim(X)[2]
  alpha=params[1:p]
  gamma=params[(p+1):(2*p)]

  nSpace=dim(initCov)[1]
  initParams=params[(2*p+1):length(params)]

  init=NULL
  if (!kernel){
    if (KF){
      init=1/(1+exp(-as.numeric(initCov%*%initParams)))
    }
    if (!KF){
      init=exp(as.numeric(initCov%*%initParams))
    }
  }

  if (kernel){
    long=initParams[1]
    lat=initParams[2]
    distances=sqrt((initCov[,1]-long)^2+(initCov[,2]-lat)^2)
    sigma=exp(initParams[3])
    if (KF){
      kappa=1/(1+exp(-initParams[4]))
      init=exp(-distances^2/sigma^2)
      init=init/max(init)*kappa
    }
    if (!KF){
      kappa=exp(initParams[4])
      init=exp(-distances/sigma^2)
      init=init*kappa
    }
  }

  U=implicit::lapTest(c(alpha,gamma),X,init,neighborMatrix,nSpace,timeSteps,KF,tol)
  return(U)
}


# computeDiffusion<-function(params,X,X_spline,neighborMatrix,timeSteps,logistic=FALSE,timeHomogeneous=FALSE){
#   nSpace=dim(X_spline)[1]
#   nSpline=dim(X_spline)[2]
# 
#   p=ncol(X)
#   alpha=params[1:p]
#   gamma=params[(p+1):(2*p)]
# 
#   splineCoefficients=params[(2*p+1):(2*p+nSpline)]
#   carryingCapacity=1
# 
#   if (logistic){
#     carryingCapacity=1/(1+exp(-params[length(params)]))
#   }
# 
#   columnSums=Matrix::colSums(neighborMatrix)
#   internalPoints=which(columnSums==4)
#   boundaryPoints=which(columnSums<4)
# 
# 
#   U=matrix(0,nSpace,timeSteps)
#   U[,1]=1/(1+exp(as.numeric(-X_spline%*%splineCoefficients)))
#   U[boundaryPoints,1]=0
# 
#   mu = exp(X[,]%*%as.matrix(alpha))
#   lambda = (X[,]%*%as.matrix(gamma))
# 
#   mu=matrix(mu,ncol=timeSteps)
#   lambda=exp(matrix(lambda,ncol=timeSteps))
# 
#   diffOperator=NULL
#   I_lap=NULL
#   chol_Ilap=NULL
#   if (timeHomogeneous){
#     laplacian=laplaceMatrix(neighborMatrix=neighborMatrix,mu=mu[,1])
#     I_lap=Matrix::Diagonal(n=dim(neighborMatrix)[1],x=1)-laplacian
#   }
#   for (i in 1:(timeSteps-1)){
# 
#     explicit=FALSE
#     if (explicit){
#       laplacian=laplaceMatrix(neighborMatrix=neighborMatrix,mu=mu[,i])
#       U[,i+1]=as.numeric(laplacian%*%U[,i]+lambda[,i]*U[,i]*(1-U[,i]/carryingCapacity))
#     }
#     semiImplicit=TRUE
# 
#     if (semiImplicit){
#       rhs=U[,i]+lambda[,i]*U[,i]*(1-U[,i]/carryingCapacity)
#       UNew=NULL
#       if (timeHomogeneous){
#         UNew=tryCatch(Matrix::solve(I_lap,rhs),error=function(e) return(TRUE))
#         if (is.logical(UNew)){
#           return(-Inf)
#         }
#         U[,i+1]=as.numeric(UNew)
#         U[boundaryPoints,i+1]=0
#       }
#       if (!timeHomogeneous){
#         laplacian=laplaceMatrix(neighborMatrix=neighborMatrix,mu=mu[,i])
#         I_lap=Matrix::Diagonal(n=dim(neighborMatrix)[1],x=1)-laplacian
#         UNew=tryCatch(Matrix::solve(I_lap,rhs),error=function(e) return(TRUE))
#         if (is.logical(UNew)){
#           return(-Inf)
#         }
#         U[,i+1]=as.numeric(UNew)
#         U[boundaryPoints,i+1]=0
#         U[,i+1][U[,i+1]>1]=1
#       }
#     }
#   }
#   return(U)
# }
# 
# 
# laplaceMatrix<-function(neighborMatrix,mu){
# 
#   diffMat=neighborMatrix
#   nNeighbors=neighborMatrix@p[2:(dim(neighborMatrix)[1]+1)]-neighborMatrix@p[1:dim(neighborMatrix)[1]]
#   diagVec=rep(0,dim(neighborMatrix)[1])
#   diffVec=diffMat@x
# 
#   for (i in 1:ncol(neighborMatrix)){
# 
#     neighborStart=neighborMatrix@p[i]+1
#     neighborEnd=neighborMatrix@p[i+1]
#     neighbors=diffMat@i[neighborStart:neighborEnd]+1
# 
#     for (j in 1:nNeighbors[i]){
#       if (nNeighbors[neighbors[j]]>3){
#         if (neighbors[j]<i){
#           diffVec[neighborStart+(j-1)]=mu[neighbors[j]]
#           diagVec[neighbors[j]]=diagVec[neighbors[j]]-mu[neighbors[j]]
#         }
#         if (neighbors[j]>i){
#           diffVec[neighborStart+(j-1)]=mu[i]
#           diagVec[neighbors[j]]=diagVec[neighbors[j]]-mu[i]
#         }
#       }
#     }
#   }
#   diffMat@x=diffVec
#   diagMat=Matrix::Diagonal(n=length(diagVec),x=diagVec)
#   val=diffMat+diagMat
#   return(mat)
# }

computeDiffusion_wrapper<-function(params,X,initCov,neighborMatrix,timeSteps,KF=TRUE,y,yIndices,tol,method=NULL,kernel=TRUE,explicit=FALSE){
  p=dim(X)[2]
  alpha=params[1:p]
  gamma=params[(p+1):(2*p)]

  nSpace=dim(initCov)[1]
  initParams=params[(2*p+1):length(params)]

  init=NULL
  if (!kernel){
    if (KF){
      init=1/(1+exp(-as.numeric(initCov%*%initParams)))
    }
    if (!KF){
      init=exp(as.numeric(initCov%*%initParams))
    }
  }

  if (kernel){
    long=initParams[1]
    lat=initParams[2]
    distances=sqrt((initCov[,1]-long)^2+(initCov[,2]-lat)^2)
    sigma=exp(initParams[3])
    if (KF){
      kappa=exp(initParams[4])
      init=1/(1+kappa*exp(distances^2/sigma^2))
    }
    if (!KF){
      kappa=exp(initParams[4])
      init=exp(-distances^2/sigma^2)
      init=init*kappa
    }
  }
  U=implicit::computeDiffusion(c(alpha,gamma),X,init,neighborMatrix,nSpace,timeSteps,KF,tol,explicit)
  return(U)
}
