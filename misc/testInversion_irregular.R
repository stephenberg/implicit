library(Matrix)
library(igraph)
library(implicit)
library(microbenchmark)

timeReg=1:4
timeDST=1:4
size=1:4
for (i in 1:4){
  set.seed(33)
  
  rows=round(sqrt(5000*2^i))
  cols=round(sqrt(5000*2^i))
  
  size[i]=rows*cols
  A=get.adjacency(make_lattice(c(rows,cols)))
  
  
  internalFull=(colSums(A)==4)
  A=A-Matrix::Diagonal(n=dim(A)[1],x=Matrix::colSums(A))
  
  intMat=matrix(1:sum(internalFull),rows-2,cols-2)
  internalSampled=unique(c(intMat[round(0.1*rows):round(0.9*rows),][,round(0.1*cols):round(0.9*cols)]))
  #internalSampled=1:((rows-2)*(cols-2))
  
  internalSampled=internalSampled[order(internalSampled)]
  notSampled=which(!((1:((rows-2)*(cols-2)))%in%internalSampled))
  internalMatrix=matrix(0,rows-2,cols-2)
  internalMatrix[internalSampled]=1
  image(Matrix(c(internalMatrix),ncol=cols-2))
  
  
  b=rnorm((rows-2)*(cols-2))
  b[notSampled]=0
  d=runif(rows*cols,min=0.001,max=1)
  d=rep(1,rows*cols)
  
  internalPoints=rep(0,(rows-2)*(cols-2))
  internalPoints[internalSampled]=1
  
  A=A[which(internalFull)[internalSampled],][,which(internalFull)[internalSampled]]
  
  
  
  t1=implicit::invert_Irregular(rows = rows-2,
                                cols = cols-2,
                                mu = d,
                                x = rep(0,(rows-2)*(cols-2)),
                                b = b,
                                diffusionType = 1,
                                tol=10^-10,
                                nIter = 100,
                                internalPoints=internalPoints,
                                dirichlet = TRUE,
                                lengthX = 1,
                                lengthY = 1,
                                debug=TRUE)
  
  identityMatrix=Matrix::Diagonal(n=dim(A)[1],x=1)
  AD=A%*%Matrix::Diagonal(n=dim(A)[1],x=d[which(internalFull)][internalSampled])*(rows-1)^2
  
  I_AD=identityMatrix-AD
  t1_0=t1[internalSampled]
  t2=as.numeric(Matrix::solve(I_AD,b[internalSampled]))
  print(max(abs(t1_0-t2)))
  
  
  timeDST[i]=mean(microbenchmark(implicit::invert_Irregular(rows = rows-2,
                                                            cols = cols-2,
                                                            mu = d,
                                                            x = rep(0,(rows-2)*(cols-2)),
                                                            b = b,
                                                            diffusionType = 1,
                                                            tol=10^-10,
                                                            nIter = 100,
                                                            internalPoints=internalPoints,
                                                            lengthX = 1,
                                                            lengthY = 1)
                                 ,times=10)$time)
  timeReg[i]=mean(microbenchmark(solve(1.0*I_AD,b[internalSampled]),times=10)$time)
}


plot(log(size),log(timeReg),ylim=c(15,23))
points(log(size),log(timeDST),pch="*")

lm(log(timeReg)~log(size))
lm(log(timeDST)~log(size))