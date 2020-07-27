diffusion<-function(params,
         X,
         rows,
         cols,
         geom,
         lengthX=1,
         lengthY=1,
         nTime){
  
  p=dim(X)[2]
  mu=params[1]
  gamma=params[(2):(p+1)]
  
  coords=st_coordinates(st_centroid(geom))
  initParams=params[(p+2):length(params)]
  
  init=NULL
  long=initParams[1]
  lat=initParams[2]
  distances=sqrt((coords[,1]-long)^2+(coords[,2]-lat)^2)
  sigma=initParams[3]
  scale=1/(1+exp(-initParams[4]))
  init=scale*exp(-distances^2/(2*sigma^2))
  
  U=matrix(0,rows*cols,nTime)
  U[,1]=c(matrix(init,ncol=cols+2)[2:(rows+1),][,2:(cols+1)])

  
  #make laplacian
  hx=lengthX/(cols+1)
  hy=lengthY/(rows+1)
  mx=matrix(0,cols,cols)
  diag(mx)=-2/hx^2
  for (i in 1:(cols-1)){
    mx[i,i+1]=1/hx^2
    mx[i+1,i]=1/hx^2
  }
  mx=Matrix(mx,sparse=TRUE)
  I_y=Matrix::Diagonal(rows,x=1)
  
  my=matrix(0,rows,rows)
  diag(my)=-2/hy^2
  for (i in 1:(rows-1)){
    my[i,i+1]=1/hy^2
    my[i+1,i]=1/hy^2
  }
  I_x=Matrix::Diagonal(cols,x=1)
  
  L=kronecker(mx,I_y)+kronecker(I_x,my)
  
  H=Matrix::kronecker(I_x,I_y)-mu*L
  
  #growth coefficients
  lambda=exp(X%*%gamma)
  
  for (i in 2:nTime){
    lambda_i=matrix(lambda,ncol=nTime)[,i-1]
    lambda_i=c(matrix(lambda_i,ncol=cols+2)[2:(rows+1),][,2:(cols+1)])
    rhs=U[,i-1]+lambda_i*U[,i-1]*(1-U[,i-1])
    U[,i]=as.numeric(solve(H,rhs))
  }
  Ufull=matrix(0,(rows+2)*(cols+2),nTime)
  for (i in 1:nTime){
    temp=matrix(0,rows+2,cols+2)
    temp[2:(rows+1),][,2:(cols+1)]=matrix(U[,i],ncol=cols)
    Ufull[,i]=c(temp)
  }
  return(Ufull)
}