library(Matrix)
library(igraph)
makeLaplacian<-function(rows,cols){
  A=get.adjacency(make_lattice(c(rows,cols)))
  internal=which(colSums(A)==4)
  
  A=A-Matrix::Diagonal(n=dim(A)[1],x=colSums(A))
  
  A=A[internal,][,internal]
  
}

makeLaplacian2<-function(rows,cols){
  U=matrix(1:(rows*cols),ncol=cols)
  Dx=matrix(0,rows*cols,rows*cols)
  count=1
  for (j in 1:cols){
    for (i in 1:rows){
      if ((j<cols)& (j>1)){
        Dx[count,U[i,j+1]]=1
        Dx[count,U[i,j-1]]=-1
      }
      if (j==cols){
        Dx[count,U[i,j]]=1
        Dx[count,U[i,j-1]]=-1
      }
      if (j==1){
        Dx[count,U[i,j]]=-1
        Dx[count,U[i,j+1]]=1
      }
      count=count+1
    }
  }
  
  Dy=matrix(0,rows*cols,rows*cols)
  count=1
  for (j in 1:cols){
    for (i in 1:rows){
      if ((i<rows)& (i>1)){
        Dy[count,U[i-1,j]]=-1
        Dy[count,U[i+1,j]]=1
      }
      if (i==rows){
        Dy[count,U[i,j]]=1
        Dy[count,U[i-1,j]]=-1
      }
      if (i==1){
        Dy[count,U[i,j]]=-1
        Dy[count,U[i+1,j]]=1
      }
      count=count+1
    }
  }
  
  t1=t(Dx)%*%Dx+t(Dy)%*%Dy
  
  A=get.adjacency(make_lattice(c(rows,cols)))
  internal=which(colSums(A)==4)
  t1=t1[internal,][,internal]
}


rows=5
cols=6
Dx=-diag(cols)
for (i in 2:cols){
  Dx[i,i-1]=1
}
Dy=-diag(rows)
for (i in 2:rows){
  Dy[i-1,i]=1
}
t1=Matrix(-t(kronecker(t(Dx),diag(rows)))%*%kronecker(t(Dx),diag(rows))+
  -t(kronecker(diag(cols),Dy))%*%kronecker(diag(cols),Dy))
A=get.adjacency(make_lattice(c(rows+2,cols+2)))
internal=which(colSums(A)==4)
A=A-Matrix::Diagonal(n=dim(A)[1],x=colSums(A))
t2=A[internal,][,internal]
