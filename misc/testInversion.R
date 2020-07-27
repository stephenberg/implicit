library(Matrix)
library(igraph)
library(implicit)
library(microbenchmark)

timeReg=1:4
timeDST=1:4
size=1:4
for (i in 1:4){
  set.seed(33)
  
  rows=round(sqrt(2000*2^i))
  cols=round(sqrt(2000*2^i))
  size[i]=rows*cols
  A=get.adjacency(make_lattice(c(rows,cols)))
  internal=which(colSums(A)==4)
  
  A=A-Matrix::Diagonal(n=dim(A)[1],x=colSums(A))
  
  A=A[internal,][,internal]
  
  b=rnorm((rows-2)*(cols-2))
  d=runif(rows*cols,min=0.001,max=1)
  #d=rep(10,rows*cols)
  # Eigen::VectorXd invert(int rows,
  #                        int cols,
  #                        Eigen::VectorXd& mu,
  #                        Eigen::VectorXd& x,
  #                        Eigen::VectorXd& b,
  #                        int diffusionType,
  #                        int nIter,
  #                        double lengthX=1,
  #                        double lengthY=1)
  t1=implicit::invert(rows = rows-2,
            cols = cols-2,
            mu = d,
            x = rep(0,(rows-2)*(cols-2)),
            b = b,
            diffusionType = 1,
            nIter = 100,
            lengthX = 1,
            lengthY = 1)
  
  I_AD=Matrix::Diagonal(n=(rows-2)*(cols-2),x=1)-A%*%Matrix::Diagonal(n=(rows-2)*(cols-2),x=d[internal])*(rows-1)*(rows-1)
  t2=solve(I_AD,b)
  print(max(abs(t2-t1)))
  
  timeDST[i]=mean(microbenchmark(implicit::invert(rows = rows-2,
                        cols = cols-2,
                        mu = d,
                        x = rep(0,(rows-2)*(cols-2)),
                        b = b,
                        diffusionType = 1,
                        nIter = 100,
                        lengthX = 1,
                        lengthY = 1),times=10)$time)
  timeReg[i]=mean(microbenchmark(solve(1.0*I_AD,b),times=10)$time)
}

plot(log(size),log(timeReg),ylim=c(15,23))
points(log(size),log(timeDST),pch="*")
lm(log(timeReg)~log(size))
lm(log(timeDST)~log(size))
