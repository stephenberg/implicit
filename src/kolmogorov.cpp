#include <RcppEigen.h>
#include <math.h>
#include "fftPlan.h"

//#include "C:/Users/saber/OneDrive/Documents/R/win-library/3.6/RcppEigen/include/unsupported/Eigen/src/AutoDiff/AutoDiffScalar.h"
//#include "C:/Users/saber/OneDrive/Documents/R/win-library/3.6/RcppEigen/include/unsupported/Eigen/src/AutoDiff/AutoDiffVector.h"
//#include "C:/Users/saber/OneDrive/Documents/R/win-library/3.6/RcppEigen/include/unsupported/Eigen/src/AutoDiff/AutoDiffJacobian.h"


using namespace Rcpp;
// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]



// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// typedef Eigen::VectorXd vec;
// typedef Eigen::MatrixXd mat;
// typedef Eigen::SparseMatrix<double> spMat;

//[[Rcpp::export]]
Eigen::VectorXd solveBiCGSTAB(Eigen::SparseMatrix<double> &lhs,
                              Eigen::VectorXd& rhs,
                              int maxit,
                              double tol,
                              bool print=false){
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
  solver.compute(lhs);
  Eigen::VectorXd x(rhs.size());
  solver.setMaxIterations(maxit);
  solver.setTolerance(tol);
  x=solver.solve(rhs);
  if (print){
    Rcout << "#iterations:     " << solver.iterations() << std::endl;
    Rcout << "estimated error: " << solver.error()      << std::endl;
  }
  return x;
}

Eigen::VectorXd rowSums(Eigen::SparseMatrix<double>& L){
  Eigen::VectorXd rowsums(L.rows());
  rowsums.setZero(L.rows());
  for (int colInd=0; colInd<L.outerSize(); ++colInd){
    for (Eigen::SparseMatrix<double>::InnerIterator it(L,colInd); it; ++it)
    {
      int row=it.row();
      rowsums(row)+=it.value();
    }
  }
  return rowsums;
}


//[[Rcpp::export]]
Eigen::VectorXd solveJacobi(Eigen::SparseMatrix<double>& L,
                            Eigen::VectorXd& rhs,
                            Eigen::VectorXd guess,
                            double tol){
  //function solves (I+diag(rowSums(L))-L)x=rhs for x
  
  Eigen::VectorXd solution,diagVec;
  
  solution.setZero(guess.size());
  diagVec.setZero(guess.size());
  for (int colInd=0; colInd<L.outerSize(); ++colInd){
    for (Eigen::SparseMatrix<double>::InnerIterator it(L,colInd); it; ++it)
    {
      int row=it.row();
      diagVec(row)+=it.value();
    }
  }
  diagVec=diagVec.array()+1.0;
  
  
  solution=guess;
  Eigen::VectorXd oldSolution(solution);
  
  
  
  Eigen::VectorXd Lx;
  Lx.setZero(guess.size());
  Lx=L*solution;
  Eigen::VectorXd errorVec;
  errorVec.setZero(guess.size());
  
  double diff=1;
  int count=0;
  while ((diff>tol)& (count<10000)){
    solution=diagVec.asDiagonal().inverse()*(rhs+Lx);
    Lx=L*solution;
    errorVec=diagVec.asDiagonal()*solution-Lx-rhs;
    diff=(oldSolution-solution).array().abs().array().maxCoeff();
    oldSolution=solution;
    count+=1;
  }
  return solution;
}


// [[Rcpp::export]]
void makeLaplacian(Eigen::SparseMatrix<double>& L,
                   Eigen::VectorXd mu){
  
  for (int colInd=0; colInd<L.outerSize(); ++colInd){
    for (Eigen::SparseMatrix<double>::InnerIterator it(L,colInd); it; ++it)
    {
      int row=it.row();
      int nNonzero=L.col(row).nonZeros();
      if (nNonzero>3){
        if (row>colInd){
          it.valueRef()=mu(colInd);
        }
        else{
          it.valueRef()=mu(row);
        }
      }
      else{
        it.valueRef()=0;
      }
    }
  }
}


//[[Rcpp::export]]
void make_I_Laplacian(Eigen::SparseMatrix<double>& I_L,
                      Eigen::VectorXd mu){
  Eigen::VectorXd diagVec;
  diagVec.setZero(mu.size());
  for (int colInd=0; colInd<I_L.outerSize(); ++colInd){
    for (Eigen::SparseMatrix<double>::InnerIterator it(I_L,colInd); it; ++it)
    {
      int row=it.row();
      int nNonzero=I_L.col(row).nonZeros();
      //the diagonal has been filled in for I_L............
      if (nNonzero>4){
        if (row>colInd){
          it.valueRef()=-mu(colInd);
        }
        else{
          it.valueRef()=-mu(row);
        }
      }
      else{
        it.valueRef()=0;
      }
    }
  }
  for (int colInd=0; colInd<I_L.outerSize(); ++colInd){
    for (Eigen::SparseMatrix<double>::InnerIterator it(I_L,colInd); it; ++it)
    {
      int row=it.row();
      if (row!=colInd){
        diagVec(row)-=it.value();
      }
    }
  }
  diagVec=diagVec.array()+1.0;
  for (int colInd=0;colInd<I_L.outerSize();++colInd){
    for (Eigen::SparseMatrix<double>::InnerIterator it(I_L,colInd); it; ++it)
    {
      int row=it.row();
      if (row==colInd){
        it.valueRef()=diagVec(colInd);
      }
    }
  }
}


// [[Rcpp::export]]
std::vector<Eigen::SparseMatrix<double> > lapTest(Eigen::VectorXd params,
                                                  Eigen::MatrixXd X,
                                                  Eigen::VectorXd init,
                                                  Eigen::SparseMatrix<double> neighborMatrix,
                                                  int nSpace,
                                                  int nTime,
                                                  bool KF,
                                                  double tol,
                                                  bool useExplicit=false){
  int p=X.cols();
  
  Eigen::VectorXd alpha,gamma;
  alpha.setZero(p);
  gamma.setZero(p);
  
  alpha=params.segment(0,p);
  gamma=params.segment(p,p);
  
  int nBoundary=0;
  for (int k=0; k<neighborMatrix.outerSize(); ++k){
    if (neighborMatrix.col(k).nonZeros()<4){
      nBoundary+=1;
    }
  }
  Eigen::VectorXd boundaryPoints;
  boundaryPoints.setZero(nBoundary);
  int count=0;
  for (int k=0; k<neighborMatrix.outerSize(); ++k){
    if (neighborMatrix.col(k).nonZeros()<4){
      boundaryPoints(count)=k;
      count+=1;
    }
  }
  
  Eigen::MatrixXd U;
  U.setZero(nSpace,nTime);
  U.col(0)=init;
  
  for (int k=0;k<boundaryPoints.size();k++){
    U.col(0)(boundaryPoints(k))=0;
  }
  
  Eigen::VectorXd muVec, lambdaVec;
  muVec.setZero(X.rows());
  lambdaVec.setZero(X.rows());
  
  muVec=(X*alpha).array().exp();
  lambdaVec = (X*gamma).array().exp();
  
  Eigen::Map<Eigen::MatrixXd> mu(muVec.data(),nSpace,nTime);
  Eigen::Map<Eigen::MatrixXd> lambda(lambdaVec.data(),nSpace,nTime);
  
  // //contains off diagonals of laplacian, same sparsity structure as neighborMatrix
  Eigen::SparseMatrix<double> L(neighborMatrix);
  
  //set fullLap to false to compare to Jacobi iteration
  ///////////
  Eigen::SparseMatrix<double> I_L(L);
  for (int i=0;i<I_L.rows();i++){
    I_L.insert(i,i)=0;
  }
  //////////
  
  Eigen::VectorXd rhs(nSpace);
  makeLaplacian(L,mu.col(0));
  make_I_Laplacian(I_L,mu.col(0));
  std::vector<Eigen::SparseMatrix<double> > val;
  val.push_back(L);
  val.push_back(I_L);
  return(val);
}


// // [[Rcpp::export]]
// double adTest(double xIn){
//   typedef Eigen::AutoDiffScalar<Eigen::VectorXd> adScalar;
//   typedef Eigen::Matrix<adScalar,Eigen::Dynamic,1> adVec;
//   adVec A_(5, 1);
//   for (int i = 0; i < 5; i++) A_(i).value() = i;
//   const int derivative_num = A_.size();
//   int derivative_idx = 0;
//   for (int i = 0; i < A_.size(); i++) {
//     A_(i).derivatives() = Eigen::VectorXd::Unit(derivative_num, derivative_idx);
//     derivative_idx++;
//   }
//   adScalar A_inv=A_.array().pow(2.0).array().sum();
//   Rcout<<A_inv.derivatives()<<std::endl;
//   //
//   // for (int i = 0; i < A_.rows(); i++) {
//   //   Rcout<<A_(i,1).derivatives()<<std::endl;
//   // }
//   return 0;
// }

//[[Rcpp::export]]
double logLikelihood(Eigen::MatrixXd& negatives,
                     Eigen::MatrixXd& positives,
                     Eigen::MatrixXd& probabilities){
  int nSpace=negatives.rows();
  int nTime=negatives.cols();
  double logLike=0;
  for (int i=0;i<nTime;i++){
    for (int j=0;j<nSpace;j++){
      if (negatives(j,i)>0){
        logLike+=(negatives(j,i)*std::log(1-probabilities(j,i)));
      }
      if (positives(j,i)>0){
        logLike+=(positives(j,i)*std::log(probabilities(j,i)));
      }
    }
  }
  return logLike;
}

//[[Rcpp::export]]
double logLikelihood_general(Eigen::MatrixXd& negatives,
                     Eigen::MatrixXd& positives,
                     Eigen::MatrixXd& probabilities,
                     Eigen::VectorXi& internalPoints,
                     int rows,
                     int cols){
  int nTime=negatives.cols();
  Eigen::Map<Eigen::MatrixXi> internal(internalPoints.data(),rows,cols);
  
  double logLike=0;
  
  for (int i=0;i<nTime;i++){
    Eigen::Map<Eigen::MatrixXd> probs_i(probabilities.col(i).data(),rows+2,cols+2);
    Eigen::Map<Eigen::MatrixXd> pos_i(positives.col(i).data(),rows+2,cols+2);
    Eigen::Map<Eigen::MatrixXd> neg_i(negatives.col(i).data(),rows+2,cols+2);
    for (int colInd=0;colInd<cols;colInd++){
      for (int rowInd=0;rowInd<rows;rowInd++){
        if (internal(rowInd,colInd)!=0){
          logLike+=pos_i(rowInd+1,colInd+1)*std::log(probs_i(rowInd+1,colInd+1));
          logLike+=neg_i(rowInd+1,colInd+1)*std::log(1-probs_i(rowInd+1,colInd+1));
        }
      }
    }
  }
  return logLike;
}

//[[Rcpp::export]]
Eigen::MatrixXd computeDiffusion(Eigen::VectorXd params,
                                 Eigen::MatrixXd X,
                                 Eigen::VectorXd init,
                                 Eigen::SparseMatrix<double> neighborMatrix,
                                 int nSpace,
                                 int nTime,
                                 bool KF,
                                 double tol,
                                 bool useExplicit=false){
  
  int p=X.cols();
  
  Eigen::VectorXd alpha,gamma;
  alpha.setZero(p);
  gamma.setZero(p);
  
  alpha=params.segment(0,p);
  gamma=params.segment(p,p);
  
  int nBoundary=0;
  for (int k=0; k<neighborMatrix.outerSize(); ++k){
    if (neighborMatrix.col(k).nonZeros()<4){
      nBoundary+=1;
    }
  }
  Eigen::VectorXd boundaryPoints;
  boundaryPoints.setZero(nBoundary);
  int count=0;
  for (int k=0; k<neighborMatrix.outerSize(); ++k){
    if (neighborMatrix.col(k).nonZeros()<4){
      boundaryPoints(count)=k;
      count+=1;
    }
  }
  
  Eigen::MatrixXd U;
  U.setZero(nSpace,nTime);
  U.col(0)=init;
  
  for (int k=0;k<boundaryPoints.size();k++){
    U.col(0)(boundaryPoints(k))=0;
  }
  
  Eigen::VectorXd muVec, lambdaVec;
  muVec.setZero(X.rows());
  lambdaVec.setZero(X.rows());
  
  muVec=(X*alpha).array().exp();
  lambdaVec = (X*gamma).array().exp();
  
  Eigen::Map<Eigen::MatrixXd> mu(muVec.data(),nSpace,nTime);
  Eigen::Map<Eigen::MatrixXd> lambda(lambdaVec.data(),nSpace,nTime);
  
  // //contains off diagonals of laplacian, same sparsity structure as neighborMatrix
  Eigen::SparseMatrix<double> L(neighborMatrix);
  
  //set fullLap to false to compare to the Jacobi iteration
  ///////////
  bool fullLap=false;
  Eigen::SparseMatrix<double> I_L(L);
  if (fullLap){
    for (int i=0;i<I_L.rows();i++){
      I_L.insert(i,i)=0;
    }
  }
  //////////
  
  Eigen::VectorXd rhs(nSpace);
  
  for (int i=1;i<nTime;i++){
    //Kolmogorov-Fisher
    if (KF){
      //U2-U1=HU2+Lambda*U1(1-U1)
      //(I-H)U2=U1+Lambda*U1(1-U1)
      rhs=U.col(i-1).array()+lambda.col(i-1).array()*U.col(i-1).array()*(1-U.col(i-1).array());
    }
    //exponential growth
    else{
      //U2-U1=HU2+Lambda*U1
      //(I-H)U2=U1+Lambda*U1
      rhs=U.col(i-1).array()+lambda.col(i-1).array()*U.col(i-1).array();
    }
    if (!fullLap){
      makeLaplacian(L,mu.col(i));
      U.col(i)=solveJacobi(L,rhs,U.col(i-1),tol);
    }
    else{
      make_I_Laplacian(I_L,mu.col(i));
      
      ///////BiCGSTAB
      Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
      solver.setTolerance(tol);
      solver.compute(I_L);
      ///////
      
      ///////Conjugate gradient
      // Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > solver;
      // solver.compute(I_L);
      ///////
      
      ///////Sparse LU
      // Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >   solver;
      // //Compute the ordering permutation vector from the structural pattern of A
      // solver.analyzePattern(I_L);
      // //Compute the numerical factorization
      // solver.factorize(I_L);
      ///////
      
      Eigen::VectorXd guess(nSpace);
      guess=U.col(i-1);
      U.col(i)= solver.solveWithGuess(rhs,guess);
      //change to solve(rhs) for non iterative solvers
    }
    
    // //for comparison only, often unstable
    if (useExplicit){
      Eigen::VectorXd diagVec;
      diagVec.setZero(nSpace);
      for (int colInd=0; colInd<L.outerSize(); ++colInd){
        for (Eigen::SparseMatrix<double>::InnerIterator it(L,colInd); it; ++it)
        {
          int row=it.row();
          diagVec(row)+=it.value();
        }
      }
      U.col(i)=U.col(i-1)+L*U.col(i-1)-diagVec.asDiagonal()*U.col(i-1);
      U.col(i)=U.col(i).array()+lambda.col(i-1).array()*U.col(i-1).array()*(1-U.col(i-1).array());
      for (int j=0;j<nSpace;j++){
        if (U.col(i)(j)<0){
          U.col(i)(j)=0;
        }
        if (U.col(i)(j)>1){
          U.col(i)(j)=1;
        }
      }
    }
    
    
    for (int k=0;k<boundaryPoints.size();k++){
      U.col(i)(boundaryPoints(k))=0;
    }
    if (KF){
      for (int k=0;k<nSpace;k++){
        if (U.col(i)(k)>1){
          U.col(i)(k)=1;
        }
      }
    }
  }
  return U;
}

//discrete sine transformation implementation
//primarily for testing purposes
//[[Rcpp::export]]
Eigen::MatrixXd dst1(Eigen::MatrixXd U){
  int rows=U.rows();
  int cols=U.cols();
  Eigen::MatrixXd V(rows,cols);
  fftPlan p(rows,
            cols,
            U.data(),
            V.data(),
            U.data());
  p.dst_1();
  return(V);
}

//implicit matrix vector multiplication with inhomogeneous diffusion coefficients
//options are 1: diffusion coefficient in Laplacian
//0: diffusion coefficient in middle of Laplacian (Fickian diffusion)
//-1: diffusion coefficient before Laplacian
//[[Rcpp::export]]
Eigen::VectorXd L_f(Eigen::VectorXd& mu_,
                    Eigen::VectorXd& f_,
                    int rows,
                    int cols,
                    int diffusionType,
                    double lengthX=1,
                    double lengthY=1){
  Eigen::VectorXd b_;
  b_.setZero(rows*cols);
  
  double h_y,h_x;
  h_y=lengthY/(rows+1);
  h_x=lengthX/(cols+1);
  
  //use value of diffusion coefficient for boundary cells
  Map<Eigen::MatrixXd> mu(mu_.data(),rows+2,cols+2);
  
  //f and b implicitly defined as 0 on boundary points
  Map<Eigen::MatrixXd> f(f_.data(),rows,cols);
  Map<Eigen::MatrixXd> b(b_.data(),rows,cols);
  
  //diffusion type=1: return \Lambda\cdot\Lambda(\mu f)
  if (diffusionType==1){
    
    //Eigen::Block<Eigen::MatrixXd> internalMu(mu,1,1,rows,cols);
    Ref<Eigen::Map<Eigen::MatrixXd> > internalMu=mu.block(1,1,rows,cols);
    
    //d2dx2
    for (int j=0;j<cols;j++){
      for (int i=0;i<rows;i++){
        if (j==0){
          b(i,j)+=(internalMu(i,j+1)*f(i,j+1)-2*internalMu(i,j)*f(i,j))/std::pow(h_x,2);
        }
        if (j==(cols-1)){
          b(i,j)+=(internalMu(i,j-1)*f(i,j-1)-2*internalMu(i,j)*f(i,j))/std::pow(h_x,2);
        }
        if ((j>0)&(j<(cols-1))){
          b(i,j)+=(internalMu(i,j+1)*f(i,j+1)+
            internalMu(i,j-1)*f(i,j-1)-2*internalMu(i,j)*f(i,j))/std::pow(h_x,2);
        }
      }
    }
    
    //d2dy2
    for (int j=0;j<cols;j++){
      for (int i=0;i<rows;i++){
        if (i==0){
          b(i,j)+=(internalMu(i+1,j)*f(i+1,j)-2*internalMu(i,j)*f(i,j))/std::pow(h_y,2);
        }
        if (i==(rows-1)){
          b(i,j)+=(internalMu(i-1,j)*f(i-1,j)-2*internalMu(i,j)*f(i,j))/std::pow(h_y,2);
        }
        if ((i>0)&(i<(rows-1))){
          b(i,j)+=(internalMu(i+1,j)*f(i+1,j)+
            internalMu(i-1,j)*f(i-1,j)-2*internalMu(i,j)*f(i,j))/std::pow(h_y,2);
        }
      }
    }
  }
  
  //Fickian diffusion
  //diffusion type=0: return \Lambda\cdot(\mu\Lambda f)
  if (diffusionType==0){
    
    //d2dx2
    for (int j=0;j<cols;j++){
      for (int i=0;i<rows;i++){
        double leftMu,rightMu;
        
        //note: mu contains the boundary points, in particular
        //mu(i+1,j+1) corresponds to the point f(i,j)
        leftMu=(mu(i+1,j)+mu(i+1,j+1))/2;
        rightMu=(mu(i+1,j+1)+mu(i+1,j+2))/2;
        if (j==0){
          b(i,j)+=(rightMu*f(i,j+1)-(leftMu+rightMu)*f(i,j))/std::pow(h_x,2);
        }
        if (j==(cols-1)){
          b(i,j)+=(leftMu*f(i,j-1)-(leftMu+rightMu)*f(i,j))/std::pow(h_x,2);
        }
        if ((j>0)&(j<(cols-1))){
          b(i,j)+=(rightMu*f(i,j+1)+leftMu*f(i,j-1)-(leftMu+rightMu)*f(i,j))/std::pow(h_x,2);
        }
      }
    }
    
    //d2ydy2
    for (int j=0;j<cols;j++){
      for (int i=0;i<rows;i++){
        double upperMu,lowerMu;
        upperMu=(mu(i,j+1)+mu(i+1,j+1))/2;
        lowerMu=(mu(i+1,j+1)+mu(i+2,j+1))/2;
        if (i==0){
          b(i,j)+=(lowerMu*f(i+1,j)-(lowerMu+upperMu)*f(i,j))/std::pow(h_y,2);
        }
        if (i==(rows-1)){
          b(i,j)+=(upperMu*f(i-1,j)-(lowerMu+upperMu)*f(i,j))/std::pow(h_y,2);
        }
        if ((i>0)& (i<(rows-1))){
          b(i,j)+=(lowerMu*f(i+1,j)+upperMu*f(i-1,j)-(lowerMu+upperMu)*f(i,j))/std::pow(h_y,2);
        }
      }
    }
  }
  
  //diffusion type=-1: return \mu\Lambda\cdot\Lambda(f)
  if (diffusionType==-1){
    Ref<Eigen::Map<Eigen::MatrixXd> > internalMu=mu.block(1,1,rows,cols);
    
    //d2dx2
    for (int j=0;j<cols;j++){
      for (int i=0;i<rows;i++){
        if (j==0){
          b(i,j)+=internalMu(i,j)*(f(i,j+1)-2*f(i,j))/std::pow(h_x,2);
        }
        if (j==(cols-1)){
          b(i,j)+=internalMu(i,j)*(f(i,j-1)-2*f(i,j))/std::pow(h_x,2);
        }
        if ((j>0)&(j<(cols-1))){
          b(i,j)+=internalMu(i,j)*(f(i,j+1)+f(i,j-1)-2*f(i,j))/std::pow(h_x,2);
        }
      }
    }
    
    //d2dy2
    for (int j=0;j<cols;j++){
      for (int i=0;i<rows;i++){
        if (i==0){
          b(i,j)+=internalMu(i,j)*(f(i+1,j)-2*f(i,j))/std::pow(h_y,2);
        }
        if (i==(rows-1)){
          b(i,j)+=internalMu(i,j)*(f(i-1,j)-2*f(i,j))/std::pow(h_y,2);
        }
        if ((i>0)&(i<(rows-1))){
          b(i,j)+=internalMu(i,j)*(f(i+1,j)+f(i-1,j)-2*f(i,j))/std::pow(h_y,2);
        }
      }
    }
  }
  return b_;
}


//implicit multiplication by the Laplacian with homogeneous diffusion
//[[Rcpp::export]]
Eigen::VectorXd homogeneous_L_f(Eigen::VectorXd& f_,
                                int rows,
                                int cols,
                                double lengthX=1,
                                double lengthY=1){
  Eigen::VectorXd b_;
  b_.setZero(rows*cols);
  
  double h_y,h_x;
  h_y=lengthY/(rows+1);
  h_x=lengthX/(cols+1);
  
  
  //f and b implicitly defined as 0 on boundary points
  Map<Eigen::MatrixXd> f(f_.data(),rows,cols);
  Map<Eigen::MatrixXd> b(b_.data(),rows,cols);
  
  //d2dx2
  for (int j=0;j<cols;j++){
    for (int i=0;i<rows;i++){
      if (j==0){
        b(i,j)+=(f(i,j+1)-2*f(i,j))/std::pow(h_x,2);
      }
      if (j==(cols-1)){
        b(i,j)+=(f(i,j-1)-2*f(i,j))/std::pow(h_x,2);
      }
      if ((j>0)&(j<(cols-1))){
        b(i,j)+=(f(i,j+1)+f(i,j-1)-2*f(i,j))/std::pow(h_x,2);
      }
    }
  }
  
  //d2dy2
  for (int j=0;j<cols;j++){
    for (int i=0;i<rows;i++){
      if (i==0){
        b(i,j)+=(f(i+1,j)-2*f(i,j))/std::pow(h_y,2);
      }
      if (i==(rows-1)){
        b(i,j)+=(f(i-1,j)-2*f(i,j))/std::pow(h_y,2);
      }
      if ((i>0)&(i<(rows-1))){
        b(i,j)+=(f(i+1,j)+f(i-1,j)-2*f(i,j))/std::pow(h_y,2);
      }
    }
  }
  return b_;
}

//implicit multiplication by the Laplacian with homogeneous diffusion
//[[Rcpp::export]]
Eigen::VectorXd irregular_Homogeneous_L_f(Eigen::VectorXd& f_,
                                          int rows,
                                          int cols,
                                          Eigen::VectorXi& boundaryPoints_,
                                          double lengthX=1,
                                          double lengthY=1){
  Eigen::VectorXd b_;
  b_.setZero(rows*cols);
  
  Eigen::Map<MatrixXi> boundary(boundaryPoints_.data(),rows,cols);
  
  double h_y,h_x;
  h_y=lengthY/(rows+1);
  h_x=lengthX/(cols+1);
  
  
  //f and b implicitly defined as 0 on boundary points
  Map<Eigen::MatrixXd> f(f_.data(),rows,cols);
  Map<Eigen::MatrixXd> b(b_.data(),rows,cols);
  
  //d2dx2
  for (int j=0;j<cols;j++){
    for (int i=0;i<rows;i++){
      if (j==0){
        b(i,j)+=boundary(i,j)*(boundary(i,j+1)*f(i,j+1)-2*f(i,j))/std::pow(h_x,2);
      }
      if (j==(cols-1)){
        b(i,j)+=boundary(i,j)*(boundary(i,j-1)*f(i,j-1)-2*f(i,j))/std::pow(h_x,2);
      }
      if ((j>0)&(j<(cols-1))){
        b(i,j)+=boundary(i,j)*(boundary(i,j+1)*f(i,j+1)+boundary(i,j-1)*f(i,j-1)-2*f(i,j))/std::pow(h_x,2);
      }
    }
  }
  
  //d2dy2
  for (int j=0;j<cols;j++){
    for (int i=0;i<rows;i++){
      if (i==0){
        b(i,j)+= boundary(i, j)*(boundary(i+1,j)*f(i+1,j)-2*f(i,j))/std::pow(h_y,2);
      }
      if (i==(rows-1)){
        b(i,j)+= boundary(i, j)*(boundary(i-1,j)*f(i-1,j)-2*f(i,j))/std::pow(h_y,2);
      }
      if ((i>0)&(i<(rows-1))){
        b(i,j)+= boundary(i, j)*(boundary(i+1,j)*f(i+1,j)+boundary(i-1,j)*f(i-1,j)-2*f(i,j))/std::pow(h_y,2);
      }
    }
  }
  return b_;
}


//implicit multiplication by the Laplacian with homogeneous diffusion
//[[Rcpp::export]]
Eigen::VectorXd general_Homogeneous_Lf(Eigen::VectorXd& f_,
                                       int rows,
                                       int cols,
                                       Eigen::VectorXi& boundaryPoints_,
                                       bool dirichlet=true,
                                       double lengthX=1,
                                       double lengthY=1){
  Eigen::VectorXd b_;
  b_.setZero(rows*cols);
  
  Eigen::Map<MatrixXi> boundary(boundaryPoints_.data(),rows,cols);
  
  double h_y,h_x;
  h_y=lengthY/(rows+1);
  h_x=lengthX/(cols+1);
  
  
  //f and b implicitly defined as 0 on boundary points (Dirichlet) or with normal derivative 0 (Neumann)
  Map<Eigen::MatrixXd> f(f_.data(),rows,cols);
  Map<Eigen::MatrixXd> b(b_.data(),rows,cols);
  
  //d2dx2
  for (int j=0;j<cols;j++){
    for (int i=0;i<rows;i++){
      if (j==0){
        //Dirichlet
        if (dirichlet){
          b(i,j)+=boundary(i,j)*(boundary(i,j+1)*f(i,j+1)-2*f(i,j))/std::pow(h_x,2);
        }
        //Neumann
        else{
          b(i,j)+=boundary(i,j)*(boundary(i,j+1)*(f(i,j+1)-f(i,j)))/std::pow(h_x,2);
        }
      }
      if (j==(cols-1)){
        //Dirichlet
        if (dirichlet){
          b(i,j)+=boundary(i,j)*(boundary(i,j-1)*f(i,j-1)-2*f(i,j))/std::pow(h_x,2);
        }
        //Neumann
        else{
          b(i,j)+=boundary(i,j)*(boundary(i,j-1)*(f(i,j-1)-f(i,j)))/std::pow(h_x,2);
        }
      }
      if ((j>0)&(j<(cols-1))){
        //Dirichlet
        if (dirichlet){
          b(i,j)+=boundary(i,j)*(boundary(i,j+1)*f(i,j+1)+boundary(i,j-1)*f(i,j-1)-2*f(i,j))/std::pow(h_x,2);
        }
        //Neumann
        else{
          b(i,j)+=boundary(i,j)*(boundary(i,j+1)*(f(i,j+1)-f(i,j))-boundary(i,j-1)*(f(i,j)-f(i,j-1)))/std::pow(h_x,2);
        }
      }
    }
  }
  
  //d2dy2
  for (int j=0;j<cols;j++){
    for (int i=0;i<rows;i++){
      if (i==0){
        //Dirichlet
        if (dirichlet){
          b(i,j)+= boundary(i, j)*(boundary(i+1,j)*f(i+1,j)-2*f(i,j))/std::pow(h_y,2);
        }
        //Neumann
        else{
          b(i,j)+=boundary(i,j)*(boundary(i+1,j)*(f(i+1,j)-f(i,j)))/std::pow(h_y,2);
        }
      }
      if (i==(rows-1)){
        //Dirichlet
        if (dirichlet){
          b(i,j)+= boundary(i, j)*(boundary(i-1,j)*f(i-1,j)-2*f(i,j))/std::pow(h_y,2);
        }
        else{
          b(i,j)+=boundary(i,j)*(boundary(i-1,j)*(f(i-1,j)-f(i,j)))/std::pow(h_y,2);
        }
      }
      if ((i>0)&(i<(rows-1))){
        //Dirichlet
        if (dirichlet){
          b(i,j)+= boundary(i, j)*(boundary(i+1,j)*f(i+1,j)+boundary(i-1,j)*f(i-1,j)-2*f(i,j))/std::pow(h_y,2);
        }
        //Neumann
        else{
          b(i,j)+=boundary(i,j)*(boundary(i-1,j)*(f(i-1,j)-f(i,j))-boundary(i+1,j)*(f(i,j)-f(i+1,j)))/std::pow(h_y,2);
        }
      }
    }
  }
  return b_;
}

//solve (I-Laplacian)x=b for x, where the Laplacian corresponds to the diffusion
//type as described in L_f; the Laplacian corresponds to a rows by cols grid 
//of *internal points* surrounded implicitly by 
//0's (Dirichlet boundary)
//nIter is the number of iterations, the argument x is an initial guess
//[[Rcpp::export]]
Eigen::VectorXd invert(int rows,
                       int cols,
                       Eigen::VectorXd& mu,
                       Eigen::VectorXd& x,
                       Eigen::VectorXd& b,
                       int diffusionType,
                       int nIter,
                       double lengthX=1,
                       double lengthY=1,
                       int preconditionerType=1,
                       bool debug=false){
  
  int N=rows*cols;
  
  Eigen::VectorXd r(N),p(N),y(N),Ap(N);
  
  
  
  Eigen::Map<MatrixXd> muMatrix(mu.data(),rows+2,cols+2);
  Eigen::MatrixXd internalMuMatrix(muMatrix.block(1,1,rows,cols).array().pow(-1.0));
  Eigen::Map<VectorXd> internalMu_inverse(internalMuMatrix.data(),rows*cols);
  double delta=std::pow(internalMu_inverse.array().sum()/(rows*cols),-1);
  if (debug) Rcout<<"delta="<<delta<<"\n";
  
  fftPlan plan(rows,
               cols,
               r.data(),
               Ap.data(),
               y.data());
  if (diffusionType==1){
    r=internalMu_inverse.asDiagonal()*x
    -homogeneous_L_f(x,rows,cols,lengthX,lengthY)-b;
  }
  if (diffusionType==0){
    r=x-L_f(mu,x,rows,cols,diffusionType,lengthX,lengthY)-b;
  }
  if (diffusionType==-1){
    r=internalMu_inverse.asDiagonal()*x-
      homogeneous_L_f(x,rows,cols,lengthX,lengthY)-
      internalMu_inverse.asDiagonal()*b;
  }
  

  if (preconditionerType==0){
    if (debug) Rcout<<"Preconditioning by L^{-1}"<<"\n";
    plan.multiplyBy_pow_neg_A(-1.0,lengthX,lengthY);
  }
  else{
    if (debug) Rcout<<"Preconditioning by (I-delta\\cdot L)^{-1}"<<"\n";
    plan.multiplyBy_I_A_inverse(lengthX,lengthY,delta);
  }
  
  p=-y;
  
  double r_norm=std::pow(r.array().pow(2).array().sum(),0.5)/r.size();
  double alphaNum=(r.array()*y.array()).array().sum();
  int count=0;
  while ((r_norm>std::pow(10,-10.0))& (count<nIter)){
    if (diffusionType==1){
      Ap=internalMu_inverse.asDiagonal()*p
      -homogeneous_L_f(p,rows,cols,lengthX,lengthY);
    }
    if (diffusionType==0){
      Ap=p-L_f(mu,p,rows,cols,diffusionType,lengthX,lengthY);
    }
    if (diffusionType==-1){
      Ap=internalMu_inverse.asDiagonal()*p
      -homogeneous_L_f(p,rows,cols,lengthX,lengthY);
    }
    double alpha=alphaNum/((p.adjoint()*Ap)(0));
    x=x+alpha*p;
    r=r+alpha*Ap;
    
    if (preconditionerType==0){
      if (debug) Rcout<<"Preconditioning by L^{-1}"<<"\n";
      plan.multiplyBy_pow_neg_A(-1.0,lengthX,lengthY);
    }
    else{
      if (debug) Rcout<<"Preconditioning by (I-delta\\cdot L)^{-1}"<<"\n";
      plan.multiplyBy_I_A_inverse(lengthX,lengthY,delta);
    }
    
    double beta=(r.array()*y.array()).array().sum()/alphaNum;
    p=-y+beta*p;
    alphaNum=alphaNum*beta;
    r_norm=std::pow(r.array().pow(2).array().sum(),0.5)/r.size();
    count+=1;
    
    if (debug){
      Rcout<<"iter "<<count<<"\n";
      Rcout<<"norm "<<r_norm<<"\n\n";
    }
  }
  if (diffusionType==1){
    return internalMu_inverse.asDiagonal()*x;
  }
  return x;
}


//[[Rcpp::export]]
Eigen::VectorXd fixBoundary(Eigen::VectorXd x,
                            Eigen::VectorXi& internalPoints) {
  for (int i = 0; i < x.size(); i++) {
    x(i) = x(i) * internalPoints(i);
  }
  return(x);
}

//solve (I-Laplacian)x=b for x, where the Laplacian corresponds to the diffusion
//type as described in L_f; the Laplacian corresponds to a rows by cols grid 
//of *internal points* surrounded implicitly by 
//0's (Dirichlet boundary)
//nIter is the number of iterations, the argument x is an initial guess
//[[Rcpp::export]]
Eigen::VectorXd invert_Irregular(int rows,
                                 int cols,
                                 Eigen::VectorXd& mu,
                                 Eigen::VectorXd& x,
                                 Eigen::VectorXd& b,
                                 int diffusionType,
                                 double tol,
                                 int nIter,
                                 Eigen::VectorXi& internalPoints,
                                 bool dirichlet=true,
                                 double lengthX=1,
                                 double lengthY=1,
                                 int preconditionerType=1,
                                 bool debug=false){
  int N=rows*cols;
  
  Eigen::VectorXd r(N),p(N),y(N),Ap(N);
  
  
  
  Eigen::Map<MatrixXd> muMatrix(mu.data(),rows+2,cols+2);
  Eigen::MatrixXd internalMuMatrix(muMatrix.block(1,1,rows,cols).array().pow(-1.0));
  Eigen::Map<VectorXd> internalMu_inverse(internalMuMatrix.data(),rows*cols);

  //make sure we start out with properly zero padded mu, x, b
  //the remaining operations will not affect the zero padding?
  internalMu_inverse = fixBoundary(internalMu_inverse,internalPoints);
  double delta=std::pow(internalMu_inverse.array().sum()/(internalPoints.array().sum()),-1);
  if (debug) Rcout<<"delta="<<delta<<"\n";
  
  x = fixBoundary(x,internalPoints);
  b = fixBoundary(b,internalPoints);
  
  fftPlan plan(rows,
               cols,
               r.data(),
               Ap.data(),
               y.data());
  if (diffusionType==1){
    r=internalMu_inverse.asDiagonal()*x
    -general_Homogeneous_Lf(x,rows,cols,internalPoints,dirichlet,lengthX,lengthY)-b;
  }
  if (diffusionType==0){
    r=x-fixBoundary(L_f(mu,x,rows,cols,diffusionType,lengthX,lengthY),internalPoints)-b;
  }
  if (diffusionType==-1){
    r=internalMu_inverse.asDiagonal()*x-
      irregular_Homogeneous_L_f(x,rows,cols,internalPoints,lengthX,lengthY)-
      internalMu_inverse.asDiagonal()*b;
  }
  
  
  if (preconditionerType==0){
    if (debug) Rcout<<"Preconditioning by L^{-1}"<<"\n";
    plan.multiplyBy_pow_neg_A(-1.0,lengthX,lengthY);
  }
  else{
    if (debug) Rcout<<"Preconditioning by (I-delta\\cdot L)^{-1}"<<"\n";
    plan.multiplyBy_I_A_inverse(lengthX,lengthY,delta);
  }
  
  //fixBoundary_reference(y,internalPoints);
  p=-y;
  
  double r_norm=std::pow(r.array().pow(2).array().sum(),0.5)/r.size();
  double alphaNum=(r.array()*y.array()).array().sum();
  int count=0;
  while ((r_norm>tol)& (count<nIter)){
    if (diffusionType==1){
      Ap=internalMu_inverse.asDiagonal()*p
      -general_Homogeneous_Lf(p,rows,cols,internalPoints,dirichlet,lengthX,lengthY);
    }
    if (diffusionType==0){
      Ap=p-fixBoundary(L_f(mu,p,rows,cols,diffusionType,lengthX,lengthY),internalPoints);
    }
    if (diffusionType==-1){
      Ap=internalMu_inverse.asDiagonal()*p
      - irregular_Homogeneous_L_f(p, rows, cols, internalPoints, lengthX, lengthY);
    }
    double alpha=alphaNum/((p.adjoint()*Ap)(0));
    x=x+alpha*p;
    r=r+alpha*Ap;
    
    if (preconditionerType==0){
      if (debug) Rcout<<"Preconditioning by L^{-1}"<<"\n";
      plan.multiplyBy_pow_neg_A(-1.0,lengthX,lengthY);
    }
    else{
      if (debug) Rcout<<"Preconditioning by (I-delta\\cdot L)^{-1}"<<"\n";
      plan.multiplyBy_I_A_inverse(lengthX,lengthY,delta);
    }
    
    double beta=(r.array()*y.array()).array().sum()/alphaNum;
    p=-y+beta*p;
    alphaNum=alphaNum*beta;
    r_norm=std::pow(r.array().pow(2).array().sum(),0.5)/r.size();
    // double r_norm_2=0;
    // for (int i=0;i<r.size();i++){
    //   r_norm_2+=internalPoints(i)*std::pow(r(i),2);
    // }
    // r_norm_2=std::pow(r_norm_2,0.5)/r.size();
    count+=1;
    
    if (debug){
      Rcout<<"iter "<<count<<"\n";
      Rcout<<"norm "<<r_norm<<"\n\n";
    }
  }
  
  if (diffusionType==1){
    return internalMu_inverse.asDiagonal()*x;
  }
  return x;
}


//[[Rcpp::export]]
Eigen::VectorXd irregularLaplacianTest(int rows,
                                       int cols,
                                       Eigen::VectorXd& x,
                                       Eigen::VectorXi& internalPoints,
                                       double lengthX=1,
                                       double lengthY=1){
  return(irregular_Homogeneous_L_f(x,rows,cols,internalPoints,lengthX,lengthY));
}

//[[Rcpp::export]]
Eigen::MatrixXd computeDiffusion2(Eigen::VectorXd params,
                                  Eigen::MatrixXd X,
                                  Eigen::VectorXd init,
                                  int rows,
                                  int cols,
                                  int nTime,
                                  int diffusionType,
                                  double tol,
                                  int nIter,
                                  double lengthX=1,
                                  double lengthY=1,
                                  bool pad=true){
  
  int p=X.cols();
  
  Eigen::VectorXd alpha,gamma;
  alpha.setZero(p);
  gamma.setZero(p);
  
  alpha=params.segment(0,p);
  gamma=params.segment(p,p);
  
  
  //number of internal points, full grid points
  int nInternal=rows*cols;
  int nFull=(rows+2)*(cols+2);
  
  Eigen::MatrixXd U;
  
  U.setZero(nInternal,nTime);
  Eigen::Map<Eigen::MatrixXd> initMap(init.data(), rows + 2, cols + 2);
  Eigen::Map<Eigen::MatrixXd> U0_map(U.col(0).data(), rows, cols);
  U0_map=initMap.block(1,1,rows,cols);
  
  Eigen::VectorXd muVec, lambdaVec;
  muVec.setZero(X.rows());
  lambdaVec.setZero(X.rows());
  
  muVec=(X*alpha).array().exp();
  lambdaVec = (X*gamma).array().exp();
  
  
  Eigen::Map<Eigen::MatrixXd> mu(muVec.data(),nFull,nTime);
  Eigen::Map<Eigen::MatrixXd> lambda(lambdaVec.data(),nFull,nTime);
  
  Eigen::VectorXd rhs(nInternal);
  
  for (int i=1;i<nTime;i++){
    //Kolmogorov-Fisher
    //U2-U1=HU2+Lambda*U1(1-U1)
    //(I-H)U2=U1+Lambda*U1(1-U1)
    rhs=U.col(i-1).array()+lambda.col(i-1).array()*U.col(i-1).array()*(1-U.col(i-1).array());
    Eigen::VectorXd guess(U.col(i-1));
    Eigen::VectorXd mu_i(mu.col(i-1));
    U.col(i)=invert(rows,
          cols,
          mu_i,
          guess,
          rhs,
          diffusionType,
          nIter,
          lengthX,
          lengthY);
  }
  
  if (pad) {
    Eigen::MatrixXd U_full;
    U_full.setZero(nFull, nTime);
    for (int i = 0; i < nTime; i++) {
      Eigen::Map<Eigen::MatrixXd> U_i(U.col(i).data(), rows, cols);
      Eigen::Map<Eigen::MatrixXd> U_full_i(U_full.col(i).data(), rows + 2, cols + 2);
      U_full_i.block(1, 1, rows, cols) = U_i;
    }
    return U_full;
  }
  return U;
}

//[[Rcpp::export]]
Eigen::MatrixXd computeDiffusion_general(Eigen::VectorXd params,
                                         Eigen::MatrixXd X,
                                         Eigen::VectorXd init,
                                         int rows,
                                         int cols,
                                         int nTime,
                                         double tol,
                                         int nIter,
                                         Eigen::VectorXi internalPoints,
                                         bool dirichlet=true,
                                         double lengthX=1,
                                         double lengthY=1,
                                         bool pad=true,
                                         bool iterative=true,
                                         int diffusionType=1){
  int p=X.cols();
  
  Eigen::VectorXd alpha,gamma;
  alpha.setZero(p);
  gamma.setZero(p);
  
  alpha=params.segment(0,p);
  gamma=params.segment(p,p);
  
  
  //number of internal points, full grid points
  int nInternal=rows*cols;
  int nFull=(rows+2)*(cols+2);
  
  Eigen::MatrixXd U;
  
  U.setZero(nInternal,nTime);
  Eigen::Map<Eigen::MatrixXd> initMap(init.data(), rows + 2, cols + 2);
  Eigen::Map<Eigen::MatrixXd> U0_map(U.col(0).data(), rows, cols);
  
  U0_map=initMap.block(1,1,rows,cols);
  U.col(0)=fixBoundary(U.col(0),internalPoints);
  
  Eigen::VectorXd muVec, lambdaVec;
  muVec.setZero(X.rows());
  lambdaVec.setZero(X.rows());
  
  muVec=(X*alpha).array().exp();
  lambdaVec = (X*gamma).array().exp();
  
  
  Eigen::Map<Eigen::MatrixXd> mu(muVec.data(),nFull,nTime);
  Eigen::Map<Eigen::MatrixXd> lambda(lambdaVec.data(),nFull,nTime);
  
  Eigen::VectorXd rhs(nInternal);
  
  for (int i=1;i<nTime;i++){
    //Kolmogorov-Fisher
    //U2-U1=HU2+Lambda*U1(1-U1)
    //(I-H)U2=U1+Lambda*U1(1-U1)
    rhs=U.col(i-1).array()+lambda.col(i-1).array()*U.col(i-1).array()*(1-U.col(i-1).array());
    Eigen::VectorXd guess(U.col(i-1));
    Eigen::VectorXd mu_i(mu.col(i-1));
    U.col(i)=invert_Irregular(rows,
          cols,
          mu_i,
          guess,
          rhs,
          diffusionType,
          tol,
          nIter,
          internalPoints,
          dirichlet,
          lengthX,
          lengthY);
  }
  
  if (pad) {
    Eigen::MatrixXd U_full;
    U_full.setZero(nFull, nTime);
    for (int i = 0; i < nTime; i++) {
      Eigen::Map<Eigen::MatrixXd> U_i(U.col(i).data(), rows, cols);
      Eigen::Map<Eigen::MatrixXd> U_full_i(U_full.col(i).data(), rows + 2, cols + 2);
      U_full_i.block(1, 1, rows, cols) = U_i;
    }
    return U_full;
  }
  return U;
}

//[[Rcpp::export]]
Eigen::MatrixXd computeDiffusion_guess(Eigen::VectorXd params,
                                         Eigen::MatrixXd X,
                                         Eigen::VectorXd init,
                                         int rows,
                                         int cols,
                                         int nTime,
                                         double tol,
                                         int nIter,
                                         Eigen::VectorXi internalPoints,
                                         Eigen::MatrixXd U_guess,
                                         bool dirichlet=true,
                                         double lengthX=1,
                                         double lengthY=1,
                                         bool pad=true,
                                         bool iterative=true,
                                         int diffusionType=1){
  int p=X.cols();
  
  Eigen::VectorXd alpha,gamma;
  alpha.setZero(p);
  gamma.setZero(p);
  
  alpha=params.segment(0,p);
  gamma=params.segment(p,p);
  
  
  //number of internal points, full grid points
  int nInternal=rows*cols;
  int nFull=(rows+2)*(cols+2);
  
  Eigen::MatrixXd U;
  
  U.setZero(nInternal,nTime);
  Eigen::Map<Eigen::MatrixXd> initMap(init.data(), rows + 2, cols + 2);
  Eigen::Map<Eigen::MatrixXd> U0_map(U.col(0).data(), rows, cols);
  
  U0_map=initMap.block(1,1,rows,cols);
  U.col(0)=fixBoundary(U.col(0),internalPoints);
  
  Eigen::VectorXd muVec, lambdaVec;
  muVec.setZero(X.rows());
  lambdaVec.setZero(X.rows());
  
  muVec=(X*alpha).array().exp();
  lambdaVec = (X*gamma).array().exp();
  
  
  Eigen::Map<Eigen::MatrixXd> mu(muVec.data(),nFull,nTime);
  Eigen::Map<Eigen::MatrixXd> lambda(lambdaVec.data(),nFull,nTime);
  
  Eigen::VectorXd rhs(nInternal);
  
  for (int i=1;i<nTime;i++){
    //Kolmogorov-Fisher
    //U2-U1=HU2+Lambda*U1(1-U1)
    //(I-H)U2=U1+Lambda*U1(1-U1)
    rhs=U.col(i-1).array()+lambda.col(i-1).array()*U.col(i-1).array()*(1-U.col(i-1).array());
    
    Eigen::VectorXd guess(rows*cols);
    Eigen::Map<Eigen::MatrixXd> guess_i(guess.data(),rows,cols);
    Eigen::Map<Eigen::MatrixXd> U_guess_i(U_guess.data(),rows+2,cols+2);
    guess_i=U_guess_i.block(1,1,rows,cols);
    
    Eigen::VectorXd mu_i(mu.col(i-1));
    U.col(i)=invert_Irregular(rows,
          cols,
          mu_i,
          guess,
          rhs,
          diffusionType,
          tol,
          nIter,
          internalPoints,
          dirichlet,
          lengthX,
          lengthY);
  }
  
  if (pad) {
    Eigen::MatrixXd U_full;
    U_full.setZero(nFull, nTime);
    for (int i = 0; i < nTime; i++) {
      Eigen::Map<Eigen::MatrixXd> U_i(U.col(i).data(), rows, cols);
      Eigen::Map<Eigen::MatrixXd> U_full_i(U_full.col(i).data(), rows + 2, cols + 2);
      U_full_i.block(1, 1, rows, cols) = U_i;
    }
    return U_full;
  }
  return U;
}