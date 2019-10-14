#include <RcppEigen.h>
#include "C:/Users/saber/OneDrive/Documents/R/win-library/3.6/RcppEigen/include/unsupported/Eigen/src/AutoDiff/AutoDiffScalar.h"
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


//[[Rcpp::export]]
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


// [[Rcpp::export]]
double adTest(double xIn){
  typedef Eigen::AutoDiffScalar<Eigen::VectorXd> adScalar;
  typedef Eigen::Matrix<adScalar,Eigen::Dynamic,1> adVec;
  adVec A_(5, 1);
  for (int i = 0; i < 5; i++) A_(i).value() = i;
  const int derivative_num = A_.size();
  int derivative_idx = 0;
  for (int i = 0; i < A_.size(); i++) {
    A_(i).derivatives() = Eigen::VectorXd::Unit(derivative_num, derivative_idx);
    derivative_idx++;
  }
  adScalar A_inv=A_.array().pow(2.0).array().sum();
  Rcout<<A_inv.derivatives()<<std::endl;
  //
  // for (int i = 0; i < A_.rows(); i++) {
  //   Rcout<<A_(i,1).derivatives()<<std::endl;
  // }
  return 0;
}


// [[Rcpp::export]]
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
