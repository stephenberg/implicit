#include <RcppEigen.h>
#include "AutoDiffScalar.h"

using namespace Rcpp;
typedef Eigen::AutoDiffScalar<Eigen::VectorXd> adScalar;
typedef Eigen::Matrix<adScalar,Eigen::Dynamic,Eigen::Dynamic> adMat;
typedef Eigen::SparseMatrix<adScalar> spAdmat;


adScalar adZero(int nDer){
  Eigen::VectorXd zeroVec(nDer);
  zeroVec.setZero(nDer);
  adScalar temp;
  temp.value()=0;
  temp.derivatives()=zeroVec;
  return temp;
}
void setZero_ad(adMat &A, int nDer){
  for (int i = 0; i < A.rows(); i++){
    for (int j=0;j<A.cols();j++){
      A(i,j)=adZero(nDer);
    }
  }
}
// 
adMat solveJacobi(Eigen::SparseMatrix<adScalar>& L,
                  adMat& rhs,
                  adMat guess,
                  double tol){
  //function solves (I+diag(rowSums(L))-L)x=rhs for x
  int nDer=rhs(0,0).derivatives().size();
  adMat solution(guess.rows(),1),diagVec(guess.rows(),1);
  setZero_ad(solution,nDer);
  setZero_ad(diagVec,nDer);
  for (int colInd=0; colInd<L.outerSize(); ++colInd){
    for (Eigen::SparseMatrix<adScalar>::InnerIterator it(L,colInd); it; ++it)
    {
      int row=it.row();
      diagVec(row)+=it.value();
    }
  }
  diagVec=diagVec.array()+1.0;


  solution=guess;
  adMat oldSolution(solution);



  adMat Lx(guess.rows(),1);
  setZero_ad(Lx,nDer);

  //Lx=L*solution;
  for (int i=0;i<solution.rows();i++){
    for (Eigen::SparseMatrix<adScalar>::InnerIterator it(L,i); it; ++it)
    {
      //Rcout<<it.value().derivatives()<<std::endl;
      int row=it.row();
      Lx(row)+=it.value()*solution(i);
    }
  }
  adMat errorVec(guess.rows(),1);
  setZero_ad(errorVec,nDer);
  // // 
  double diff=1;
  int count=0;
  while ((diff>tol)& (count<1000)){
  //   //solution=diagVec.array().pow(-1.0)*(rhs+Lx);
    for (int i=0;i<guess.rows();i++){
      solution(i)=(rhs(i)+Lx(i))/diagVec(i);
    }
    //Lx=L*solution;
    setZero_ad(Lx,nDer);
    for (int i=0;i<solution.rows();i++){
      for (Eigen::SparseMatrix<adScalar>::InnerIterator it(L,i); it; ++it)
      {
        //Rcout<<it.value().derivatives()<<std::endl;
        int row=it.row();
        Lx(row)+=it.value()*solution(i);
      }
    }
  // 
  //   //errorVec=diagVec.array()*solution.array()-Lx-rhs;
    for (int i=0;i<guess.rows();i++){
      errorVec(i)=diagVec(i)*solution(i)-Lx(i)-rhs(i);
    }
  // 
    diff=0;
    for (int i=0;i<errorVec.rows();i++){
      if (std::abs(errorVec(i,0).value())>diff){
        diff=std::abs(errorVec(i,0).value());
      }
    }
  //   oldSolution=solution;
    count+=1;
  }
  return solution;
}
// //
adMat matrixVectorMultiply(Eigen::MatrixXd X, adMat y){
  adMat val(X.rows(),1);
  setZero_ad(val,y(0).derivatives().size());
  //Rcout<<"start"<<val(0,0).derivatives()<<std::endl;
  for (int i=0;i<X.cols();i++){
    //Rcout<<"matrix vector deriv "<<y(i,0).derivatives()<<std::endl;
    //Rcout<<"matrix vector deriv "<<(y(i,0)*1.0).derivatives()<<std::endl;
    
    for (int j=0;j<X.rows();j++){
      val(j,i)+=X(j,i)*y(i,0);
    }
  }
  //Rcout<<"matrix vector product "<<val(0,0).derivatives()<<std::endl;
  return val;
}

void zero_I_L(spAdmat& I_L, Eigen::SparseMatrix<double> neighborMatrix, int nDer,bool diagonal){
  for (int colInd=0; colInd<neighborMatrix.outerSize(); ++colInd){
    for (Eigen::SparseMatrix<double>::InnerIterator it(neighborMatrix,colInd); it; ++it)
    {
      int row=it.row();
      int col=it.col();
      I_L.insert(row,col)=adZero(nDer);
    }
  }
  if (diagonal){
    for (int i=0;i<neighborMatrix.outerSize();i++){
      I_L.insert(i,i)=adZero(nDer);
    }
  }
}

void make_Laplacian(Eigen::SparseMatrix<adScalar>& I_L,
                   adMat mu,
                   int nDer,
                   bool diagonal){
  adMat diagVec(mu.rows(),1);
  setZero_ad(diagVec,nDer);
  for (int colInd=0; colInd<I_L.outerSize(); ++colInd){
    for (Eigen::SparseMatrix<adScalar>::InnerIterator it(I_L,colInd); it; ++it)
    {
      int row=it.row();
      int nNonzero=I_L.col(row).nonZeros();
      //the diagonal has been filled in for I_L............
      if (nNonzero>3){
        if (row>colInd){
          it.valueRef()=mu(colInd);
        }
        else{
          it.valueRef()=mu(row);
        }
      }
      else{
        it.valueRef()=adZero(nDer);
      }
    }
  }
  if (diagonal){
    for (int colInd=0; colInd<I_L.outerSize(); ++colInd){
      for (Eigen::SparseMatrix<adScalar>::InnerIterator it(I_L,colInd); it; ++it)
      {
        int row=it.row();
        if (row!=colInd){
          diagVec(row,0)-=it.value();
        }
      }
    }
    diagVec=diagVec.array()+1.0;
    for (int colInd=0;colInd<I_L.outerSize();++colInd){
      for (Eigen::SparseMatrix<adScalar>::InnerIterator it(I_L,colInd); it; ++it)
      {
        int row=it.row();
        if (row==colInd){
          it.valueRef()=diagVec(colInd,0);
        }
      }
    }
  }
}



adMat computeDiffusion(adMat params,
                       Eigen::MatrixXd X,
                       adMat init,
                       Eigen::SparseMatrix<double> neighborMatrix,
                       int nSpace,
                       int nTime,
                       double tol,
                       int nDer){

  int p=X.cols();

  adMat alpha(p,1),gamma(p,1);

  alpha=params.block(0,0,p,1);
  gamma=params.block(p,0,p,1);
  
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

  adMat U(nSpace,nTime);
  setZero_ad(U,nDer);
  U.col(0)=init;

  for (int k=0;k<boundaryPoints.size();k++){
    U.col(0)(boundaryPoints(k))=0*U.col(0)(boundaryPoints(k));
  }

  adMat muVec(X.rows(),1), lambdaVec(X.rows(),1);
  setZero_ad(muVec,nDer);
  setZero_ad(lambdaVec,nDer);

  muVec=matrixVectorMultiply(X,alpha).array().exp();
  lambdaVec = matrixVectorMultiply(X,gamma).array().exp();

  Eigen::Map<adMat> mu(muVec.data(),nSpace,nTime);
  Eigen::Map<adMat> lambda(lambdaVec.data(),nSpace,nTime);
  

  Eigen::SparseMatrix<adScalar> I_L(neighborMatrix.rows(),neighborMatrix.cols());
  
  bool diagonal=false;
  zero_I_L(I_L,neighborMatrix,nDer,diagonal);
  adMat rhs(nSpace,1);
  setZero_ad(rhs,nDer);

  for (int i=1;i<nTime;i++){
    //Kolmogorov-Fisher
    // U2-U1=HU2+Lambda*U1(1-U1)
    // (I-H)U2=U1+Lambda*U1(1-U1)
    rhs=U.col(i-1).array()+lambda.col(i-1).array()*U.col(i-1).array()*(1-U.col(i-1).array());
    make_Laplacian(I_L,mu.col(i),nDer,diagonal);

    //use BiCGSTAB
    if (diagonal){
      
      Eigen::BiCGSTAB<Eigen::SparseMatrix<adScalar> > solver;
      solver.setTolerance(tol);
      solver.compute(I_L);
      adMat guess(nSpace,1);
      setZero_ad(guess,nDer);
      guess=U.col(i-1);
      U.col(i)= solver.solveWithGuess(rhs,guess);
    }
    //use Jacobi iteration
    if (!diagonal){
      U.col(i)=solveJacobi(I_L,
          rhs,
          U.col(i-1),
          tol);
    }

    for (int k=0;k<boundaryPoints.size();k++){
      U.col(i)(boundaryPoints(k))=0*U.col(i)(boundaryPoints(k));
    }

    for (int k=0;k<nSpace;k++){
      if (U.col(i)(k).value()>1){
        U.col(i)(k).value()=1;
      }
    }
  }
  return U;
}

//
Eigen::MatrixXd getDiffusionFromAdmat(adMat U_ad){
  Eigen::MatrixXd U(U_ad.rows(),U_ad.cols());
  for (int i=0;i<U.cols();i++){
    for (int j=0;j<U.rows();j++){
      U(j,i)=U_ad(j,i).value();
    }
  }
  return U;
}

//[[Rcpp::export]]
Rcpp::List logLikelihood_ad(Eigen::VectorXd params,
                     Eigen::MatrixXd X,
                     Eigen::MatrixXd initCov,
                     Eigen::SparseMatrix<double> neighborMatrix,
                     int timeSteps,
                     double tol,
                     Eigen::VectorXi y,
                     Eigen::VectorXi yIndices){
  int p=X.cols();
  int nSpace=initCov.rows();
  adMat params_ad(params.size(),1),init_ad(nSpace,1);
  int nDer=params.size();
  for (int i=0;i<params.size();i++){
    params_ad(i).value()=params(i);
    params_ad(i).derivatives()=Eigen::VectorXd::Unit(nDer,i);
  }
  //
  // adScalar long_=params_ad(2*p);
  // adScalar lat_=params_ad(2*p+1);
  // adScalar sigma=exp(params_ad(2*p+2));
  // adScalar kappa=exp(params_ad(2*p+3));
  adMat distances(nSpace,1);
  setZero_ad(distances,nDer);
  setZero_ad(init_ad,nDer);
  // //
  distances=(initCov.col(0).array()-params_ad(2*p)).array().pow(2)+(initCov.col(1).array()-params_ad(2*p+1)).array().pow(2);
  init_ad=(distances/pow(exp(params_ad(2*p+2)),2)).array().exp();
  init_ad=1.0/(1.0+exp(params_ad(2*p+3))*init_ad.array());
  // adMat params,
  // Eigen::MatrixXd X,
  // Eigen::VectorXd init,
  // Eigen::SparseMatrix<double> neighborMatrix,
  // int nSpace,
  // int nTime,
  // bool KF,
  // double tol,
  // int nDer)
  adMat U(nSpace,timeSteps);
  U=computeDiffusion(params_ad.block(0,0,2*p,1),
                     X,
                     init_ad,
                     neighborMatrix,
                     nSpace,
                     timeSteps,
                     tol,
                     params_ad.rows());
  adScalar logLike;
  logLike=adZero(params_ad.rows());
  for (int i=0;i<yIndices.size();i++){
    if ((yIndices(i)-1)>(U.rows()*U.cols())-1){
      //Rcout<<yIndices(i)<<std::endl;
    }
    if ((yIndices(i)-1)<=(U.rows()*U.cols()-1)){
      //Rcout<<logLike.value()<<std::endl;
      if (U(yIndices(i)-1)<0.001){
        //Rcout<<U(yIndices(i)-1).value()<<std::endl;
      }
      //if (U(yIndices(i)-1)>0){
        logLike=logLike+y(i)*log(U(yIndices(i)-1))+(1-y(i))*log(1-U(yIndices(i)-1));
      //}
    }
  }
  return List::create(Named("logLike")=logLike.value(),
                      Named("derivative")=logLike.derivatives(),
                      Named("diffusion")=getDiffusionFromAdmat(U));
}
