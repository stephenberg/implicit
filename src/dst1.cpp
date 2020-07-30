<<<<<<< HEAD
#include <Eigen>
=======
#include <Rcpp.h>
#include <RcppEigen.h>
>>>>>>> 263e79fd2f86d4863e510a160792143bb747d30a
#include "fftPlan.h"

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

