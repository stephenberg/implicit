//#include <Rcpp.h>
#include <Eigen>
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

int main() {

}

