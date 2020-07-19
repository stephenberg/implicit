#include <Rcpp.h>
#include <fftw3.h>
#include <RcppEigen.h>
#include "matrixFree.h"
using namespace Rcpp;
using namespace Eigen;

class fftPlan
{
private:

  fftw_plan p1,p2;
  double *input1,*output1,*output2;
public:
  int rows;
  int cols;
  fftPlan(int rows_,
          int cols_,
          double* input1_,
          double* output1_,
          double* output2_){
    rows=rows_;
    cols=cols_;
    input1=input1_;
    output1=output1_;
    output2=output2_;
    // fftw_plan fftw_plan_r2r_2d(int n0, int n1, double *in, double *out,
    //                            fftw_r2r_kind kind0, fftw_r2r_kind kind1,
    //                            unsigned flags);
    p1 = fftw_plan_r2r_2d(cols,
                         rows,
                         input1,
                         output1,
                         FFTW_RODFT00, 
                         FFTW_RODFT00,
                         FFTW_ESTIMATE);
    p2 = fftw_plan_r2r_2d(cols,
                          rows,
                          output1,
                          output2,
                          FFTW_RODFT00, 
                          FFTW_RODFT00,
                          FFTW_ESTIMATE);
  }
  

  
  //run discrete sine transform and rescale the output
  void dst_1(){
    fftw_execute(p1);
    Map<MatrixXd> y(output1,rows,cols);
    y=y/std::sqrt(2*2*(rows+1)*(cols+1));
  }
  
  void dst_1_reverse(){
    fftw_execute(p2);
    Map<MatrixXd> y(output2,rows,cols);
    y=y/std::sqrt(2*2*(rows+1)*(cols+1));
  }
  
  void multiplyBy_I_A_inverse(double lengthX=1,
                              double lengthY=1,
                              double delta=1){
    this->dst_1();
    
    Map<VectorXd> y(output1,rows*cols);
    int count=0;
    for (int j=0;j<cols;j++){
      for (int i=0;i<rows;i++){
        y(count)=y(count)*std::pow(1-delta*laplacianEigenvalue(i+1,
                   j+1,
                   rows,
                   cols,
                   lengthX,
                   lengthY),-1.0);
        count+=1;
      }
    }
    this->dst_1_reverse();
  }
  
  void multiplyBy_pow_neg_A(double power,
                            double lengthX=1,
                            double lengthY=1){
    this->dst_1();
    
    Map<VectorXd> y(output1,rows*cols);
    int count=0;
    for (int j=0;j<cols;j++){
      for (int i=0;i<rows;i++){
        y(count)=y(count)*std::pow(-laplacianEigenvalue(i+1,
                   j+1,
                   rows,
                   cols,
                   lengthX,
                   lengthY),power);
        count+=1;
      }
    }
    this->dst_1_reverse();
  }
  
  inline double laplacianEigenvalue(int i,
                                    int j,
                                    int rows,
                                    int cols,
                                    double lengthX,
                                    double lengthY){
    double lambda_i,lambda_j,lambda_ij;
    double h_y,h_x;
    h_y=lengthY/(rows+1);
    h_x=lengthX/(cols+1);
    lambda_i=-4/std::pow(h_y,2)*std::pow(std::sin(M_PI*i/(2*(rows+1))),2);
    lambda_j=-4/std::pow(h_x,2)*std::pow(std::sin(M_PI*j/(2*(cols+1))),2);
    lambda_ij=lambda_i+lambda_j;
    return lambda_ij;
  }
  
};





