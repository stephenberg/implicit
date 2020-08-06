#include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>

//this file contains:

//fick_L_f(): multiply by fickian diffusion laplacian for rectangular grid

//homogeneous_L_f(): multiply by spatially homogeneous diffusion laplacian 
//for rectangular grid

//general_Homogeneous_L_f(): multiply by spatially homogeneous diffusion laplacian 
//for non-rectangular grid; allows no-flux boundary condition

//to add:
//general_fick_L_f(): multiply by fickian diffusion laplacian for non-rectangular grid,
//allow no-flux boundary condition



//implicit matrix vector multiplication with inhomogeneous diffusion coefficients
//diffusion coefficient in middle of Laplacian (Fickian diffusion)
Eigen::VectorXd fick_L_f(const Eigen::Ref<const Eigen::VectorXd>& mu_,
                         Eigen::VectorXd& f_,
                         int rows,
                         int cols,
                         bool dirichlet=true,
                         double lengthX=1,
                         double lengthY=1){
  VectorXd b_;
  b_.setZero(rows*cols);
  
  double h_y,h_x;
  h_y=lengthY/(rows+1);
  h_x=lengthX/(cols+1);
  
  //use value of diffusion coefficient for boundary cells
  const Map<const MatrixXd> mu(mu_.data(),rows+2,cols+2);
  
  //f and b implicitly defined as 0 on boundary points
  Map<MatrixXd> f(f_.data(),rows,cols);
  Map<MatrixXd> b(b_.data(),rows,cols);
  
  //Fickian diffusion
  
  //d2dx2
  for (int j=0;j<cols;j++){
    for (int i=0;i<rows;i++){
      double leftMu,rightMu;
      
      //note: mu contains the boundary points, in particular
      //mu(i+1,j+1) corresponds to the point f(i,j)
      leftMu=(mu(i+1,j)+mu(i+1,j+1))/2;
      rightMu=(mu(i+1,j+1)+mu(i+1,j+2))/2;
      if (j==0){
        if (dirichlet){
          b(i,j)+=(rightMu*f(i,j+1)-(leftMu+rightMu)*f(i,j))/std::pow(h_x,2);
        }
        else{
          b(i,j)+=(rightMu*f(i,j+1)-(rightMu)*f(i,j))/std::pow(h_x,2);          
        }
      }
      if (j==(cols-1)){
        if (dirichlet){
          b(i,j)+=(leftMu*f(i,j-1)-(leftMu+rightMu)*f(i,j))/std::pow(h_x,2);
        }
        else{
          b(i,j)+=(leftMu*f(i,j-1)-(leftMu)*f(i,j))/std::pow(h_x,2);
        }
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
        if (dirichlet){
          b(i,j)+=(lowerMu*f(i+1,j)-(lowerMu+upperMu)*f(i,j))/std::pow(h_y,2);
        }
        else{
          b(i,j)+=(lowerMu*f(i+1,j)-(lowerMu)*f(i,j))/std::pow(h_y,2);
        }
      }
      if (i==(rows-1)){
        if (dirichlet){
          b(i,j)+=(upperMu*f(i-1,j)-(lowerMu+upperMu)*f(i,j))/std::pow(h_y,2);
        }
        else{
          b(i,j)+=(upperMu*f(i-1,j)-(upperMu)*f(i,j))/std::pow(h_y,2);
        }
      }
      if ((i>0)& (i<(rows-1))){
        b(i,j)+=(lowerMu*f(i+1,j)+upperMu*f(i-1,j)-(lowerMu+upperMu)*f(i,j))/std::pow(h_y,2);
      }
    }
  }
  return b_;
}


//implicit multiplication by the Laplacian with homogeneous diffusion
Eigen::VectorXd homogeneous_L_f(const Ref<const Eigen::VectorXd>& f_,
                                int rows,
                                int cols,
                                bool dirichlet=true,
                                double lengthX=1,
                                double lengthY=1){
  
  double boundaryMultiplier;
  if (dirichlet){
    boundaryMultiplier=2.0;
  }
  else{
    boundaryMultiplier=1.0;
  }
  
  VectorXd b_;
  b_.setZero(rows*cols);
  
  double h_y,h_x;
  h_y=lengthY/(rows+1);
  h_x=lengthX/(cols+1);
  
  
  //f and b implicitly defined as 0 on boundary points
  const Map<const MatrixXd> f(f_.data(),rows,cols);
  Map<MatrixXd> b(b_.data(),rows,cols);
  
  //d2dx2
  for (int j=0;j<cols;j++){
    for (int i=0;i<rows;i++){
      if (j==0){
        b(i,j)+=(f(i,j+1)-boundaryMultiplier*f(i,j))/std::pow(h_x,2);
      }
      if (j==(cols-1)){
        b(i,j)+=(f(i,j-1)-boundaryMultiplier*f(i,j))/std::pow(h_x,2);
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
        b(i,j)+=(f(i+1,j)-boundaryMultiplier*f(i,j))/std::pow(h_y,2);
      }
      if (i==(rows-1)){
        b(i,j)+=(f(i-1,j)-boundaryMultiplier*f(i,j))/std::pow(h_y,2);
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
Eigen::VectorXd general_Homogeneous_Lf(Eigen::VectorXd& f_,
                                       int rows,
                                       int cols,
                                       Eigen::VectorXi& internalPoints_,
                                       bool dirichlet=true,
                                       double lengthX=1,
                                       double lengthY=1){
  VectorXd b_;
  b_.setZero(rows*cols);
  
  Map<MatrixXi> internal(internalPoints_.data(),rows,cols);
  
  double h_y,h_x;
  h_y=lengthY/(rows+1);
  h_x=lengthX/(cols+1);
  
  
  //f and b implicitly defined as 0 on boundary points (Dirichlet) or with normal derivative 0 (Neumann)
  Map<MatrixXd> f(f_.data(),rows,cols);
  Map<MatrixXd> b(b_.data(),rows,cols);
  
  //d2dx2
  for (int j=0;j<cols;j++){
    for (int i=0;i<rows;i++){
      if (j==0){
        //Dirichlet
        if (dirichlet){
          b(i,j)+=internal(i,j)*(internal(i,j+1)*f(i,j+1)-2*f(i,j))/std::pow(h_x,2);
        }
        //Neumann
        else{
          b(i,j)+=internal(i,j)*(internal(i,j+1)*(f(i,j+1)-f(i,j)))/std::pow(h_x,2);
        }
      }
      if (j==(cols-1)){
        //Dirichlet
        if (dirichlet){
          b(i,j)+=internal(i,j)*(internal(i,j-1)*f(i,j-1)-2*f(i,j))/std::pow(h_x,2);
        }
        //Neumann
        else{
          b(i,j)+=internal(i,j)*(internal(i,j-1)*(f(i,j-1)-f(i,j)))/std::pow(h_x,2);
        }
      }
      if ((j>0)&(j<(cols-1))){
        //Dirichlet
        if (dirichlet){
          b(i,j)+=internal(i,j)*(internal(i,j+1)*f(i,j+1)+internal(i,j-1)*f(i,j-1)-2*f(i,j))/std::pow(h_x,2);
        }
        //Neumann
        else{
          b(i,j)+=internal(i,j)*(internal(i,j+1)*(f(i,j+1)-f(i,j))-internal(i,j-1)*(f(i,j)-f(i,j-1)))/std::pow(h_x,2);
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
          b(i,j)+= internal(i, j)*(internal(i+1,j)*f(i+1,j)-2*f(i,j))/std::pow(h_y,2);
        }
        //Neumann
        else{
          b(i,j)+=internal(i,j)*(internal(i+1,j)*(f(i+1,j)-f(i,j)))/std::pow(h_y,2);
        }
      }
      if (i==(rows-1)){
        //Dirichlet
        if (dirichlet){
          b(i,j)+= internal(i, j)*(internal(i-1,j)*f(i-1,j)-2*f(i,j))/std::pow(h_y,2);
        }
        else{
          b(i,j)+=internal(i,j)*(internal(i-1,j)*(f(i-1,j)-f(i,j)))/std::pow(h_y,2);
        }
      }
      if ((i>0)&(i<(rows-1))){
        //Dirichlet
        if (dirichlet){
          b(i,j)+= internal(i, j)*(internal(i+1,j)*f(i+1,j)+internal(i-1,j)*f(i-1,j)-2*f(i,j))/std::pow(h_y,2);
        }
        //Neumann
        else{
          b(i,j)+=internal(i,j)*(internal(i-1,j)*(f(i-1,j)-f(i,j))-internal(i+1,j)*(f(i,j)-f(i+1,j)))/std::pow(h_y,2);
        }
      }
    }
  }
  return b_;
}


//implicit multiplication by the Laplacian with homogeneous diffusion
//[[Rcpp::export]]
Eigen::VectorXd general_Fick_Lf(Eigen::VectorXd& mu_,
                                Eigen::VectorXd& f_,
                                int rows,
                                int cols,
                                Eigen::VectorXi& internalPoints_,
                                bool dirichlet=true,
                                double lengthX=1,
                                double lengthY=1){
  VectorXd b_;
  b_.setZero(rows*cols);
  
  
  Map<MatrixXi> internal(internalPoints_.data(),rows,cols);
  
  double h_y,h_x;
  h_y=lengthY/(rows+1);
  h_x=lengthX/(cols+1);
  
  //use value of diffusion coefficient for boundary cells
  Map<MatrixXd> mu(mu_.data(),rows+2,cols+2);
  
  //f and b implicitly defined as 0 on boundary points
  Map<MatrixXd> f(f_.data(),rows,cols);
  Map<MatrixXd> b(b_.data(),rows,cols);
  
  //Fickian diffusion
  
  //d2dx2
  for (int j=0;j<cols;j++){
    for (int i=0;i<rows;i++){
      double leftMu,rightMu;
      
      //note: mu contains the boundary points, in particular
      //mu(i+1,j+1) corresponds to the point f(i,j)
      leftMu=(mu(i+1,j)+mu(i+1,j+1))/2;
      rightMu=(mu(i+1,j+1)+mu(i+1,j+2))/2;
      
      if (j==0){
        //Dirichlet
        if (dirichlet){
          b(i,j)+=internal(i,j)*(internal(i,j+1)*rightMu*f(i,j+1)-(leftMu+rightMu)*f(i,j))/std::pow(h_x,2);
        }
        //Neumann
        else{
          b(i,j)+=internal(i,j)*(internal(i,j+1)*rightMu*(f(i,j+1)-f(i,j)))/std::pow(h_x,2);
        }
      }
      if (j==(cols-1)){
        //Dirichlet
        if (dirichlet){
          b(i,j)+=internal(i,j)*(internal(i,j-1)*leftMu*f(i,j-1)-(leftMu+rightMu)*f(i,j))/std::pow(h_x,2);
        }
        //Neumann
        else{
          b(i,j)+=internal(i,j)*(internal(i,j-1)*leftMu*(f(i,j-1)-f(i,j)))/std::pow(h_x,2);
        }
      }
      if ((j>0)&(j<(cols-1))){
        //Dirichlet
        if (dirichlet){
          b(i,j)+=internal(i,j)*(internal(i,j+1)*rightMu*f(i,j+1)+internal(i,j-1)*leftMu*f(i,j-1)-(leftMu+rightMu)*f(i,j))/std::pow(h_x,2);
        }
        //Neumann
        else{
          b(i,j)+=internal(i,j)*(internal(i,j+1)*rightMu*(f(i,j+1)-f(i,j))-internal(i,j-1)*leftMu*(f(i,j)-f(i,j-1)))/std::pow(h_x,2);
        }
      }
    }
  }
  
  //d2dy2
  for (int j=0;j<cols;j++){
    for (int i=0;i<rows;i++){
      
      double upperMu,lowerMu;
      upperMu=(mu(i,j+1)+mu(i+1,j+1))/2;
      lowerMu=(mu(i+1,j+1)+mu(i+2,j+1))/2;
      
      if (i==0){
        //Dirichlet
        if (dirichlet){
          b(i,j)+= internal(i, j)*(internal(i+1,j)*lowerMu*f(i+1,j)-(lowerMu+upperMu)*f(i,j))/std::pow(h_y,2);
        }
        //Neumann
        else{
          b(i,j)+=internal(i,j)*(internal(i+1,j)*lowerMu*(f(i+1,j)-f(i,j)))/std::pow(h_y,2);
        }
      }
      if (i==(rows-1)){
        //Dirichlet
        if (dirichlet){
          b(i,j)+= internal(i, j)*(internal(i-1,j)*upperMu*f(i-1,j)-(lowerMu+upperMu)*f(i,j))/std::pow(h_y,2);
        }
        else{
          b(i,j)+=internal(i,j)*(internal(i-1,j)*upperMu*(f(i-1,j)-f(i,j)))/std::pow(h_y,2);
        }
      }
      if ((i>0)&(i<(rows-1))){
        //Dirichlet
        if (dirichlet){
          b(i,j)+= internal(i, j)*(internal(i+1,j)*lowerMu*f(i+1,j)+internal(i-1,j)*upperMu*f(i-1,j)-(lowerMu+upperMu)*f(i,j))/std::pow(h_y,2);
        }
        //Neumann
        else{
          b(i,j)+=internal(i,j)*(internal(i-1,j)*upperMu*(f(i-1,j)-f(i,j))-internal(i+1,j)*lowerMu*(f(i,j)-f(i+1,j)))/std::pow(h_y,2);
        }
      }
    }
  }
  return b_;
}
