<<<<<<< HEAD
#include <Eigen>
=======
#include <RcppEigen.h>
>>>>>>> 263e79fd2f86d4863e510a160792143bb747d30a
#include <math.h>
#ifndef __CONVENIENCE__
#define __CONVENIENCE__

inline double logit(double x){
  return std::log(x/(1-x));
}

inline double expit(double x){
  return (1/(1+std::exp(-x)));
}

<<<<<<< HEAD
#endif
=======
#endif
>>>>>>> 263e79fd2f86d4863e510a160792143bb747d30a
