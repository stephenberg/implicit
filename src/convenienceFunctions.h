#include <Eigen>
#include <math.h>
#ifndef __CONVENIENCE__
#define __CONVENIENCE__

inline double logit(double x){
  return std::log(x/(1-x));
}

inline double expit(double x){
  return (1/(1+std::exp(-x)));
}

#endif