grep -RiIl '<RcppEigen.h>' src/* | xargs sed -i 's/<RcppEigen.h>/<Eigen>/g'
grep -RiIl '^#include <Rcpp.h>' src/* | xargs sed -i 's/^#include <Rcpp.h>/\/\/#include <Rcpp.h>/g'
grep -RiIl '^using namespace Rcpp;' src/* | xargs sed -i 's/^using namespace Rcpp;/\/\/using namespace Rcpp;/g'
grep -RiIl 'Rcout' src/* | xargs sed -i 's/Rcout/std::cout/g'










