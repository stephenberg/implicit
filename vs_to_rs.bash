grep -RiIl '<Eigen>' src/* | xargs sed -i 's/<Eigen>/<RcppEigen.h>/g'
grep -RiIl '^\/\/#include <Rcpp.h>' src/* | xargs sed -i 's/^\/\/#include <Rcpp.h>/#include <Rcpp.h>/g'
grep -RiIl '^\/\/using namespace Rcpp;' src/* | xargs sed -i 's/^\/\/using namespace Rcpp;/using namespace Rcpp;/g'
grep -RiIl 'std::cout' src/* | xargs sed -i 's/std::cout/Rcout/g'










