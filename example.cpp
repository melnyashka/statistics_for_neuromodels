#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
NumericVector mult_rcpp(NumericVector x, NumericVector y) {
  return x * y;
}

// [[Rcpp::export]]
NumericMatrix add_rcpp(NumericMatrix A, NumericMatrix B) {
  const Map<MatrixXd> Ae(as<Map<MatrixXd>>(A));
  const Map<MatrixXd> Be(as<Map<MatrixXd>>(B));
  MatrixXd sum = Ae + Be;
  return(wrap(sum));
}

// [[Rcpp::export]]
double edet(NumericMatrix A){
  // determinant of the square matrix, as required in statistics
  const Map<MatrixXd> Ae(as<Map<MatrixXd>>(A));
  MatrixXd prod = Ae*Ae;
  return(prod.determinant());
}

// [[Rcpp::export]]
double sum_rcpp(NumericVector xx) {
  const Map<VectorXd> x(as<Map<VectorXd>>(xx));
  return x.sum();
}
                      
// [[Rcpp::export]]
NumericVector add_rcpp_vec(NumericVector x, NumericVector y) {
  return x + y;
}

/*** R

cat("It's finished!")

*/
