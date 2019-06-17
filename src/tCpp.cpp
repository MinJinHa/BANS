#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
Eigen::MatrixXd tCpp(const Eigen::MatrixXd A) {
  unsigned int n = A.rows(), k = A.cols();
  Eigen::MatrixXd Y = MatrixXd::Zero(n,k);
  Y = A.transpose();
  return Y;
}
