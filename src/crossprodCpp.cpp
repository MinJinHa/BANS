#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
Eigen::MatrixXd crossprodCpp(const Eigen::MatrixXd B, const Eigen::MatrixXd C) {
  unsigned int n = B.cols(), k = C.cols();
  Eigen::MatrixXd Y = MatrixXd::Zero(n,k);
  Y = B.adjoint() * C;
  return Y;
}
  
