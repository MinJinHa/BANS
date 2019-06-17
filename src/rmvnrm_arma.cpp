#include "RcppArmadillo.h"
#include "RcppEigen.h"
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]


arma::mat rmvnrm_arma(int n,
                      arma::vec mu,
                      arma::mat sigma) {
    int ncols = sigma.n_cols;
    arma::mat Y = arma::randn(n, ncols);
    return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}
