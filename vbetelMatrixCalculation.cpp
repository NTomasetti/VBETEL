// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::List matrixCalculations(arma::mat g, arma::mat hessian, arma::vec lambdaHat, arma::vec exponent, arma::cube dgdt){
  using std::pow; using std::log;
  int n = g.n_rows;
  int d = g.n_cols;

  
  arma::mat dh2dtlam(d, d, arma::fill::zeros);
  for(int i = 0; i < n; ++i){
    dh2dtlam += 1.0 / n * (exponent(i) * (dgdt.slice(i) + g.row(i).t() * lambdaHat.t() * dgdt.slice(i)));
  }
    
  arma::mat dlamdt = - hessian.i() * dh2dtlam;
  arma::mat productRule(d, n, arma::fill::zeros);
  for(int i = 0; i < n; ++i){
    productRule.col(i) = dlamdt.t() * g.row(i).t() + dgdt.slice(i).t() * lambdaHat;
  }
  
  arma::vec dpdt(d, arma::fill::zeros);
  for(int j = 0; j < d; ++j){
    dpdt(j) = sum(productRule.row(j)) - n / sum(exponent) * as_scalar(productRule.row(j) * exponent);
  }
  double logp = sum(log(exponent)) - n * log(sum(exponent));
  return Rcpp::List::create(Rcpp::Named("grad") = dpdt,
                            Rcpp::Named("val") = logp);
  
}
