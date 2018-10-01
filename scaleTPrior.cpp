// [[Rcpp::depends(rstan)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <stan/math.hpp>
#include <Eigen/Dense>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace arma;
using namespace std;
using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::MatrixXd; 
using Eigen::Map;

struct scaleT {
  const vec scale;
  const vec df;
  scaleT(const vec scIn, const vec dfIn) :
   scale(scIn), df(dfIn) {}
  template <typename T> //
  T operator ()(const Matrix<T, Dynamic, 1>& theta)
    const{
    using std::log; using std::exp; using std::pow; using std::sqrt; using std::fabs; using std::lgamma;
    int dim = theta.rows();
    
    T logDens = 0;
    
    for(int i = 0; i < dim; ++i){
      logDens += - log(scale(i)) + lgamma((df(i) + 1)/2) - 0.5 * log(df(i) * 3.14159) - lgamma(df(i)/2) - 
        ((df(i) + 1) / 2) * log(1 + pow(theta(i) / scale(i), 2) / df(i));
    }
    return logDens;
  }
};

// [[Rcpp::export]]
Rcpp::List gradScaleT(vec scale, vec df, Rcpp::NumericMatrix thetaIn){
  Map<MatrixXd> theta(Rcpp::as<Map<MatrixXd> >(thetaIn));
  double eval;
  int dim = theta.rows();
  Matrix<double, Dynamic, 1>  grad(dim);
  
  scaleT logp(scale, df);
  stan::math::set_zero_all_adjoints();
  stan::math::gradient(logp, theta, eval, grad);
  return Rcpp::List::create(Rcpp::Named("grad") = grad,
                            Rcpp::Named("val") = eval);
}
