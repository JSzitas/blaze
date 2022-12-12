#ifndef XREG_SOLVER
#define XREG_SOLVER

#include "third_party/eigen.h"
#include "utils/utils.h"

template <typename U=double> struct lm_coef {
  lm_coef<U>(){
    this->coef = std::vector<U>(0);
    // this->covar= std::vector<U>(0);
    this->intercept = false;
  };
  lm_coef<U>(int xreg_ncol, bool intercept){
    this->coef = std::vector<U>(xreg_ncol + intercept);
    // this->covar= std::vector<U>(xreg_ncol + intercept);
    this->intercept = intercept;
  };
  // move coefficients when creating, copy intercept
  lm_coef<U>( std::vector<U> coef,
              bool intercept ) :
    coef(std::move(coef)),
    intercept(intercept) {};
  const int size() const {
    return this->coef.size();
  }
  const std::vector<U> data() const {
    return this->coef();
  }
  const bool has_intercept() const {
    return this->intercept;
  };
  U operator[](int i) {
    return coef[i];
  }
  const U operator[](int i) const {
    return coef[i];
  }
private:
  std::vector<U> coef;
  bool intercept = false;
};

template <typename U=double> std::vector<U> est_variance(
  Eigen::Matrix<U, Eigen::Dynamic, 1> y,
  Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> X,
  Eigen::Matrix<U, Eigen::Dynamic, 1> coef,
  Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> covar_mat ) {

  U sigma_sqr = 0;
  auto y_hat = X * coef;
  for( int i=0 ; i < y.size(); i++ ) {
    sigma_sqr += pow(y[i] - y_hat[i], 2 );
  }
  sigma_sqr /= (U)(X.rows() - X.cols());

  // take diagonal elements of covariance matrix
  auto covars = covar_mat.diagonal().array();
  for( auto &covar:covars ) {
    covar *= sigma_sqr;
    covar = sqrt(covar);
  }

  std::vector<U> result(X.cols());
  for( int i=0; i < result.size(); i++) {
    result[i] = covars[i];
  }
  return result;
}


// map 2 vectors to Eigen matrices and call solve
template <typename U=double> lm_coef<U> xreg_coef(
  std::vector<U> &y,
  std::vector<U> & xreg,
  bool use_intercept = true) {

  const int n = y.size();
  Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> new_mat = Eigen::Map<
    Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
    >(xreg.data(), n, xreg.size()/n);
  if( use_intercept ) {
    new_mat.conservativeResize(Eigen::NoChange, new_mat.cols()+1);
    new_mat.col(new_mat.cols()-1) = Eigen::Matrix<U, Eigen::Dynamic, 1>::Constant(n, 1, 1);
  }
  Eigen::Matrix<U, Eigen::Dynamic, 1> new_vec = Eigen::Map<Eigen::Matrix<U, Eigen::Dynamic, 1>>(y.data(), n , 1);

  auto decomp = new_mat.completeOrthogonalDecomposition();
  Eigen::Matrix<U, Eigen::Dynamic, 1> res = decomp.solve(new_vec);
  std::vector<double> result(res.data(), res.data() + res.rows() * res.cols());
  // put intercept at the start if it is used
  if( use_intercept ) {
    auto intercept = result[result.size()-1];
    result.resize(result.size()-1);
    std::reverse(result.begin(), result.end());
    result.push_back(intercept);
    std::reverse(result.begin(), result.end());
  }
  // compute standard errors - this is mildly annoying - perhaps it is faster to
  // only invert once and then do the multiplies maybe?
  // Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> covar_mat = (new_mat.transpose() * new_mat).inverse();
  // auto covariances = est_variance(new_vec, new_mat, res, covar_mat);

  lm_coef<U> final( result, use_intercept );
  return final;
}

template <typename U=double> std::vector<U> predict( lm_coef<U> coef,
                                                     std::vector<U> & new_xreg ) {

  int n_cols = coef.size();
  int n = new_xreg.size()/(n_cols - coef.intercept);
  Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> new_mat;

  if( coef.intercept ){
    new_mat = Eigen::Map<Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>>(new_xreg.data(), n, n_cols-1);
    new_mat.conservativeResize(Eigen::NoChange, new_mat.cols()+1);
    new_mat.col(new_mat.cols()-1) = Eigen::Matrix<U, Eigen::Dynamic, 1>::Constant(n, 1, 1);
  }
  else {
    new_mat = Eigen::Map<Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>>(new_xreg.data(), n, n_cols);
  }
  Eigen::Matrix<U, Eigen::Dynamic, 1> coef_ = Eigen::Map<Eigen::Matrix<U, Eigen::Dynamic, 1>>(coef, n_cols , 1);
  // predictions
  Eigen::Matrix<U, Eigen::Dynamic, 1> result = new_mat * coef;

  std::vector<double> res(result.data(), result.data() + result.rows() );
  return res;
}

template <typename U=double> struct fixed_xreg{
  fixed_xreg<double>( std::vector<U> & xreg, const int n, const bool use_intercept = true) {
    Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> new_mat = Eigen::Map<
      Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
      >(xreg.data(), n, xreg.size()/n);
      if( use_intercept ) {
        new_mat.conservativeResize(Eigen::NoChange, new_mat.cols()+1);
        new_mat.col(new_mat.cols()-1) = Eigen::Matrix<U, Eigen::Dynamic, 1>::Constant(n, 1, 1);
      }
      this->fixed = new_mat;
  }
  fixed_xreg<U>(Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> &X) {
    this->fixed = X;
  }
  Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> fixed;
};

template <typename U=double> std::vector<U> operator*( fixed_xreg<U> xreg,
                                                       lm_coef<U> xreg_pars) {
  Eigen::Matrix<U, Eigen::Dynamic, 1> coef_vec = Eigen::Map<Eigen::Matrix<U, Eigen::Dynamic, 1>>(xreg_pars.data(), xreg_pars.size(), 1);
  auto result = xreg.fixed * coef_vec;
  return std::vector<double>(result.data(), result.data() + result.rows() );
}

#endif
