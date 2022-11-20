#ifndef XREG_SOLVER
#define XREG_SOLVER

#include <RcppEigen.h>
#include "Eigen/Dense"
#include "utils.h"

// [[Rcpp::depends(RcppEigen)]]

template <typename U=double> struct lm_coef {
  // move coefficients when creating, copy intercept
  lm_coef<U>( std::vector<U> coef, bool intercept ) :
    coef(std::move(coef)),
    intercept(intercept) {};
  std::vector<U> coef;
  bool intercept = false;
};

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

  Eigen::Matrix<U, Eigen::Dynamic, 1> res = new_mat.completeOrthogonalDecomposition().solve(new_vec);
  std::vector<double> result(res.data(), res.data() + res.rows() * res.cols());
  // put intercept at the start if it is used
  if( use_intercept ) {
    // auto a = result[0];
    auto intercept = result[result.size()-1];
    result.resize(result.size()-1);
    std::reverse(result.begin(), result.end());
    result.push_back(intercept);
    std::reverse(result.begin(), result.end());
  }
  lm_coef<U> final( result, use_intercept );
  return final;
}

template <typename U=double>std::vector<U> predict( lm_coef<U> coef,
                                                    std::vector<U> & new_xreg ) {

  int n_cols = coef.coef.size();
  int n = new_xreg/n_cols - coef.intercept;
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





#endif
