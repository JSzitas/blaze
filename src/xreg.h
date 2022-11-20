#ifndef XREG_SOLVER
#define XREG_SOLVER

#include <RcppEigen.h>
#include "Eigen/Dense"
#include "utils.h"

// [[Rcpp::depends(RcppEigen)]]

// map 2 vectors to Eigen matrices and call solve
template <typename U=double> std::vector<U> xreg_coef(
  std::vector<U> &y,
  std::vector<U> & xreg,
  bool use_intercept = true) {

  const int n = y.size();
  Eigen::MatrixXd new_mat = Eigen::Map<Eigen::MatrixXd>(xreg.data(), n, xreg.size()/n);
  if( use_intercept ) {
    new_mat.conservativeResize(Eigen::NoChange, new_mat.cols()+1);
    new_mat.col(new_mat.cols()-1) = Eigen::Matrix<U, Eigen::Dynamic, 1>::Constant(n, 1, 1);
  }
  Eigen::VectorXd new_vec = Eigen::Map<Eigen::VectorXd>(y.data(), n , 1);

  Eigen::VectorXd res = new_mat.completeOrthogonalDecomposition().solve(new_vec);
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
  return result;
}

#endif
