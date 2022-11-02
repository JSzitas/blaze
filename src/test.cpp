#include <Rcpp.h>
using namespace Rcpp;
using namespace RcppEigen;

#include <RcppEigen.h>
// #include <Eigen/Dense>
// #include <Eigen/Map>

// [[Rcpp::depends(RcppEigen)]]

// map 2 vectors to Eigen matrices and call solve
// [[Rcpp::export]]
std::vector<double> solve_mat_vec( std::vector<double> &mat,
                                   std::vector<double> &vec ) {
  const int n = vec.size();
  // double *mat_ptr = mat.data();
  // double *vec_ptr = vec.data();
  Eigen::MatrixXd new_mat = Eigen::Map<Eigen::MatrixXd>(mat.data(), n, n);
  Eigen::VectorXd new_vec = Eigen::Map<Eigen::VectorXd>(vec.data(), n , 1);

  Eigen::VectorXd res = new_mat.completeOrthogonalDecomposition().solve(new_vec);

  std::vector<double> result(res.data(), res.data() + res.rows() * res.cols());
  return result;
}
