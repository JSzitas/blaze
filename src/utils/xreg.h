#ifndef XREG_SOLVER
#define XREG_SOLVER

#include "third_party/eigen.h"
#include "utils/utils.h"

template <typename scalar_t> using EigVec = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;
template <typename scalar_t> using EigMat = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;

template <typename scalar_t = double> void add_constant(EigMat<scalar_t> & mat) {
  mat.conservativeResize(Eigen::NoChange, mat.cols() + 1);
  mat.col(mat.cols() - 1) = EigVec<scalar_t>::Constant(mat.rows(), 1, mat.rows());
}

template <typename scalar_t = double> void add_drift(EigMat<scalar_t> & mat) {
  mat.conservativeResize(Eigen::NoChange, mat.cols() + 1);
  mat.col(mat.cols() - 1) = EigVec<scalar_t>::LinSpaced(mat.rows(), 1, mat.rows());
}

template <typename scalar_t = double> void put_last_item_first( std::vector<scalar_t> & x) {
  const auto last_elem = x.size() - 1;
  const auto item = x[last_elem];
  x.resize(last_elem);
  std::reverse(x.begin(), x.end());
  x.push_back(item);
  std::reverse(x.begin(), x.end());
}

template <typename scalar_t = double> std::vector<scalar_t> solve_ortho_decomp(
  EigMat<scalar_t> &mat,
  EigVec<scalar_t> &vec) {
  auto decomp = mat.completeOrthogonalDecomposition();
  EigVec<scalar_t> res = decomp.solve(vec);
  std::vector<double> result(res.data(), res.data() + res.rows() * res.cols());
  return result;
}


// map 2 vectors to Eigen matrices and call solve
template <typename scalar_t = double>
std::vector<scalar_t> xreg_coef(
    std::vector<scalar_t> &y,
    std::vector<scalar_t> &xreg,
    const bool use_intercept = true,
    const bool use_drift = false) {
  const auto n = y.size();
  EigMat<scalar_t> mat = Eigen::Map<EigMat<scalar_t>>( xreg.data(), n, xreg.size() / n);
  if (use_drift) add_drift(mat);
  if (use_intercept) add_constant(mat);
  EigVec<scalar_t> vec = Eigen::Map<EigVec<scalar_t>>(y.data(), n, 1);
  // solve using orthogonal decomposition
  auto result = solve_ortho_decomp(mat, vec);
  // put intercept at the start if it is used
  if (use_intercept) put_last_item_first(result);
  return result;
}

// map 2 vectors to Eigen matrices and call solve
template <typename scalar_t = double>
std::vector<scalar_t> xreg_coef(
    std::vector<scalar_t> &y,
    std::vector<std::vector<scalar_t>> &xreg,
    const bool use_intercept = true,
    const bool use_drift = false) {

  const auto n = y.size();
  const auto ncol = xreg.size();
  auto _xreg = flatten_vec(xreg);
  EigMat<scalar_t> mat = Eigen::Map<EigMat<scalar_t>>(_xreg.data(), n, ncol);
  if (use_drift) add_drift(mat);
  if (use_intercept) add_constant(mat);
  EigVec<scalar_t> vec = Eigen::Map<EigVec<scalar_t>>(y.data(), n, 1);
  auto result = solve_ortho_decomp(mat, vec);
  // put intercept at the start if it is used
  if (use_intercept) put_last_item_first(result);
  return result;
}

template<typename scalar_t=double> std::vector<scalar_t> predict(
  const size_t n,
  std::vector<scalar_t> &coef,
  const bool use_intercept,
  const bool use_drift,
  std::vector<std::vector<scalar_t>> &xreg ) {

  std::vector<scalar_t> _xreg = flatten_vec(xreg);
  EigMat<scalar_t> mat = Eigen::Map<EigMat<scalar_t>>(_xreg.data(), n, xreg.size());
  if (use_drift) add_drift(mat);
  if (use_intercept) add_constant(mat);
  EigVec<scalar_t> _coef = Eigen::Map<EigVec<scalar_t>>( coef.data(), coef.size(), 1);
  EigVec<scalar_t> res = mat * _coef;
  std::vector<scalar_t> result(res.size());
  for(size_t i=0; i< result.size(); i++) result[i] = res[i];
  return result;
}

template < typename scalar_t = double> EigMat<scalar_t> vec_to_mat(
  const std::vector<std::vector<scalar_t>> &xreg,
  const size_t n,
  const bool use_intercept,
  const bool use_drift
  ) {
  std::vector<double> _xreg = flatten_vec(xreg);
  EigMat<scalar_t> mat = Eigen::Map<EigMat<scalar_t>>( _xreg.data(), n, xreg.size());
  if (use_drift) add_drift(mat);
  if (use_intercept) add_constant(mat);
  return mat;
}

#endif
