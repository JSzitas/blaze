#ifndef XREG_SOLVER
#define XREG_SOLVER

#include "third_party/eigen.h"
#include "utils/utils.h"

template <typename scalar_t> using EigVec = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;
template <typename scalar_t> using EigMat = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;

template <typename U = double> struct lm_coef {
  lm_coef<U>() {
    this->coef = std::vector<U>(0);
    this->intercept = false;
  };
  lm_coef<U>(int xreg_ncol, bool intercept) {
    this->coef = std::vector<U>(xreg_ncol + intercept);
    this->intercept = intercept;
  };
  // move coefficients when creating, copy intercept
  lm_coef<U>(std::vector<U> coef, bool intercept)
      : coef(std::move(coef)), intercept(intercept){};
  const int size() const { return this->coef.size(); }
  std::vector<U> data() const { return this->coef; }
  const bool has_intercept() const { return this->intercept; };
  U &operator[](int i) { return coef[i]; }
  const U get_intercept() const { return coef.back(); }
  // private:
  std::vector<U> coef;
  bool intercept = false;
};

// map 2 vectors to Eigen matrices and call solve
template <typename U = double>
lm_coef<U> xreg_coef(std::vector<U> &y, std::vector<U> &xreg,
                     bool use_intercept = true) {

  const int n = y.size();
  EigMat<U> new_mat = Eigen::Map<EigMat<U>>( xreg.data(), n, xreg.size() / n);
  if (use_intercept) {
    new_mat.conservativeResize(Eigen::NoChange, new_mat.cols() + 1);
    new_mat.col(new_mat.cols() - 1) = EigVec<U>::Constant(n, 1, 1);
  }
  EigVec<U> new_vec = Eigen::Map<EigVec<U>>(y.data(), n, 1);

  auto decomp = new_mat.completeOrthogonalDecomposition();
  EigVec<U> res = decomp.solve(new_vec);
  std::vector<double> result(res.data(), res.data() + res.rows() * res.cols());
  // put intercept at the start if it is used
  if (use_intercept) {
    auto intercept = result[result.size() - 1];
    result.resize(result.size() - 1);
    std::reverse(result.begin(), result.end());
    result.push_back(intercept);
    std::reverse(result.begin(), result.end());
  }

  lm_coef<U> final(result, use_intercept);
  return final;
}

// map 2 vectors to Eigen matrices and call solve
template <typename U = double>
lm_coef<U> xreg_coef(std::vector<U> &y, std::vector<std::vector<U>> &xreg,
                     bool use_intercept = true) {

  const int n = y.size();
  const int ncol = xreg.size();
  auto _xreg = flatten_vec(xreg);
  EigMat<U> new_mat = Eigen::Map<EigMat<U>>(_xreg.data(), n, ncol);
  if (use_intercept) {
    new_mat.conservativeResize(Eigen::NoChange, new_mat.cols() + 1);
    new_mat.col(new_mat.cols() - 1) = EigVec<U>::Constant(n, 1, 1);
  }
  EigVec<U> new_vec = Eigen::Map<EigVec<U>>(y.data(), n, 1);

  auto decomp = new_mat.completeOrthogonalDecomposition();
  EigVec<U> res = decomp.solve(new_vec);
  std::vector<double> result(res.data(), res.data() + res.rows() * res.cols());
  // put intercept at the start if it is used
  if (use_intercept) {
    auto intercept = result[result.size() - 1];
    result.resize(result.size() - 1);
    std::reverse(result.begin(), result.end());
    result.push_back(intercept);
    std::reverse(result.begin(), result.end());
  }

  lm_coef<U> final(result, use_intercept);
  return final;
}

template<typename U=double> std::vector<U> predict(
  const int n,
  lm_coef<U> coef,
  std::vector<std::vector<U>> xreg ) {

  std::vector<U> _xreg = flatten_vec(xreg);
  EigMat<U> new_mat = Eigen::Map<EigMat<U>>(_xreg.data(), n, xreg.size());
  if (coef.has_intercept()) {
    new_mat.conservativeResize(Eigen::NoChange, new_mat.cols() + 1);
    new_mat.col(new_mat.cols() - 1) = EigVec<U>::Constant(n, 1, 1);
  }
  EigVec<U> _coef = Eigen::Map<EigVec<U>>( coef.data().data(), coef.size(), 1);
  EigVec<U> res = new_mat * _coef;
  std::vector<U> result(res.size());
  for(int i=0; i< result.size(); i++) {
    result[i] = res[i];
  }
  return result;
}

template < typename U = double> EigMat<U> vec_to_mat(
  const std::vector<std::vector<U>> &xreg,
  const size_t n,
  const bool intercept
  ) {

  std::vector<double> _xreg = flatten_vec(xreg);
  EigMat<U> new_mat = Eigen::Map<EigMat<U>>( _xreg.data(), n, xreg.size());
  if (intercept) {
    new_mat.conservativeResize(Eigen::NoChange, new_mat.cols() + 1);
    new_mat.col(new_mat.cols() - 1) = EigMat<U>::Constant(n, 1, 1);
  }
  return new_mat;
}


#endif
