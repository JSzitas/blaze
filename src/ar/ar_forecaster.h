#ifndef AR_FORECAST
#define AR_FORECAST

#include "common/forecast_result.h"
#include "third_party/eigen.h"
#include "utils/utils.h"
#include "utils/xreg.h"

template <typename scalar_t> std::vector<scalar_t>
ar_std_err(const std::vector<scalar_t> &coef,
           const size_t p,
           const size_t h,
           const scalar_t sigma2) {
  std::vector<scalar_t> psi(h + p + 1 + 1, 0.0);
  psi[0] = 1.0;
  for(size_t i = 1; i < p+1; i++) psi[i] = coef[i-1];
  for(size_t i = 0; i < h; i++) {
    for (size_t j = 0; j < p; j++) {
      psi[i + j + 1 + 1] += coef[j] * psi[i+1];
    }
  } 
  psi.resize(h);
  // square all elements in psi
  for( auto& elem:psi ) elem *= elem;
  // run cummulative sum
  for(size_t j = 1; j < h; j++) psi[j] += psi[j-1];
  return psi;
}

template <typename scalar_t> forecast_result<scalar_t>
forecast_ar( const size_t h,
             std::vector<std::vector<scalar_t>> &newxreg, 
             const std::vector<scalar_t> &prev_y, 
             std::vector<scalar_t> &coef,
             const scalar_t sigma2,
             const size_t p,
             const bool intercept = true, 
             const bool drift = false,
             const size_t drift_offset = 1) {

  using EigVec = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;
  using EigMat = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;
  using EigRow = Eigen::Matrix<scalar_t, 1, Eigen::Dynamic>;
    
  std::vector<scalar_t> _xreg = flatten_vec(newxreg);
  EigMat xreg_mat = Eigen::Map<EigMat>(_xreg.data(), h, newxreg.size());
  if (drift) add_drift(xreg_mat, drift_offset);
  if (intercept) add_constant(xreg_mat);
  // add lags
  EigRow lag_mat = EigRow::Zero(p);
  for(size_t i = 0; i < p; i++) {
    lag_mat(i) = prev_y[prev_y.size()-1-i];
  }
  // set up coefficients   
  EigVec _coef = Eigen::Map<EigVec>(coef.data(), coef.size(), 1);
  EigVec result = EigVec::Zero(h);
    // prediction recursion 
  for(size_t i = 0; i < h; i++) {
    result[i] = lag_mat * _coef.head(p);
    result[i] += xreg_mat.row(i) * _coef.tail(_coef.size()-p);
    // update lag_mat 
    for( size_t j = 0; j < p; j++ ) {
      lag_mat(p-j-1) = lag_mat(p-j-2);
    }
    lag_mat(0) = result[i];
  }
  std::vector<scalar_t> means(result.size());
  for(size_t i=0; i< result.size(); i++) means[i] = result[i];
  std::vector<scalar_t> std_errs = ar_std_err(coef, p, h, sigma2);
  return forecast_result<scalar_t>(means, std_errs);
}

#endif
