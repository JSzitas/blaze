#ifndef AR_FITTER
#define AR_FITTER

#include "third_party/eigen.h"

#include "utils/utils.h"
#include "utils/xreg.h"

template <typename scalar_t> scalar_t fit_ar(
    std::vector<scalar_t> & coef,
    std::vector<scalar_t> & fitted, 
    std::vector<scalar_t> & residuals,
    const std::vector<scalar_t> &y,
    const std::vector<std::vector<scalar_t>> & xreg,
    const size_t p,
    const bool intercept,
    const bool drift) {
 
 using EigVec = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;
 using EigMat = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;
 
 const size_t n = y.size();
 
 EigMat xregs = vec_to_mat(xreg, n, intercept, drift);
 auto lags = get_lags(y, p);
 // recover y with lags reflected 
 EigVec y_ = Eigen::Map<EigVec>(lags[0].data(), lags[0].size(), 1);
 // get the lags as an eigen matrix
 EigMat lag_mat = EigMat::Zero(n-p, p);
 for (size_t i = 1; i < lags.size(); i++) {
   lag_mat.col(i-1) += Eigen::Map<EigVec>(lags[i].data(), n-p, 1);
 }
 // reflect this in xreg
 xregs = xregs.bottomRows(n-p);
 // concatenate all matrices together
 EigMat all_features(lag_mat.rows(), lag_mat.cols() + xregs.cols());
 all_features << lag_mat, xregs;
 // estimate parameters
 coef = solve_ortho_decomp(all_features, y_);
 // print_vector(coef);
 auto eig_coef = Eigen::Map<EigVec>(coef.data(), coef.size(), 1);
 // compute in-sample fit 
 auto fitted_vals = all_features * eig_coef;
 for(size_t i = p; i < n; i++) {
   fitted[i] = fitted_vals[i-p];
 }
 scalar_t ssq = 0.0;
 // only values from p up are valid and should be accounted for in the residual
 // sum of squares
 for(size_t i = p; i < n; i++) {
   scalar_t resid = y[i] - fitted[i];
   residuals[i] = resid;
   ssq += resid * resid;
 }
 return ssq/(n-p);
}

#endif
