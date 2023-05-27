#ifndef BLAZE_AUTO_AR
#define BLAZE_AUTO_AR

// #include "ar/ar.h"
#include "utils/xreg.h"
#include "common/forecast_result.h"
#include "ar/ar_forecaster.h"

#include "utils/stopwatch.h"

enum AutoARMethod {
  AIC,
  AICc,
  BIC
};

template <typename scalar_t> using EigVec = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;
template <typename scalar_t> using EigMat = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;

template <typename scalar_t> scalar_t fit_ar_spec(
    std::vector<scalar_t> &coef,
    const std::vector<scalar_t> &y,
    EigVec<scalar_t> &fitted, 
    EigMat<scalar_t> &xregs, 
    const size_t p) {
  using EVec  = EigVec<scalar_t>;
  using EMat  = EigMat<scalar_t>;
  
  const size_t n = y.size();
  auto lags = get_lags(y, p);
  // recover y with lags reflected 
  EVec y_ = Eigen::Map<EVec>(lags[0].data(), lags[0].size(), 1);
  // get the lags as an eigen matrix
  EMat lag_mat = EMat::Zero(n-p, p);
  for (size_t i = 1; i < lags.size(); i++) {
    lag_mat.col(i-1) += Eigen::Map<EVec>(lags[i].data(), n-p, 1);
  }
  // concatenate all matrices together
  EMat all_features(lag_mat.rows(), lag_mat.cols() + xregs.cols());
  // reflect p offset in xreg
  all_features << lag_mat, xregs.bottomRows(n-p);
  // estimate parameters
  coef.resize(p);
  coef = solve_system<1000>(all_features, y_, fitted);
  scalar_t ssq = 0.0;
  // // only values from p up are valid and should be accounted for in the residual
  // // sum of squares
  for(size_t i = p; i < n; i++) {
    ssq += std::pow(y[i] - fitted[i-p],2);
  }
  return ssq/(n-p);
}

template <typename scalar_t, typename Scaler = StandardScaler<scalar_t>> class AutoAR {
private:
  // data
  std::vector<scalar_t> y;
  size_t min_p, max_p, p;
  bool intercept, drift;
  AutoARMethod method;
  std::vector<Scaler> scalers;
  // 
  std::vector<scalar_t> residuals, prev_y, coef, raw_coef, fitted_vals;
  //
  EigMat<scalar_t> xregs;
  EigVec<scalar_t> eig_fitted_vals;
  scalar_t sigma2, aic, bic, aicc;
  bool fitted;
public:
  AutoAR<scalar_t, Scaler>(){};
  AutoAR<scalar_t, Scaler>(
      const std::vector<scalar_t> &y,
      const size_t min_p,
      const size_t max_p,
      const std::vector<std::vector<scalar_t>> &xreg = {{}},
      const bool intercept = true, const bool drift = true,
      const bool standardize = true, 
      AutoARMethod method = AIC) : y(y), min_p(min_p), max_p(max_p), p(min_p),
      intercept(intercept), drift(drift), method(method),
      scalers(std::vector<Scaler>(standardize * (1 + xreg.size()))) {
    //
    this->residuals = std::vector<scalar_t>(y.size(), 0.0);    
    this->coef = std::vector<scalar_t>(
      xreg.size() + this->max_p + this->intercept + this->drift, 0.0);
    this->raw_coef = std::vector<scalar_t>(this->coef.size());
    this->fitted_vals = std::vector<scalar_t>(this->y.size(), 0.0);
    //
    // take xreg copy so we do not modify user data
    auto _xreg = xreg;
    if( this->scalers.size() > 0 ) {
      // first scaler used for target
      this->scalers[0] = Scaler(this->y);
      this->scalers[0].scale(this->y);
      size_t i = 1;
      while(i < this->scalers.size()) {
        this->scalers[i] = Scaler(_xreg[i-1]);
        this->scalers[i].scale(_xreg[i-1]);
        i++;
      }
    }
    //
    this->xregs = vec_to_mat(_xreg, y.size(), intercept, drift);
    this->eig_fitted_vals = EigVec<scalar_t>::Zero(y.size());
    // 
    this->fitted = false;
  };
  void fit() {
    // if we are asked to fit a model with p,q, P,Q == 0, we simply throw 
    if(max_p == 0 || min_p > max_p) {
      std::cout << "You specified a model with no parameters (max_p equal" <<
        " to zero) or a model where min_p > max_p - we will not fit a model " <<
          " in this case. Please use the .respecify() method and provide  " <<
            "valid 'min_p' and 'max_p' parameters." <<
              std::endl;
      return;
    }
    if(this->method == AutoARMethod::AIC) this->fit_w_ic<AutoARMethod::AIC>();
    else if(this->method == AutoARMethod::BIC) this->fit_w_ic<AutoARMethod::BIC>();
    else if(this->method == AutoARMethod::AICc) this->fit_w_ic<AutoARMethod::AICc>();
    else {
      std::cout << "Invalid fitting method - please use respecify to provide " <<
        "a valid method, ideally one of " << 
        "AIC, AICc, BIC" << std::endl;
    }
    // compute residuals 
    for(size_t i = p; i < this->y.size(); i++) {
      residuals[i] = y[i] - fitted_vals[i];
    }
    this->prev_y = std::vector<scalar_t>(p+1);
    for(size_t i = 0; i < this->prev_y.size(); i++) {
      const size_t index = this->y.size() - this->prev_y.size() + i;
      this->prev_y[i] = this->y[index];
    }
    //
    if( this->scalers.size() > 0 ) {
      // first scaler used for target
      this->scalers[0].rescale(this->y);
      // write back raw coefficients - we will need them later
      this->raw_coef = this->coef;
      for(size_t i = 1; i < this->scalers.size(); i++) {
        scalar_t temp = this->scalers[i].rescale_val(coef[this->p+i-1]);
        this->coef[p+i-1] = temp;
      }
      if( this->intercept ) {
        // the intercept also has to be rescaled to have any meaning
        scalar_t temp = scalers[0].rescale_val_w_sd(this->coef.back());
        this->coef.back() = temp;
      }
      // rescale fitted values 
      this->scalers[0].rescale(this->fitted_vals);
      // rescale residuals to be on the original scale
      this->scalers[0].rescale_w_sd(this->residuals);
    }
    this->sigma2 = sum_of_squares(this->residuals)/(this->y.size()-this->p);
    // recompute all information criteria since this is not particularly heavy
    compute_all_ic(coef.size());
    this->fitted = true;
  }
  template <const AutoARMethod method> void fit_w_ic() {
    scalar_t best_ic = 1e12; 
    scalar_t current_ic = best_ic;
    std::vector<scalar_t> best_coef(coef.size());
    for(size_t p_ = this->max_p; p_ >= this->min_p; --p_) {
      this->eig_fitted_vals.resize(this->y.size()-p_);
      // create candidate model
      this->sigma2 = fit_ar_spec(coef, this->y, this->eig_fitted_vals,
                                 this->xregs, p_);
      if constexpr(method == AutoARMethod::AIC) {
        current_ic = compute_aic(p_);
      }
      if constexpr(method == AutoARMethod::BIC) {
        current_ic = compute_bic(p_);
      }
      if constexpr(method == AutoARMethod::AICc) {
        current_ic = compute_aicc(p_);
      }
      // check aic
      if(current_ic < best_ic) {
        best_coef.resize(p_);
        // if best, set the fitted model to candidate 
        best_coef = this->coef;
        // and update best aic
        best_ic = current_ic;
        this->p = p_;
        // write back fitted values
        for(size_t i = p_; i < this->y.size(); i++) {
          this->fitted_vals[i] = this->eig_fitted_vals[i-p_];
        }
      }
    }
    this->coef.resize(this->p);
    this->coef = best_coef;
  };
  const scalar_t compute_aic(const size_t p) {
    return this->y.size() * std::log(this->sigma2) + (2 * this->coef.size());  
  }
  const scalar_t compute_bic(const size_t p) {
    // these need to be verified (they might be slightly off)
    return this->y.size() * std::log(this->sigma2) + 
          (2 * this->coef.size()) + (p + 1) * (log(this->y.size()) - 2);
  }
  const scalar_t compute_bic(const scalar_t aic, const size_t p) {
    return aic + (p + 1) * (log(this->y.size()) - 2);
  }
  const scalar_t compute_aicc(const size_t p) {
    return this->y.size() * std::log(this->sigma2) + 
          (2 * this->coef.size()) + (2*(std::pow(p, 2) +
          this->p)/(this->y.size() - p - 1)); 
  }
  const scalar_t compute_aicc(const scalar_t aic, const size_t p) {
    return aic + (2*(std::pow(p, 2) + p)/(this->y.size() - p - 1)); 
  }
  void compute_all_ic(const size_t p) {
    this->aic = compute_aic(p);
    this->bic = compute_bic(this->aic, p);
    this->aicc = compute_aicc(this->aic, p);
  }
  void respecify(const size_t min_p, const size_t max_p) {
    this->min_p = min_p;
    this->max_p = max_p;
  }
  void respecify(AutoARMethod method) {
    this->method = method;
  }
  forecast_result<scalar_t> forecast(
      const size_t h = 10,
      std::vector<std::vector<scalar_t>> &newxreg = {{}}) {
    // validate xreg length
    if (!this->fitted || newxreg.size() != (this->xregs.cols() - intercept - drift))
      return forecast_result<scalar_t>(0);
    // otherwise run forecast
    // scale xreg as needed
    if(newxreg.size() > 0 || this->intercept) {
      if(this->scalers.size() > 0) {
        size_t i = 1;
        for(auto & xreg_val:newxreg) {
          this->scalers[i].scale(xreg_val);
          i++;
        }
      }
    }
    // run actual forecast
    auto res = forecast_ar(h, newxreg, this->prev_y, this->raw_coef, this->sigma2,
                           this->p, this->intercept, this->drift,
                           this->y.size());
    // scale by sigma2 and take square root
    for(size_t i = 0; i < h; i++) {
      res.std_err[i] = sqrt(res.std_err[i] * this->sigma2);
    }
    // rescale forecasts
    if(this->scalers.size() > 0) {
      this->scalers[0].rescale(res.forecast);
      this->scalers[0].rescale_w_sd(res.std_err);
    }
    return res;
  };
  const std::vector<scalar_t> get_coef() const { return this->coef; }
  const std::vector<scalar_t> get_residuals() const { return this->residuals; }
  const std::vector<scalar_t> get_fitted() const { return this->fitted_vals; }
  const scalar_t get_sigma2() const { return this->sigma2; }
  const scalar_t get_aic() const { return this->aic; }
  const scalar_t get_bic() const { return this->bic; }
  const scalar_t get_aicc() const { return this->aicc; }
  const size_t get_order() const { return this->p;}
  const bool is_fitted() const { return this->fitted; }
};

#endif
