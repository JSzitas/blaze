#ifndef BLAZE_ONLINE_ARIMA
#define BLAZE_ONLINE_ARIMA

#include <cmath>

#include "common/forecast_result.h"
#include "structures/circulant.h"
#include "utils/utils.h"
#include "utils/stopwatch.h"

template <typename scalar_t> std::vector<scalar_t>
ar_std_err(const std::vector<scalar_t> &coef,
           const size_t h,
           const scalar_t sigma2) {
  const size_t p = coef.size();
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

// note that this can probably be made even more efficient - matrix A 
// is symmetric, so in theory only the upper/lower triangular part need be computed
template <typename scalar_t> void 
  sherman_morrison_update(
    std::vector<scalar_t>& A, 
    std::vector<scalar_t>& g) {
    // update an inverse matrix A using updated gradient g
    const size_t p = g.size();
    for(size_t i = 0; i < p; i++) {
      // this is so we avoid repeated dereferences - this should be a very cheap, 
      // as the temporary is created by dereferencing once, and is put on the stack
      // rather than dereferencing many times
      const scalar_t g_i = g[i];
      for(size_t j = 0; j < p; j++) {
        // see above
        const scalar_t a = A[(i*p) + j], g_j = g[j];
        A[(i*p) + j] -= (a * g_i * g_j * a)/(1 + g_i * a * g_j);
      }
    }
}

enum OArimaType{
  OGD,
  OND,
  // mixture is equal weight average of both OGD and OND
  MIXTURE
};

template <typename scalar_t,
          const OArimaType type = OArimaType::OND> class OnlineARIMA {
private:
  // data
  circulant<scalar_t> y, diff_y;
  // coef_temp and A_inv only necessary for OGD
  std::vector<scalar_t> coef, grad, coef_temp, A_inv;
  // parameters 
  size_t mk, d, n_samples = 0;
  scalar_t eta;
  // performance stats 
  scalar_t sigma2, aic, bic, aicc;
  // scaling variables
  std::array<scalar_t,2> scaling = {0.0, 0.0};
public:
  OnlineARIMA<scalar_t, type>(){}
  OnlineARIMA<scalar_t, type>(
      const size_t mk = 10,
      const size_t d = 1,
      const scalar_t eta = 0.03,
      const scalar_t eps = 0.5) : mk(mk), d(d), eta(eta),
      y(circulant<scalar_t>(mk)), diff_y(circulant<scalar_t>(mk*(d>0))) {
    //
    this->n_samples = 0;
    // coef ordered from oldest to newest (i.e. opposite to the common 
    // ordering of AR_1 then AR_2 till AR_p)
    this->coef = std::vector<scalar_t>(this->mk, 1.0);
    // initialize accordingly - the oldest coefficients should be smallest 
    for(size_t i = 0; i < mk; i++) {
      this->coef[i] /= mk-i;
    }
    this->grad = std::vector<scalar_t>(this->mk);
    if constexpr(type == OArimaType::OND) {
      this->coef_temp = std::vector<scalar_t>(this->mk, 0.0);
      this->A_inv = std::vector<scalar_t>(this->mk * this->mk, 0.0);
      // initialize A_inv to identity matrix
      for(size_t i = 0; i < this->mk; i++) {
        this->A_inv[(i * this->mk) + i] = eps;
      }
    }
  };
  void update(scalar_t current_y) {
    if(this->n_samples == 0) {
      this->scaling[0] = current_y;
    }
    else {
      // otherwise run welfords algorithm to reflect scaling
      welfords_algorithm(this->scaling, current_y, this->n_samples);
      // and scale
      current_y = (current_y - this->scaling[0])/std::sqrt(this->scaling[1]/this->n_samples);
    }
    this->n_samples++;
    // if we are still accumulating samples, only do this
    if(this->n_samples < (mk + this->d)) {
      return;
    }
    scalar_t tmp = 0.0;
    // otherwise run live differencing and  compute current prediction
    if( this->d == 1) {
      // difference by taking current_y vs last circulant element 
      tmp = current_y - this->y.back();
      for(size_t i = 0; i < this->mk; i++) {
        tmp -= this->diff_y[i] * this->coef[i];
      }
      // update circulant 
      this->diff_y.push_back(tmp);
    } else {
      tmp = -current_y;
      for(size_t i = 0; i < this->mk; i++) {
        tmp += this->y[i] * this->coef[i];
      }
    }
    // TODO:
    // update running 'sigma2' (this is not normalized by the number of observations!!!)
    // nor by scale - so it needs to be 'fixed'
    this->sigma2 += std::pow(tmp, 2);
    // now the REAL fun - compute gradient (actually pretty easy)
    for(size_t i = 0; i < mk; i++) {
      this->grad[i] = 2 * this->y[i] * tmp;
    }
    if constexpr(type == OArimaType::OND) {
      // update inverse using the Sherman-Morrison formula  
      sherman_morrison_update(this->A_inv, this->grad);
      // update coefficients - this must be done with a temporary, because 
      // we need to retain the coefficients while carrying out all of the vector
      // matrix multiplications (duh)
      for(size_t j = 0; j < this->mk; j++) {
        scalar_t tmp = 0.0;
        for(size_t i = 0; i < this->mk; i++) {
          // note that this is coefficient * row, for a given column of the matrix
          // hence our offset in j is constant
          tmp += this->grad[i] * this->A_inv[(this->mk*j) + i];
        }
        this->coef_temp[j] = this->eta * tmp; 
      }
      // update coefficients from coefficient update temporary
      for(size_t j = 0; j < this->mk; j++) {
        this->coef[j] -= this->coef_temp[j];
      }
    }
    else if constexpr(type == OArimaType::OGD) {
      // using the 'cooling schedule' in authors original code
      scalar_t step = eta/sqrt((this->n_samples+1) - this->mk);
      // only uses gradient to update 
      for(size_t i = 0; i < this->mk; i++) {
        this->coef[i] -= step * this->grad[i];
      }
    }
    // update last circulant element 
    this->y.push_back(current_y);
  }
  // TODO: Fix these 
  const scalar_t compute_aic() {
    return this->n_samples * std::log(
      this->sigma2/(this->n_samples - this->mk)) + (2 * coef.size());
  }
  const scalar_t compute_bic() {
    // these need to be verified (they might be slightly off)
    return this->n_samples * std::log(
      this->sigma2/(this->n_samples - this->mk)) +
      (2 * this->mk) + (this->mk + 1) * (log(this->n_samples) - 2);
  }
  const scalar_t compute_aicc() {
    return this->n_samples * std::log(
      this->sigma2/(this->n_samples - this->mk)) +
      (2 * this->mk) + (2*(std::pow(this->mk, 2) +
      this->mk)/(this->n_samples - this->mk - 1));
  }
  std::tuple<scalar_t, scalar_t> forecast_step() {
    if (this->n_samples < (this->mk + this->d))
      return std::tuple<scalar_t, scalar_t>(0.0, 1.0);
    scalar_t mean_fcst = 0;
    for(size_t i = 0; i < this->mk; i++) {
      mean_fcst += this->y[i] * this->coef[i];
    }
    return std::tuple<scalar_t, scalar_t>(
        // invert scaling 
        (mean_fcst * std::sqrt(this->scaling[1]/this->n_samples)) + this->scaling[0],
        std::sqrt(this->sigma2));
  }
  void print_summary() {
    std::cout << "Estimated coefficients (high order to low order): " << std::endl;
    print_vector(this->coef);
    std::cout << "Scaling mean: " << this->scaling[0] << " scaling sd: " << 
      std::sqrt(this->scaling[1]/this->n_samples) << std::endl;
  }
  // TODO: This has a bug that causes random segfaults, I think
  forecast_result<scalar_t> forecast(const size_t h = 10) {
    // validate xreg length
    if (this->n_samples < (this->mk + this->d))
      return forecast_result<scalar_t>(0);
    // run actual forecast recursion 
    std::vector<scalar_t> mean_fcst(h, 0.0);
    // for now pretend d is always 0 to simplify this 
    // this is actually a mixture of data and previous mean forecasts when
    // iterating forward
    for(size_t i = 0; i < h; i++) {
      scalar_t temp = 0.0;
      for(size_t j = 0; j < i; j++) {
        temp += mean_fcst[j] * this->coef[this->mk - 1 - j];
      }
      for(size_t j = i; j < this->mk-i; j++) {
        temp += this->y[j] * this->coef[j];
      }
      mean_fcst[i] = temp;
    }
    auto std_errs = ar_std_err(this->coef, h, this->sigma2);
    for(size_t i = 0; i < h; i++) {
      std_errs[i] = sqrt(std_errs[i] * this->sigma2);
    }
    return forecast_result<scalar_t>(mean_fcst, std_errs);
  };
  const std::vector<scalar_t> get_coef() const { return this->coef; }
  const scalar_t get_sigma2() const {
    return this->sigma2/(this->n_samples - this->mk);
  }
  const scalar_t get_aic() const { return compute_aic(); }
  const scalar_t get_bic() const { return compute_bic(); }
  const scalar_t get_aicc() const { return compute_aicc(); }
};

#endif
