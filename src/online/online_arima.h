#ifndef BLAZE_ONLINE_ARIMA
#define BLAZE_ONLINE_ARIMA

#include <cmath>

#include "common/forecast_result.h"
#include "structures/circulant.h"
#include "utils/utils.h"
#include "utils/stopwatch.h"

#include "utils/sherman_morisson.h"


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

enum OArimaType{
  OGD,
  OND,
  MIX
};

template <typename scalar_t,
          const OArimaType type = OArimaType::MIX,
          const size_t max_mk = 20> class OnlineARIMA {
private:
  // data
  circulant<scalar_t> y, diff_y;
  // OND data + OGD data
  std::vector<scalar_t> coef_ond = {}, A_inv, coef_ogd = {};
  // parameters 
  size_t mk, d, n_samples = 0;
  scalar_t eta;
  // performance stats 
  scalar_t sigma2, aic, bic, aicc;
  // scaling variables
  std::array<scalar_t,2> scaling = {0.0, 0.0};
public:
  OnlineARIMA<scalar_t, type, max_mk>(){}
  OnlineARIMA<scalar_t, type, max_mk>(
      const size_t mk = 10,
      const size_t d = 1,
      const scalar_t eta = 0.05) : mk(mk), d(d), eta(eta), 
      y(circulant<scalar_t>(mk)), diff_y(circulant<scalar_t>(mk*(d>0))) {
    // coef ordered from oldest to newest (i.e. opposite to the common 
    // ordering of AR_1 then AR_2 till AR_p)
    // initialize accordingly - the oldest coefficients should be smallest 
    if constexpr(type == OArimaType::OND || type == OArimaType::MIX) {
      this->coef_ond = std::vector<scalar_t>(this->mk, 1.0);
      for(size_t i = 0; i < this->mk; i++) coef_ond[i] /= (this->mk - i);
      // A_inv keeps the inverse matrix
      this->A_inv = std::vector<scalar_t>(this->mk * this->mk, 0.0);
      // initialize A_inv to identity matrix with diagonal multiplied by eps
      for(size_t i = 0; i < this->mk; i++) {
        this->A_inv[i * this->mk + i] = 1.0;
      }
    }
    if constexpr(type == OArimaType::OGD || type == OArimaType::MIX) {
      this->coef_ogd = std::vector<scalar_t>(this->mk, 1.0);
      for(size_t i = 0; i < this->mk; i++) coef_ogd[i] /= (this->mk - i);
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
    if(this->n_samples < (this->mk + 2 + this->d)) {
      if(this->n_samples > 2) {
        // at least record samples
        if(this->d == 1) this->diff_y.push_back(current_y);
        this->y.push_back(current_y);
      }
      return;
    }
    // otherwise run live differencing and  compute current prediction
    if( this->d == 1) {
      // difference by taking current_y vs last circulant element
      current_y -= this->y.back();
    } 
    if constexpr(type == OArimaType::OND || type == OArimaType::MIX) {
      ond_update(-current_y);
    }
    if constexpr(type == OArimaType::OGD || type == OArimaType::MIX) {
      ogd_update(-current_y);
    }
    // update last circulant element 
    if(this->d == 1) this->diff_y.push_back(current_y);
    this->y.push_back(current_y);
  }
private:
  void ond_update(const scalar_t current_y) {
    scalar_t grad_tmp = current_y;
    if(this->d== 1) {
      for(size_t i = 0; i < this->mk; i++) {
        grad_tmp += this->diff_y[i] * this->coef_ond[i];
      }
    } else {
      for(size_t i = 0; i < this->mk; i++) {
        grad_tmp += this->y[i] * this->coef_ond[i];
      }
    }
    // this->y.print();
    // print_vector(this->coef_ond);
    // std::cout << "Current y: " << current_y << " after filtering: " << grad_tmp << std::endl;
    this->sigma2 += std::pow(grad_tmp, 2);
    // gradient 
    std::array<scalar_t, max_mk> g;
    for(size_t i = 0; i < this->mk; i++) {
      g[i] = 2 * grad_tmp * this->y[i];
    }
    // std::cout << "Gradient: ";
    // print_vector(g, 0, this->mk);
    // carry out sherman morrison
    // print_vector(this->A_inv);
    sherman_morrison_update(this->A_inv, g);
    // print_vector(this->A_inv);
    // update coefficients - this must be done with a temporary, because 
    // we need to retain the coefficients while carrying out all of the vector
    // matrix multiplications (duh)
    for(size_t j = 0; j < this->mk; j++) {
      scalar_t tmp = 0.0;
      for(size_t i = 0; i < this->mk; i++) {
        // note that this is coefficient * row, for a given column of the matrix
        // hence our offset in j is constant
        tmp += g[i] * this->A_inv[(this->mk*i) + j];
      }
      this->coef_ond[j] -= this->eta * tmp; 
    }
    // print_vector(this->coef_ond);
  }
  void ogd_update(const scalar_t current_y) {
    scalar_t tmp = current_y;
    if(this->d == 1) {
      for(size_t i = 0; i < this->mk; i++) {
        tmp +=  this->diff_y[i] * this->coef_ogd[i];
      }
    } else {
      for(size_t i = 0; i < this->mk; i++) {
        tmp += this->y[i] * this->coef_ogd[i];
      }
    }
    // for mix this probably needs to be divided by 2 at the end
    this->sigma2 += std::pow(tmp, 2);
    // using the 'cooling schedule' in authors original code
    const scalar_t step = eta/sqrt((this->n_samples+1) - this->mk);
    // direct update - does not actually require grad 
    // we can further precompute the size of step * (2*tmp)
    // where grad would be (2* tmp * this->y[i])
    const scalar_t update = 2 * step * tmp;
    for(size_t i = 0; i < mk; i++) {
      this->coef_ogd[i] -= update * this->y[i];;
    }
  }
public:
  void print_summary() {
    std::cout << "Estimated coefficients (high order to low order): " << std::endl;
    if constexpr(type != OArimaType::OND) {
      std::cout << "OGD coef: ";
      print_vector(this->coef_ogd);
    }
    if constexpr(type != OArimaType::OGD) {
      std::cout << "OND coef: ";
      print_vector(this->coef_ond);
    }
    std::cout << "Scaling mean: " << this->scaling[0] << " scaling sd: " << 
      std::sqrt(this->scaling[1]/this->n_samples) << std::endl;
  }
  std::tuple<scalar_t, scalar_t> forecast_step() {
    if constexpr(type == OArimaType::MIX) {
      auto ogn = forecast_step_ogd();
      auto odn = forecast_step_ond();
      // std::cout << "ogn: " << std::get<0>(ogn) << " ond: " << std::get<0>(odn) <<
      //   " avg:" << (std::get<0>(ogn) + std::get<0>(odn))/2 << std::endl;
      return std::tuple<scalar_t, scalar_t>(
          (std::get<0>(ogn) + std::get<0>(odn))/2,
          (std::get<1>(ogn) + std::get<1>(odn))/2);
    }
    if constexpr(type == OArimaType::OGD) {
      return forecast_step_ogd();
    }
    if constexpr(type == OArimaType::OND) {
      return forecast_step_ond();
    }
  }
  forecast_result<scalar_t> forecast(const size_t h = 10) {
    if constexpr(type == OArimaType::MIX) {
      auto ogn = forecast_ogd(h);
      auto odn = forecast_ond(h);
      for(size_t i = 0; i < h; i++) {
        ogn.forecast[i] += odn.forecast[i];
        ogn.forecast[i] /= 2;
        
        ogn.std_err[i] += odn.std_err[i];
        ogn.std_err[i] /= 2;
      }
      return ogn;
    }
    if constexpr(type == OArimaType::OGD) return forecast_ogd(h);
    if constexpr(type == OArimaType::OND) return forecast_ond(h);
  }
private:
  std::tuple<scalar_t, scalar_t> forecast_step_ond() {
    if (this->n_samples < (this->mk + this->d))
      return std::tuple<scalar_t, scalar_t>(0.0, 1.0);
    scalar_t mean_fcst = 0.0;
    for(size_t i = 0; i < this->mk; i++) mean_fcst += this->y[i] * this->coef_ond[i];
    scalar_t eff_sigma;
    // type MIX gets averaging 
    if constexpr(type == OArimaType::MIX) eff_sigma = std::sqrt(this->sigma2/2);
    else eff_sigma = std::sqrt(this->sigma2);
    return std::tuple<scalar_t, scalar_t>(
        // invert scaling 
        (mean_fcst * std::sqrt(this->scaling[1]/this->n_samples)) + 
          this->scaling[0], eff_sigma);
  }
  std::tuple<scalar_t, scalar_t> forecast_step_ogd() {
    if (this->n_samples < (this->mk + this->d))
      return std::tuple<scalar_t, scalar_t>(0.0, 1.0);
    scalar_t mean_fcst = 0.0;
    for(size_t i = 0; i < this->mk; i++) mean_fcst += this->y[i] * this->coef_ogd[i];
    scalar_t eff_sigma;
    // type MIX gets averaging 
    if constexpr(type == OArimaType::MIX) eff_sigma = std::sqrt(this->sigma2/2);
    else eff_sigma = std::sqrt(this->sigma2);
    return std::tuple<scalar_t, scalar_t>(
        // invert scaling 
        (mean_fcst * std::sqrt(this->scaling[1]/this->n_samples)) + 
          this->scaling[0], eff_sigma);
  }
  // TODO: This has a bug that causes random segfaults, I think
  forecast_result<scalar_t> forecast_ogd(const size_t h = 10) {
    // validate xreg length
    if (this->n_samples < (this->mk + this->d))
      return forecast_result<scalar_t>(0);
    // run actual forecast recursion
    std::vector<scalar_t> mean_fcst(h, 0.0);
    // set up a temporary circulant to make the code easier to read
    fixed_circulant<scalar_t, max_mk> temp = this->y;
    // for now pretend d is always 0 to simplify this
    for(size_t i = 0; i < h; i++) {
      scalar_t tmp = 0.0;
      // size_t offset = i > this->mk ? i - this->mk : 0;
      for(size_t j = 0; j < this->mk; j++) {
        tmp += temp[j] * this->coef_ogd[j];
      }
      temp.push_back(tmp);
      // update circulant
      mean_fcst[i] = tmp;
    }
    auto std_errs = ar_std_err(this->coef_ogd, h, this->sigma2);
    for(size_t i = 0; i < h; i++) {
      std_errs[i] = sqrt(std_errs[i] * this->sigma2);
    }
    // rescale
    const scalar_t sigma_rescl = std::sqrt(this->scaling[1]/this->n_samples);
    for(auto &val:mean_fcst) {
      val *= sigma_rescl;
      val += this->scaling[0];
    }
    return forecast_result<scalar_t>(mean_fcst, std_errs);
  };
  forecast_result<scalar_t> forecast_ond(const size_t h = 10) {
    // validate xreg length
    if (this->n_samples < (this->mk + this->d))
      return forecast_result<scalar_t>(0);
    std::vector<scalar_t> mean_fcst(h, 0.0);
    fixed_circulant<scalar_t, max_mk> temp = this->y;
    // for now pretend d is always 0 to simplify this
    for(size_t i = 0; i < h; i++) {
      scalar_t tmp = 0.0;
      // size_t offset = i > this->mk ? i - this->mk : 0;
      for(size_t j = 0; j < this->mk; j++) {
        tmp += temp[j] * this->coef_ond[j];
      }
      temp.push_back(tmp);
      // update circulant
      mean_fcst[i] = tmp;
    }
    auto std_errs = ar_std_err(this->coef_ond, h, this->sigma2);
    for(size_t i = 0; i < h; i++) {
      std_errs[i] = sqrt(std_errs[i] * this->sigma2);
    }
    // rescale
    const scalar_t sigma_rescl = std::sqrt(this->scaling[1]/this->n_samples);
    for(auto &val:mean_fcst) {
      val *= sigma_rescl;
      val += this->scaling[0];
    }
    return forecast_result<scalar_t>(mean_fcst, std_errs);
  };
public:
  const std::vector<scalar_t> get_coef() const { 
    if constexpr(type == OArimaType::OGD) {
      return this->coef_ogd;
    }
    if constexpr(type == OArimaType::OND) {
      return this->coef_ond;
    }
    std::vector<scalar_t> coef(2 * this->mk);
    for(size_t i=0; i < this->mk; i++) {
      coef[i] = this->coef_ogd[i];
      coef[i+this->mk] = this->coef_ond[i];
    }
    return coef;
  }
  const scalar_t get_sigma2() const {
    return this->sigma2/(this->n_samples - this->mk);
  }
};

#endif
