#ifndef BLAZE_THETA
#define BLAZE_THETA

#include <utils/seasonality_utils.h>
#include <decomposition/decompose.h>
#include <common/forecast_result.h>

// this only impact the setting of the 'dynamic' flag and the parameter 
// constaints - generally STM and DSTM impose much stricter constraints on the 
// parameter bounds - thus we probably do not need these (or we might provide a 
// Theta )

enum ThetaType{
  STM, 
  DSTM, 
  OTM,
  DOTM
};

template <typename scalar_t,
          const DecomposeType seas_type = DecomposeType::A,
          const ThetaType type = ThetaType::DOTM>
class Theta{
  const std::vector<scalar_t> &y;
  std::vector<scalar_t> y_temp;
  const size_t n; 
  size_t seasonality;
  // initialized to 0
  std::vector<scalar_t> mu;
  scalar_t ell, A, B, mean_y;
  // initialized using data
  std::vector<scalar_t> season, fitted, residuals;
  scalar_t cte, Bn, An; 
public:
  Theta<scalar_t>(){};
  Theta<scalar_t>(const std::vector<scalar_t> &y,
                  size_t seasonality = 0,
                  const scalar_t ell = 0.0,
                  const scalar_t alpha = 0.5,
                  const scalar_t theta = 0.5) : y(y), y_temp(y), n(y.size()),
                  seasonality(seasonality),
                  mu(std::vector<scalar_t>(this->n, 0.0)), 
                  ell(ell), A(0.0), B(0.0), mean_y(0.0),
                  season(std::vector<scalar_t>(this->n, 0.0)),
                  fitted(std::vector<scalar_t>(this->n, 0.0)),
                  residuals(std::vector<scalar_t>(this->n, 0.0)) {
    // if seasonalityty is not None
    if constexpr( seas_type != DecomposeType::N ) {
      // detect seasonality if necessary
      if(seasonality == 0) {
        seasonality = find_seasonalities(y, 1)[0];
      }
      this->seasonality = seasonality < 4 ? seasonality : 0;
      // apply seasonal decomposition if it makes any sense
      // first get seasonal component 
      if constexpr(seas_type == DecomposeType::A) {
        this->season = decompose_seasonality<DecomposeType::A>(y, this->seasonality);
        for(size_t i = 0; i < this->n; i++) this->y_temp[i] -= this->season[i];
      }
      else if constexpr(seas_type == DecomposeType::M) {
        if(all_positive(y)) {
          this->season = decompose_seasonality<DecomposeType::M>(y, this->seasonality);
          for(size_t i = 0; i < this->n; i++) this->y_temp[i] /= this->season[i];
        } else {
          for(size_t i = 0; i < this->n; i++) this->season[i] = 1.0;
        }
      }
    }
    // otherwise we will not need seasonality
    // compute cte 
    scalar_t temp = 0.0;
    for(auto & val:this->y_temp) temp += std::abs(val);
    this->cte = temp/this->n;
    temp = 0.0;
    // replace temp with mean of y 
    for(auto & val:this->y_temp) temp += val;
    temp /= this->n;
    scalar_t temp2 = 0.0;
    for(size_t i = 0; i < this->n; i++) temp2 += this->y_temp[i] * (i+1);
    temp2 /= this->n;
    this->Bn = ((2 * temp2) - ((1 + n) * temp))/(6 * (std::pow(n,2) - 1));
    this->An = temp - ((n + 1) * this->Bn/2);
  }
  void fit() {
    // if(optimize) // STM does not need to be optimized iirc - it is just fixed 
    // if (opt.method == "Nelder-Mead") 
    //   opt = optim(par = par_ini, fn = sse, method = "Nelder-Mead")
    //   if (opt.method == "L-BFGS-B") 
    //     opt = optim(par = par_ini, fn = sse, method = "L-BFGS-B", 
    //                 lower = lower, upper = upper)
    //     if (opt.method == "SANN") 
    //       opt = optim(par = par_ini, fn = sse, method = "SANN")
    
    
    // opt = optim(par = par_ini, fn = sse, method = "L-BFGS-B", 
    //             lower = lower, upper = upper)
    // // if (opt.method == "SANN") 
    // // opt = optim(par = par_ini, fn = sse, method = "SANN")
    // // par = opt$par
    // // }
    // // else {
    // // par = par_ini
    // // }
    for(size_t i = 0; i < this->n;i ++) {
      const scalar_t temp;
      if constexpr(seas_type == DecomposeType::A) {
        temp = this->mu[i] + this->season[i];
      } else if constexpr(seas_type == DecomposeType::M) {
        temp = this->mu[i] * this->season[i];
      }
      this->residuals[i] = this->y[i] - temp;
      this->fitted[i] = temp;
    }
  }
  // just a nicer wrapper around the implementation
  forecast_result<scalar_t> forecast(const size_t h) {
    if constexpr(type == ThetaType::DOTM || type == ThetaType::DSTM) {
      return forecast_impl<true>(h);
    }
    return forecast_impl<false>(h); 
  }
private: 
  // this is here so we can provide a nicer wrapper around the dynamic parameter choice
  template <const bool dynamic = true> 
  forecast_result<scalar_t> forecast_impl(const size_t h) {
    // initialize forecast vector
    std::vector<scalar_t> forecast(h, 0.0);
    // avoid warnings - eventually we either keep these in the object, or keep 
    // a packed version of some sort
    const scalar_t alpha = 0.0, theta = 0.0;
    // last available index in training data is i-1
    size_t i = this->n;
    // note that these can be single values since they only need the last available 
    // value 
    scalar_t ell_temp = this->ell[i-1], temp_mean_y = this->mean_y[i-1], a = 0.0, b = 0.0;
    // first value in a, b
    if constexpr(dynamic) {
      b = this->B[i-1];
      a = this->A[i-1];
    } else {
      b = this->Bn;
      a = this->An;
    }
    // this is actually the entire forecast loop
    for(size_t j = 0; j < h; j++) {
      const scalar_t prediction = ell_temp + 
        (1 - 1/this->theta) * (a * std::pow(1 - alpha,i) + 
        b * (1 - std::pow(1 - alpha,i + 1))/alpha);
      ell_temp = alpha * prediction + (1 - alpha) * ell_temp;
      temp_mean_y = (i * temp_mean_y + prediction)/(i + 1);
      if constexpr(dynamic) {
        b = ((i - 1) * b + 6 * (prediction - temp_mean_y)/(i + 1))/(i + 2);
        a = temp_mean_y - b * (i + 2)/2;
      }
      forecast[j] = prediction;
      i++;
    }
    // TODO: fix SD 
    return forecast_result<scalar_t>(forecast, std::vector<scalar_t>(h, 1.0));
  }
  template <const bool dynamic = true> scalar_t
  loss(std::vector<scalar_t> &pars) {
    scalar_t ell0 = pars[0], alpha = pars[1], theta = pars[2];
    this->ell = alpha * this->y_temp[0] + (1 - alpha) * ell0;
    this->mean_y = this->y_temp[0];
    if constexpr(dynamic) {
      this->A = this->y_temp[0];
      this->B = 0;
      this->mu[0] = this->y_temp[0];
    }
    else {
      this->A = this->An;
      this->B = this->Bn;
      this->mu[0] = ell0 + (1 - 1/theta) * (this->An + this->Bn);
    }
    for (size_t i = 0; i < (this->n - 1); i++) {
      this->mu[i + 1] = ell + (1 - 1/theta) * (this->A * ((1 - alpha)^i) + 
        this->B * (1 - (1 - alpha)^(i + 1))/alpha);
      ell = alpha * this->y_temp[i + 1] + (1 - alpha) * ell;
      mean_y = (i * mean_y + this->y_temp[i + 1])/(i + 1);
      if constexpr(dynamic) {
        this->B = ((i - 1) * this->B + 6 * (this->y_temp[i + 1] - mean_y)/(i + 1))/(i + 2);
        this->A = mean_y - this->B * (i + 2)/2;
      }
    }
    scalar_t temp = 0.0;
    size_t i = 0; 
    if constexpr(dynamic) i = 2;
    for(; i < this->n; i++) {
      temp += std::pow((this->y_temp[i] - this->mu[i])/this->cte, 2);
    }
    return temp;
  }
};

#endif
