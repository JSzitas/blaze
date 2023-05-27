#ifndef SEASONALITY_UTILS
#define SEASONALITY_UTILS

#include "utils/utils.h"
#include "auto_ar/auto_ar.h"

template <typename scalar_t> size_t period(
    const std::vector<scalar_t>&y, const scalar_t threshold = 10,
    const size_t n_frequencies = 500) {
  const size_t n = y.size(), max_p = min(n - 1L, std::floor(10 * std::log10(n)));
  auto res = AutoAR(y, 1, max_p,
                    std::vector<std::vector<scalar_t>>(0, std::vector<scalar_t>(0)),
                    true, false, true, AutoARMethod::AICc);
  // do not forget to fit the model
  res.fit();
  auto freq = regular_sequence<scalar_t>(0.0, 0.5, n_frequencies);
  const size_t order = res.get_order();
  const scalar_t var_p = res.get_sigma2();
  auto coef = res.get_coef();

  static constexpr scalar_t two_pi = 6.283185;
  std::vector<scalar_t> spec(n_frequencies, var_p);
  if(order >= 1) {
    for(size_t i = 0; i < n_frequencies; i++) {
      scalar_t temp_cos = 0.0, temp_sin = 0.0; 
      for( size_t j = 0; j < order; j++ ) {
        scalar_t temp = two_pi * freq[i] * (j+1);
        temp_cos += std::cos(temp) * coef[j];
        temp_sin += std::sin(temp) * coef[j];
      }
      spec[i] = var_p/(std::pow(1 - temp_cos, 2) + std::pow(temp_sin,2));
    }
  }
  size_t period = 1, max_spec_index = 0, j = 0;
  scalar_t max_spec = 0.0;
  for(auto&item:spec) {
    if(item > max_spec) {
      max_spec = item;
      max_spec_index = j;
    }
    j++;
  }
  if(max_spec > threshold)
  {
    period = std::round(1.0/freq[max_spec_index]);
    // Find next local maximum
    if(std::isinf(period) || max_spec_index == 0) {
      std::vector<scalar_t> diff_spec_positive = diff(spec);
      std::vector<size_t> positive_indices;
      for(size_t j = 0; j < diff_spec_positive.size(); j++) {
        if(diff_spec_positive[j] > 0) {
          positive_indices.push_back(j);
        }
      }
      if(positive_indices.size() > 0)
      {
        // index in difference is '1 earlier' since our diff produces output of
        // size input-differences
        size_t nextmax = 1 + max_at(spec, positive_indices[0], spec.size());
        period = std::round(1/freq[nextmax]);
      }
    }
  }
  return period;
} 

template <typename scalar_t> std::vector<scalar_t> windowed_sum(
    const std::vector<scalar_t> &y,
    const size_t frequency = 1) {
  const size_t new_size = y.size()/frequency;
  std::vector<scalar_t> result(new_size);
  for(size_t j = 0; j < new_size; j++) {
    for(size_t i = 0; i < frequency; i++) {
      result[j] += y[(j*frequency) + i];
    }
  }
  return result;
}

template <typename scalar_t> std::vector<size_t> find_seasonalities(
  std::vector<scalar_t> y,
  const size_t max_iter = 5,
  const size_t upper_limit = 1500 ) {
  std::vector<size_t> periods;
  for(size_t j = 0; j < max_iter; j++) {
    size_t last_period = period(y);
    if(last_period <= 1 || (last_period > y.size())){
      break;
    }
    periods.push_back(last_period);
    y = windowed_sum(y, last_period);
  }
  std::vector<size_t> result = cummulative_product(periods);
  size_t last_valid_index = result.size();
  for(size_t j=0; j < result.size();j++) {
    if(result[j] > upper_limit) {
      last_valid_index = j; 
      break;
    }
  }
  result.resize(last_valid_index);
  return result;
}

#endif
