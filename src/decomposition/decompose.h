#ifndef BLAZE_DECOMPOSE
#define BLAZE_DECOMPOSE

#include "utils/utils.h"
#include "utils/seasonality_utils.h"
#include "utils/dot.h"

#include "decomposition/decompose_result.h"

enum DecomposeType{
  A,
  M,
  N
};

template <typename scalar_t>
std::vector<scalar_t> convolution_filter(
    const std::vector<scalar_t> &x,
    const std::vector<scalar_t> &filter)
{
  const size_t n_filter = filter.size();
  const size_t offset = n_filter/2;
  const size_t n = x.size();
  std::vector<scalar_t> result(n, std::nan("0.0"));
  for(size_t i = 0; i < (n-2*offset); i++) {
    result[i+offset] = dot<scalar_t>(x.data()+i, filter.data(), n_filter);
  }
  return result;
}

// TODO (?): untested
template <typename scalar_t>
std::vector<scalar_t> recursive_filter(
    const std::vector<scalar_t> &x,
    const std::vector<scalar_t> &filter)
{
  const size_t n_filter = filter.size(), n = x.size();
  std::vector<scalar_t> result(n + n_filter, 0.0), res(n, 0.0);;
  for(size_t i = 0; i < n; i++) {
    scalar_t sum = x[i];
    if (!std::isnan(sum)) { 
      for(size_t j = 0; j < n_filter; j++) {
        scalar_t tmp = result[n_filter + i - j - 1];
        if(!std::isnan(tmp)) sum += tmp * filter[j];
      }
      result[n_filter + i] = sum;
    }
    else {
      result[n_filter + i] = std::nan("0.0");
    }
  }
  for(size_t i=0; i < n; i++) res[i] = result[n_filter+i];
  return res;
}

enum FilterType{
  Conv,
  Rec
};

template <typename scalar_t,
          const FilterType type = FilterType::Conv>
std::vector<scalar_t> filter(
    const std::vector<scalar_t> &x,
    const std::vector<scalar_t> &filter) {
  if constexpr(type == FilterType::Conv) {
    return convolution_filter<scalar_t>(x, filter);
  }
  return recursive_filter(x, filter);
}

// The reimplementation of R stats::decompose
template <typename scalar_t,
          const DecomposeType type = DecomposeType::A>
DecompositionResult<scalar_t> decompose(
    const std::vector<scalar_t> &y,
    size_t season = 0) {
  
  static_assert(type != DecomposeType::N,
                "Decomposition of type N (for None) not permitted for decompose - please use type A (additive) or M(multiplicative)");
  const size_t n = y.size();
  // get main seasonality if not provided by the user
  if(season == 0) season = find_seasonalities(y, 1)[0];
  const bool seasonal = season % 2;
  std::vector<scalar_t> filter_coef(season + (!seasonal),
                                    1.0/static_cast<scalar_t>(season));
  // create equal weight filter coefficients 
  if (!seasonal) {
    filter_coef[0] = 0.5/static_cast<scalar_t>(season);
    filter_coef[filter_coef.size()-1] = 0.5/static_cast<scalar_t>(season);
  }
  std::vector<scalar_t> trend = filter(y, filter_coef);
  std::vector<scalar_t> detrended(n);
  if constexpr(type == DecomposeType::A) {
    for(size_t i = 0; i < n; i++) {
      detrended[i] = y[i] - trend[i];
    }
  } else if constexpr(type == DecomposeType::M) {
    for(size_t i = 0; i < n; i++) {
      detrended[i] = y[i]/trend[i];
    }
  }
  // print_vector(detrended);
  const size_t n_periods = n / season; // integer division
  std::vector<scalar_t> seasonality(n, 0.0);
  // compute seasonal averages 
  for(size_t i = 0; i < season; i++) {
    for(size_t j = 0; j < n_periods; j++) {
      const scalar_t temp = detrended[i + (j * season)];
      if(!std::isnan(temp)) {
        seasonality[i] += temp;
      }
    }
    seasonality[i] /= season;
  }
  for(size_t i = season; i <n; i++) {
    seasonality[i] = seasonality[i % season];
  }
  scalar_t seas_mean = mean(seasonality);
  if constexpr(type == DecomposeType::A) {
    for(auto &seas:seasonality) seas -= seas_mean;
  } else if constexpr(type == DecomposeType::M) {
    for(auto &seas:seasonality) seas /= seas_mean;
  }
  std::vector<scalar_t> remainder(n);
  
  for(size_t i = 0; i < n; i++) {
    if constexpr(type == DecomposeType::A) {
      remainder[i] = y[i] - seasonality[i] - trend[i];
    }
    else if constexpr(type == DecomposeType::M) {
      remainder[i] = y[i]/seasonality[i]/trend[i];
    }
  }
  return DecompositionResult<scalar_t>(trend, seasonality, remainder);
}

template <typename scalar_t,
          const DecomposeType type = DecomposeType::A>
std::vector<scalar_t> decompose_seasonality(
    const std::vector<scalar_t> &y,
    size_t season = 0) {
  
  static_assert(type != DecomposeType::N,
                "Decomposition of type N (for None) not permitted for decompose_seasonality - please use type A (additive) or M(multiplicative)");
  const size_t n = y.size();
  // get main seasonality if not provided by the user
  if(season == 0) season = find_seasonalities(y, 1)[0];
  const bool seasonal = season % 2;
  std::vector<scalar_t> filter_coef(season + (!seasonal),
                                    1.0/static_cast<scalar_t>(season));
  // create equal weight filter coefficients 
  if (!seasonal) {
    filter_coef[0] = 0.5/static_cast<scalar_t>(season);
    filter_coef[filter_coef.size()-1] = 0.5/static_cast<scalar_t>(season);
  }
  std::vector<scalar_t> trend = filter(y, filter_coef);
  const size_t n_periods = n / season; // integer division
  std::vector<scalar_t> seasonality(n, 0.0);
  // compute seasonal averages 
  for(size_t i = 0; i < season; i++) {
    scalar_t res = 0.0;
    for(size_t j = 0; j < n_periods; j++) {
      scalar_t temp = 0.0;
      if constexpr(type == DecomposeType::A) {
        temp = y[i + (j * season)] - trend[i + (j * season)];
      } else if constexpr(type == DecomposeType::M) {
        temp = y[i + (j * season)]/trend[i + (j * season)];
      }
      if(!std::isnan(temp)) res += temp;
    }
    seasonality[i] = res/season;
  }
  scalar_t seas_mean = 0.0;
  for(size_t i = 0; i < season; i++) seas_mean += seasonality[i];
  seas_mean /= season;
  for(size_t i = season; i <n; i++) seasonality[i] = seasonality[i % season];
  for(auto &seas:seasonality) {
    if constexpr(type == DecomposeType::A) seas -= seas_mean;
    else if constexpr(type == DecomposeType::M) seas /= seas_mean;
  }
  return seasonality;
}

#endif 
