#ifndef BLAZE_ACF
#define BLAZE_ACF

#include "utils/utils.h"

template <typename scalar_t,
          const bool correlation = true,
          const bool demean = true>
std::vector<scalar_t> acf(const std::vector<scalar_t> &y, 
                          size_t lag_max = 0) {
  const size_t n = y.size(),
    max_lag = lag_max == 0 ? std::floor(10 * (std::log10(n))) : lag_max;
  std::vector<scalar_t> result(max_lag + 1, 0.0), x = y; 
  if constexpr(demean) {
    const scalar_t y_mean = mean(y);
    for(auto &val:x) {
      val -= y_mean;
    }
  }
  for(size_t lag = 0; lag <= max_lag; lag++) {
    double sum = 0.0; 
    for(size_t i = lag; i < n; i++) {
      sum += x[i] * x[i-lag];
    }
    result[lag] = sum/n;
  }
  if constexpr(correlation) {
    if(n == 1) result[0] = 1.0;
    else {
      const scalar_t scale = result[0];
      for(size_t lag = 0; lag <= max_lag; lag++) {
        scalar_t scaled = result[lag]/scale;
        result[lag] = (scaled > 1.) ? 1. : ((scaled < -1.) ? -1. : scaled);
      }
    }
  }
  return result;
}

#endif
