#ifndef BLAZE_ACF
#define BLAZE_ACF

template <typename scalar_t,
          const bool correlation = true>
std::vector<scalar_t> acf(const std::vector<scalar_t> &y, 
                          size_t lag_max = 0) {
  const size_t n = y.size(),
    max_lag = lag_max == 0 ? std::floor(10 * (std::log10(n))) : lag_max;
  std::vector<scalar_t> result(lag_max, 0.0);
  for(size_t lag = 0; lag <= max_lag; lag++) {
    double sum = 0.0; size_t nu = 0;
    for(size_t i = 0; i < n-lag; i++) {
      nu++;
      sum += y[i + lag] * y[i];
    }
    result[lag] = (nu > 0) ? sum/(nu + lag) : std::nan("0.0");
  }
  if constexpr(correlation) {
    if(n == 1) result[0] = 1.0;
    else {
      const scalar_t scale = std::sqrt(result[0]);
      for(size_t lag = 0; lag <= max_lag; lag++) {
        scalar_t scaled = result[lag]/scale;
        result[lag] = (scaled > 1.) ? 1. : ((scaled < -1.) ? -1. : scaled);
      }
    }
  }
  return result;
}

#endif
