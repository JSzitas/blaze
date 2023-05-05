#ifndef SEASONALITY_UTILS
#define SEASONALITY_UTILS

#include "utils/utils.h"
#include "auto_ar/auto_ar.h"

template <typename scalar_t> const size_t period(
    const std::vector<scalar_t>&y, const scalar_t threshold = 10.0) {
  
  const size_t n = y.size(); 
  const size_t max_p = min(n - 1L, std::floor(10 * std::log10(n)));
  
  auto res = AutoAR(y, 1, max_p,
                    std::vector<std::vector<scalar_t>>(0, std::vector<scalar_t>(0)),
                    true, false, true, AutoARMethod::AIC);
  
    order <- x$order
    n.freq <- 500
    freq <- seq.int(0, 0.5, length.out = n.freq)
    if (nser == 1) {
      coh <- phase <- NULL
      var.p <- as.vector(x$var.pred)
      spec <- if (order >= 1) {
        cs <- outer(freq, 1L:order, function(x, y) cos(2 * 
          pi * x * y)) %*% x$ar
        sn <- outer(freq, 1L:order, function(x, y) sin(2 * 
          pi * x * y)) %*% x$ar
        var.p/(xfreq * ((1 - cs)^2 + sn^2))
      }
      else rep.int(var.p/xfreq, length(freq))
    }
    return spec;
  
  
  spec <- stats::spec.ar(c(x))
  if(max(spec$spec)>threshold)
  {
    period <- round(1/spec$freq[which.max(spec$spec)])
    if(period==Inf) // Find next local maximum
    {
      j <- which(diff(spec$spec)>0)
      if(length(j)>0)
      {
        nextmax <- j[1] + which.max(spec$spec[j[1]:500])
        period <- round(1/spec$freq[nextmax])
      }
      else
        period <- 1
    }
  }
  else
    period <- 1
  return(period)
} 

#endif
