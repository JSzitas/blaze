#ifndef ARIMA_HUBER_SOLVER
#define ARIMA_HUBER_SOLVER

#include "utils/utils.h"
#include "utils/xreg.h"

#include "arima/structures/structures.h"

#include "arima/solvers/css_likelihood.h"
#include "arima/solvers/state_space.h"

#include "arima/utils/transforms.h"

// include optimizer library
#include "third_party/eigen.h"
#include "third_party/optim.h"

template <typename C, typename scalar_t,
          const bool arima_only=false,
          const bool pseudo=true> scalar_t arima_css_huber(
    const C & y, const C & pars, const arima_kind kind,
    const scalar_t delta,
    const size_t n_cond, std::vector<scalar_t> resid) {
  const size_t n = y.size(), p = kind.p() + kind.period() * kind.P(),
        q = kind.q() + kind.period() * kind.Q();
  // prepare the residuals - possibly move this out and never allocate here?
  int ma_offset, nu = n-n_cond;
  scalar_t ssq = 0.0, tmp = 0.0, delta_sqr = delta * delta;
  for (size_t l = n_cond; l < n; l++) {
    ma_offset = min(l - n_cond, q);
    tmp = y[l];
    for (size_t j = 0; j < p; j++) {
      tmp -= pars[j] * y[l - j - 1];
    }
    // to offset that this is all in one vector, we need to
    // start at p and go to p + q
    for (size_t j = 0; j < ma_offset; j++) {
      tmp -= pars[p + j] * resid[l - j - 1];
    }
    resid[l] = tmp;
    if constexpr(pseudo) {
      if (!isnan(tmp)) {
        const scalar_t abs_val = std::abs(tmp);
        ssq += abs_val <= delta ? 1/2 * tmp * tmp : delta * abs_val;
      } else {
        nu--;
      }
    }
    if constexpr(!pseudo) {
      if (!isnan(tmp)) {
        ssq += delta_sqr * (std::sqrt(1 + std::pow(tmp/delta, 2))-1);
      } else {
        nu--;
      }
    }
  }
  return ssq/nu;
}

#endif
