#ifndef ARIMA_CHECKS
#define ARIMA_CHECKS

#include "utils/poly.h"
#include "arima/structures/arima_kind.h"

template <typename scalar_t=float> bool ar_check(
    const std::vector<scalar_t> &ar_coef,
    const scalar_t tol = 0.0000001) {
  // check if all inverse coefficients are zero - in that case we return
  size_t p = 0;
  for (size_t i = 0; i < ar_coef.size(); i++) {
    if (abs(-ar_coef[i]) > tol) {
      p = i;
    }
  }
  if (!p) {
    return true;
  }
  // the polyroot part, e.g.
  // all(Mod(polyroot(c(1, -ar[1L:p]))) > 1)
  // goes here instead of the false
  // then the abs(/inverse/roots) must be all greater than > 1 for
  // this to lie within the unit circle
  std::vector<scalar_t> coef(p + 1);
  coef[0] = 1;
  for (size_t i = 1; i < coef.size(); i++) {
    coef[i] = -ar_coef[i - 1];
  }
  auto res = polyroot(coef);
  // check coefficients
  for (size_t i = 0; i < res.size(); i++) {
    if (abs(res[i]) <= 1) {
      // the all(Mod(x) > 1) part
      // if any modulus from inverse coefficients is smaller than 1,
      // that implies that the coefficient lies outside the unit circle
      // and the process is not stationary
      return false;
    }
  }
  // otherwise all coefficients are ok >> return true
  return true;
}

template <typename C,
          typename scalar_t=float> std::array<bool,2> check_all_ar(
    const C & coef,
    const arima_kind & kind,
    const scalar_t tol = 0.0000001) {
  // preallocate result
  std::array<bool,2> result;
  const size_t p = kind.p(), P = kind.P();
  // first check for the non-seasonal coef
  std::vector<scalar_t> temp(p);
  for(size_t i = 0; i < p; i++) {
    temp[i] = coef[i];
  }
  result[0] = ar_check(temp);
  // if P is 0, this passes
  if( P == 0 ) {
    result[1] = true;
    return result;
  }
  // otherwise resize temp, fill with seasonal coef, and repeat
  temp.resize(P);
  const size_t q = kind.q();
  for(size_t i = 0; i < P; i++) {
    temp[i] = coef[p + q + i];
  }
  result[1] = ar_check(temp);
  return result;
}

#endif
