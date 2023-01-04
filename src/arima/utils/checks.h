#ifndef CHECKS
#define CHECKS

#include "utils/poly.h"

bool ar_check(
    const std::vector<double> &ar_coef,
    const double tol = 0.0000001) {
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
  std::vector<double> coef(p + 1);
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

std::complex<double> operator/(
    const double left,
    const std::complex<double> right) {
  /* perform /(double,complex), or div(double,complex), ie
   * div(complex, complex) where the first element has complex part == 0
   * complex division is given by formula:
   * (ac+bd)/(c2 + d2) + (bc-ad)/(c2 + d2)i.
   * which for first operand double takes b = 0
   */
  const double &a = left;
  const double &c = std::real(right);
  const double &d = std::imag(right);
  /* Thus the formula simplifies to:
   * (ac)/(c^2 + d^2) + (-ad)/(c^2 + d^2)i
   */
  std::complex<double> res(
      // real part
      (a * c) / (pow(c, 2) + pow(d, 2)),
      // complex part
      (-a * d) / (pow(c, 2) + pow(d, 2)));
  return res;
}

std::vector<double> ma_invert(
    std::vector<double> &ma_coef,
    const double tol = 0.0000001) {
  auto q = ma_coef.size();
  // analogous to ar_check >> find latest/largest order nonzero coef
  size_t q0 = 0;
  for (size_t i = 0; i < ma_coef.size(); i++) {
    if (abs(ma_coef[i]) > tol) {
      q0 = i;
    }
  }
  // if all coefs zero
  if (!q0) {
    return ma_coef;
  }
  // otherwise find roots
  std::vector<double> coef(q0 + 1);
  coef[0] = 1;
  for (size_t i = 1; i < coef.size(); i++) {
    coef[i] = ma_coef[i - 1];
  }
  auto roots = polyroot(coef);
  // find roots that are smaller than 1 - these are not inverse roots,
  // so they should be smaller, rather than larger than 1
  std::vector<bool> indices(roots.size());
  int any = 0;
  bool temp = false;
  for (size_t i = 0; i < roots.size(); i++) {
    temp = abs(roots[i]) < 1;
    indices[i] = temp;
    // this simplifies the all(!ind) part
    any += temp;
  }
  // if all are > 1, we do not need to change them
  // as this implies they are valid inverse roots already - and there is no
  // need to invert
  if (!any) {
    return ma_coef;
  }
  // otherwise if there is only 1 valid coefficient
  if (q0 == 1) {
    coef.resize(q);
    // only invert the very first coefficient and otherwise set other
    // coefficients to zero
    coef[0] = 1 / ma_coef[0];
    for (size_t i = 1; i < coef.size(); i++) {
      coef[i] = 0;
    }
    return coef;
  }
  // otherwise invert the roots which need to be inverted
  for (size_t i = 0; i < roots.size(); i++) {
    // I made this work for complex numbers :)
    if (indices[i]) {
      roots[i] = 1 / roots[i];
    }
  }
  std::vector<std::complex<double>> result(roots.size() + 1);
  // set first new coef == 1
  result[0] = 1;
  std::vector<std::complex<double>> inv_temp(roots.size() + 1);
  for (size_t i = 0; i < roots.size(); i++) {
    // take the root
    inv_temp = result;
    for (size_t k = 1; k <= i + 1; k++) {
      result[k] = inv_temp[k] - (inv_temp[k - 1] / roots[i]);
    }
  }
  coef.resize(ma_coef.size());
  for (size_t i = 1; i < result.size(); i++) {
    coef[i] = std::real(result[i]);
  }
  // if necessary pad with 0s
  for (size_t i = result.size(); i < q; i++) {
    coef[i] = 0;
  }
  return coef;
}

void ma_invert2(
    std::vector<double> &ma_coef,
    const double tol = 0.0000001) {
  auto q = ma_coef.size();
  // analogous to ar_check >> find latest/largest order nonzero coef
  size_t q0 = 0;
  for (size_t i = 0; i < ma_coef.size(); i++) {
    if (abs(ma_coef[i]) > tol) {
      q0 = i;
    }
  }
  // if all coefs zero
  if (!q0) {
    return;
  }
  // otherwise find roots
  std::vector<double> coef(q0 + 1);
  coef[0] = 1;
  for (size_t i = 1; i < coef.size(); i++) {
    coef[i] = ma_coef[i - 1];
  }
  auto roots = polyroot(coef);
  // find roots that are smaller than 1 - these are not inverse roots,
  // so they should be smaller, rather than larger than 1
  std::vector<bool> indices(roots.size());
  int any = 0;
  bool temp = false;
  for (size_t i = 0; i < roots.size(); i++) {
    temp = abs(roots[i]) < 1;
    indices[i] = temp;
    // this simplifies the all(!ind) part
    any += temp;
  }
  // if all are > 1, we do not need to change them
  // as this implies they are valid inverse roots already - and there is no
  // need to invert
  if (!any) {
    return;
  }
  // otherwise if there is only 1 valid coefficient
  if (q0 == 1) {
    // only invert the very first coefficient and otherwise set other
    // coefficients to zero
    ma_coef[0] = 1 / ma_coef[0];
    for (size_t i = 1; i < ma_coef.size(); i++) {
      ma_coef[i] = 0;
    }
    return;
  }
  // otherwise invert the roots which need to be inverted
  for (size_t i = 0; i < roots.size(); i++) {
    // I made this work for complex numbers :)
    if (indices[i]) {
      roots[i] = 1 / roots[i];
    }
  }
  std::vector<std::complex<double>> result(roots.size() + 1);
  // set first new coef == 1
  result[0] = 1;
  std::vector<std::complex<double>> inv_temp(roots.size() + 1);
  for (size_t i = 0; i < roots.size(); i++) {
    // take the root
    inv_temp = result;
    for (size_t k = 1; k <= i + 1; k++) {
      result[k] = inv_temp[k] - (inv_temp[k - 1] / roots[i]);
    }
  }
  for (size_t i = 1; i < result.size(); i++) {
    ma_coef[i] = std::real(result[i]);
  }
  // if necessary pad with 0s
  for (size_t i = result.size(); i < q; i++) {
    ma_coef[i] = 0;
  }
}

std::vector<std::complex<double>>
invert_ma_coef_from_roots( std::vector<std::complex<double>> roots ) {
  std::vector<std::complex<double>> result(roots.size() + 1);
  // set first new coef == 1
  result[0] = 1;
  std::vector<std::complex<double>> temp(roots.size() + 1);
  for (size_t i = 0; i < roots.size(); i++) {
    // take the root
    temp = result;
    for (size_t k = 1; k <= i + 1; k++) {
      result[k] = temp[k] - (temp[k - 1] / roots[i]);
    }
  }
  return result;
}

#endif
