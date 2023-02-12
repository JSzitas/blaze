#ifndef ARIMA_UTILS
#define ARIMA_UTILS

#include "cmath"

#include "utils/utils.h"
#include "utils/poly.h"

#include "arima/structures/arima_kind.h"
#include "arima/structures/structural_model.h"

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

template <typename T> void invert_ma(
    T & coef, const size_t start, const size_t end,
    const double tol = 0.0000001) {
  auto q = end - start;
  // analogous to ar_check >> find latest/largest order nonzero coef
  size_t q0 = 0;
  for (size_t i = 0; i < q; i++) {
    if (abs(coef[i]) > tol) {
      q0 = i;
    }
  }
  // if all coefs zero
  if (!q0) {
    return;
  }
  // otherwise find roots
  std::vector<double> temp(q0 + 1);
  temp[0] = 1;
  for (size_t i = 1; i < temp.size(); i++) {
    temp[i] = coef[start + i - 1];
  }
  auto roots = polyroot(temp);
  // find roots that are smaller than 1 - these are not inverse roots,
  // so they should be smaller, rather than larger than 1
  std::vector<bool> indices(roots.size());
  int any = 0;
  bool temp_id = false;
  for (size_t i = 0; i < roots.size(); i++) {
    temp_id = abs(roots[i]) < 1;
    indices[i] = temp_id;
    // this simplifies the all(!ind) part
    any += temp_id;
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
    temp[0] = 1 / temp[0];
    for (size_t i = 1; i < temp.size(); i++) {
      temp[i] = 0;
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
    coef[start + i - 1] = std::real(result[i]);
  }
  // if necessary pad with 0s
  for (size_t i = result.size(); i < q; i++) {
    coef[start + i - 1] = 0;
  }
}

std::vector<std::complex<double>> invert_ma_coef_from_roots(
    std::vector<std::complex<double>> roots ) {
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

template <class T> void parameter_transform(
    T &coef, const size_t start, const size_t end) {
  size_t j, k;
  // we solve this by using a vector rather than an array
  const size_t p = end-start;
  std::vector<double> new_par(p);
  /* Step one: map (-Inf, Inf) to (-1, 1) via tanh
   The parameters are now the pacf phi_{kk} */
  for (j = 0; j < p; j++) {
    coef[start+j] = new_par[j] = tanh(coef[start+j]);
  }
  /* Step two: run the Durbin-Levinson recursions to find phi_{j.},
   j = 2, ..., p and phi_{p.} are the autoregression coefficients */
  for (j = 0; j < p; j++) {
    for (k = 0; k < j; k++) {
      // I believe allocating a is not necessary
      new_par[k] -= coef[start + j] * coef[start + j - k - 1];
    }
    for (k = 0; k < j; k++) coef[start+k] = new_par[k];
  }
}

// invert the parameter transformation
template <typename T> void inv_parameter_transform(
    T & coef, const size_t start, const size_t end) {
  const auto p = end - start;
  std::vector<double> new_pars(p);
  for (size_t j = 0; j < p; j++) {
    // assigning work[j] = might be redundant - it gets reassigned anyways and
    // only holds intermediaries, but does not do anything for the result
    new_pars[j] = coef[start + j];
  }
  /* Run the Durbin-Levinson recursions backwards
   to find the PACF phi_{j.} from the autoregression coefficients */
  for (size_t j = p-1; j > 0; j--) {
    // this allocation is redundant, we can just take reference to new_pars[j]
    a = new_pars[j];
    for (size_t k = 0; k < j; k++) {
      coef[start + k] = (new_pars[k] + a * new_pars[j - k - 1]) / (1 - a * a);
    }
    for (size_t k = 0; k < j; k++) new_pars[k] = coef[k+start];
  }
  // revert the tanh transform from earlier
  for (size_t j = 0; j < p; j++) coef[j+start] = atanh(new_pars[j]);
}

// this just directly modifies coef
template <typename T, const bool seasonal, const bool transform>
void arima_transform_parameters(
    T &coef,
    const arima_kind &arma,
    std::vector<double> &phi,
    std::vector<double> &theta) {
  // the coefficients are all 'packed in' inside coef - so we have
  // different types of coefficients. this tells us basically how many
  // of each type there are
  const size_t mp = arma.p(), msp = arma.P(), mq = arma.q();
  if constexpr (transform) {
    if (mp > 0) parameter_transform(coef, 0, mp);
    if (msp > 0) parameter_transform(coef, mp + mq, mp + mq + msp);
  }
  if constexpr (seasonal) {
    const size_t msq = arma.Q(), ns = arma.period();
    const size_t p = mp + ns * msp;
    const size_t q = mq + ns * msq;
    size_t i, j;
    /* expand out seasonal ARMA models
     * note that the offsetting here is crucial - the original indexing was
     * into two data structures of the same combined size as our coef */
    for (i = 0; i < mp; i++) phi[i] = coef[i];
    for (i = 0; i < mq; i++) theta[i] = coef[i + mp];
    for (i = mp; i < p; i++) phi[i] = 0.0;
    for (i = mq; i < q; i++) theta[i] = 0.0;
    for (j = 0; j < msp; j++) {
      phi[(j + 1) * ns - 1] += coef[j + mp + mq];
      for (i = 0; i < mp; i++) {
        phi[(j + 1) * ns + i] -= coef[i] * coef[j + mp + mq];
      }
    }
    for (j = 0; j < msq; j++) {
      theta[(j + 1) * ns - 1] += coef[j + mp + mq + msp];
      for (i = 0; i < mq; i++) {
        theta[(j + 1) * ns + i] += coef[i + mp] * coef[j + mp + mq + msp];
      }
    }
    // the output is written back to coef
    for (i = 0; i < p; i++) coef[i] = phi[i];
    for (i = p; i < p + q; i++) coef[i] = theta[i - p];
  }
}

template <typename T> void arima_inverse_transform_parameters(
    T & coef,
    const arima_kind & kind) {
  // transform AR coefficients
  if(kind.p()) inv_parameter_transform(coef, 0, kind.p());
  // transform SAR coefficients
  if(kind.P()) inv_parameter_transform(coef,
     kind.p() + kind.q(),
     kind.p() + kind.q() + kind.P());
  // transform MA by inverting
  if(kind.q()) invert_ma(coef, kind.p(), kind.p() + kind.q());
  // // and SMA
  if(kind.Q()) invert_ma(coef, kind.p() + kind.q() + kind.P(),
                         kind.p() + kind.q()  + kind.P() + kind.Q());
}

// std::vector<double> arima_grad_transform(std::vector<double> &coef,
//                                          std::vector<int> &arma,
//                                          double eps = 1e-3) {
//   int mp = arma[0], mq = arma[1], msp = arma[2], n = coef.size();
//   std::vector<double> w1_temp(mp), w2_temp(mp), w3_temp(mp);
//
//   std::vector<double> result(n * n);
//   for (int i = 0; i < n; i++) {
//     for (int j = 0; j < n; j++) {
//       result[i + j * n] = (i == j);
//     }
//   }
//   if (mp > 0) {
//     for (int i = 0; i < mp; i++) {
//       w1_temp[i] = coef[i];
//     }
//     w2_temp = parameter_transform(w1_temp);
//     for (int i = 0; i < mp; i++) {
//       w1_temp[i] += eps;
//       w3_temp = parameter_transform(w1_temp);
//       for (int j = 0; j < mp; j++) {
//         result[i + j * n] = (w3_temp[j] - w2_temp[j]) / eps;
//       }
//       w1_temp[i] -= eps;
//     }
//   }
//   if (msp > 0) {
//     int v = mp + mq;
//     w1_temp.resize(msp);
//     w2_temp.resize(msp);
//     w3_temp.resize(msp);
//     for (int i = 0; i < msp; i++) {
//       w1_temp[i] = coef[i + v];
//     }
//     w2_temp = parameter_transform(w1_temp);
//     for (int i = 0; i < msp; i++) {
//       w1_temp[i] += eps;
//       w1_temp = parameter_transform(w3_temp);
//       for (int j = 0; j < msp; j++) {
//         result[i + v + (j + v) * n] = (w3_temp[j] - w2_temp[j]) / eps;
//       }
//       w1_temp[i] -= eps;
//     }
//   }
//   // result is basically a packed matrix
//   return result;
// }

#endif
