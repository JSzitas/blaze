#ifndef ARIMA_UTILS
#define ARIMA_UTILS

#include "cmath"

#include "utils/utils.h"

#include "arima/structures/arima_kind.h"
#include "arima/structures/structural_model.h"

#include "third_party/eigen.h"

// originally partrans
std::vector<double> parameter_transform(std::vector<double> &coef) {
  auto p = coef.size();
  int j, k;
  // we solve this by using a vector rather than an array
  std::vector<double> working_par(p);
  std::vector<double> new_par(p);
  /* Step one: map (-Inf, Inf) to (-1, 1) via tanh
   The parameters are now the pacf phi_{kk} */
  for (j = 0; j < p; j++) {
    new_par[j] = working_par[j] = tanh(coef[j]);
  }
  /* Step two: run the Durbin-Levinson recursions to find phi_{j.},
   j = 2, ..., p and phi_{p.} are the autoregression coefficients */
  for (j = 1; j < p; j++) {
    for (k = 0; k < j; k++) {
      // I believe allocating a is not necessary
      new_par[k] -= working_par[j] * working_par[j - k - 1];
    }
    for (k = 0; k < j; k++) {
      working_par[k] = new_par[k];
    }
  }
  // why not just return work directly
  return new_par;
}

// new partrans
template <class T> void parameter_transform(T &coef, int start, int end) {
  int j, k;
  // we solve this by using a vector rather than an array
  std::vector<double> new_par(coef.size());
  /* Step one: map (-Inf, Inf) to (-1, 1) via tanh
   The parameters are now the pacf phi_{kk} */
  for (j = start; j < end; j++) {
    coef[j] = new_par[j] = tanh(coef[j]);
  }
  /* Step two: run the Durbin-Levinson recursions to find phi_{j.},
   j = 2, ..., p and phi_{p.} are the autoregression coefficients */
  for (j = start; j < end; j++) {
    for (k = 0; k < j; k++) {
      // I believe allocating a is not necessary
      new_par[k] -= coef[j] * coef[j - k - 1];
    }
    for (k = 0; k < j; k++) {
      coef[k] = new_par[k];
    }
  }
}

// invert the parameter transformation
std::vector<double> inv_parameter_transform(std::vector<double> &phi) {
  auto p = phi.size();
  std::vector<double> new_pars(p);
  std::vector<double> work(p);
  for (int j = 0; j < p; j++) {
    // assigning work[j] = might be redundant - it gets reassigned anyways and
    // only holds intermediaries, but does not do anything for the result
    new_pars[j] = phi[j];
  }
  /* Run the Durbin-Levinson recursions backwards
   to find the PACF phi_{j.} from the autoregression coefficients */
  for (int j = p - 1; j > 0; j--) {
    // this allocation is redundant, we can just take reference to new_pars[j]
    // a = new_pars[j];
    for (int k = 0; k < j; k++) {
      work[k] = (new_pars[k] + new_pars[j] * new_pars[j - k - 1]) /
                (1 - new_pars[j] * new_pars[j]);
    }
    for (int k = 0; k < j; k++) {
      new_pars[k] = work[k];
    }
  }
  // revert the tanh transform from earlier
  for (int j = 0; j < p; j++) {
    new_pars[j] = atanh(new_pars[j]);
  }
  return new_pars;
}

void inv_parameter_transform2(std::vector<double> &phi) {

  auto p = phi.size();
  std::vector<double> temp(p * 2);
  // std::vector<double> new_pars(p);
  // std::vector<double> work(p);
  for (int j = 0; j < p; j++) {
    // assigning work[j] = might be redundant - it gets reassigned anyways and
    // only holds intermediaries, but does not do anything for the result
    temp[j] = phi[j];
  }
  /* Run the Durbin-Levinson recursions backwards
   to find the PACF phi_{j.} from the autoregression coefficients */
  for (int j = p - 1; j > 0; j--) {
    // this allocation is redundant, we can just take reference to new_pars[j]
    // a = new_pars[j];
    for (int k = 0; k < j; k++) {
      temp[p + k] =
          (temp[k] + temp[j] * temp[j - k - 1]) / (1 - temp[j] * temp[j]);
    }
    for (int k = 0; k < j; k++) {
      temp[k] = temp[p + k];
    }
  }
  // revert the tanh transform from earlier
  for (int j = 0; j < p; j++) {
    phi[j] = atanh(temp[j]);
  }
}

// this just directly modifies coef
template <const bool seasonal, const bool transform>
void arima_transform_parameters(Eigen::VectorXd &coef, const arima_kind &arma,
                                std::vector<double> &phi,
                                std::vector<double> &theta) {
  // the coefficients are all 'packed in' inside coef - so we have
  // different types of coefficients. this tells us basically how many
  // of each type there are
  const int mp = arma.p(), msp = arma.P(), mq = arma.q();
  if constexpr (transform) {
    if (mp > 0) {
      parameter_transform(coef, 0, mp);
    }
    if (msp > 0) {
      parameter_transform(coef, mp + mq, msp);
    }
  }
  if constexpr (seasonal) {
    const int msq = arma.Q(), ns = arma.period();
    const int p = mp + ns * msp;
    const int q = mq + ns * msq;
    int i, j;
    /* expand out seasonal ARMA models
     * note that the offsetting here is crucial - the original indexing was
     * into two data structures of the same combined size as our coef */
    for (i = 0; i < mp; i++) {
      phi[i] = coef[i];
    }
    for (i = 0; i < mq; i++) {
      theta[i] = coef[i + mp];
    }
    for (i = mp; i < p; i++) {
      phi[i] = 0.0;
    }
    for (i = mq; i < q; i++) {
      theta[i] = 0.0;
    }
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
    for (i = 0; i < p; i++) {
      coef[i] = phi[i];
    }
    for (i = p; i < p + q; i++) {
      coef[i] = theta[i - p];
    }
  }
}

// this just directly modifies coef
template <const bool seasonal, const bool transform, class T>
void arima_transform_parameters(T &coef, const arima_kind &arma) {
  const int mp = arma.p(), msp = arma.P(), mq = arma.q();
  if constexpr (transform) {
    if (mp > 0) {
      parameter_transform(coef, 0, mp);
    }
    if (msp > 0) {
      parameter_transform(coef, mp + mq, msp);
    }
  }
  if constexpr (seasonal) {
    const int msq = arma.Q(), ns = arma.period();
    const int p = mp + ns * msp;
    const int q = mq + ns * msq;

    std::vector<double> phi(p);
    std::vector<double> theta(q);
    int i, j, v;
    for (i = 0; i < mp; i++) {
      phi[i] = coef[i];
    }
    for (i = 0; i < mq; i++) {
      theta[i] = coef[i + mp];
    }
    for (i = mp; i < p; i++) {
      phi[i] = 0.0;
    }
    for (i = mq; i < q; i++) {
      theta[i] = 0.0;
    }
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
    for (i = 0; i < p; i++) {
      coef[i] = phi[i];
    }
    for (i = p; i < p + q; i++) {
      coef[i] = theta[i - p];
    }
  }
}

void arima_transform_parameters(structural_model<double> model, arima_kind kind,
                                bool transform = true) {
  // the coefficients are all 'packed in' inside coef - so we have
  // different types of coefficients. this tells us basically how many
  // of each type there are
  int mp = kind.p(), mq = kind.q(), msp = kind.P(), msq = kind.Q(),
      ns = kind.period();
  int p = mp + ns * msp;
  int q = mq + ns * msq;

  int n = mp + mq + msp + msq;
  std::vector<double> params(n);
  int i, j, v;
  for (i = 0; i < mp; i++) {
    params[i] = model.phi[i];
  }
  for (i = mp; i < mq; i++) {
    params[i] = model.theta[i];
  }
  for (i = mq; i < msp; i++) {
    params[i] = model.phi[i];
  }
  for (i = msp; i < msq; i++) {
    params[i] = model.theta[i];
  }

  if (transform) {
    std::vector<double> temp(mp);
    if (mp > 0) {
      for (i = 0; i < mp; i++) {
        temp[i] = params[i];
      }
      temp = parameter_transform(temp);
      for (i = 0; i < temp.size(); i++) {
        params[i] = temp[i];
      }
    }
    v = mp + mq;
    if (msp > 0) {
      // this is a transformation over a view
      // ie parameters v and higher
      // create a copy
      temp.resize(msp);
      // move values to a temporary
      for (i = v; i < msp; i++) {
        temp[i - v] = params[i];
      }
      // overwrite
      temp = parameter_transform(temp);
      // write back to parameters
      for (i = v; i < msp; i++) {
        params[i] = temp[i - v];
      }
    }
  }
  if (ns > 0) {
    /* expand out seasonal ARMA models */
    for (i = 0; i < mp; i++) {
      model.phi[i] = params[i];
    }
    for (i = 0; i < mq; i++) {
      model.theta[i] = params[i + mp];
    }
    for (i = mp; i < p; i++) {
      model.phi[i] = 0.0;
    }
    for (i = mq; i < q; i++) {
      model.theta[i] = 0.0;
    }
    for (j = 0; j < msp; j++) {
      model.phi[(j + 1) * ns - 1] += params[j + mp + mq];
      for (i = 0; i < mp; i++) {
        model.phi[(j + 1) * ns + i] -= params[i] * params[j + mp + mq];
      }
    }
    for (j = 0; j < msq; j++) {
      model.theta[(j + 1) * ns - 1] += params[j + mp + mq + msp];
      for (i = 0; i < mq; i++) {
        model.theta[(j + 1) * ns + i] +=
            params[i + mp] * params[j + mp + mq + msp];
      }
    }
  } else {
    for (i = 0; i < mp; i++) {
      model.phi[i] = params[i];
    }
    for (i = 0; i < mq; i++) {
      model.theta[i] = params[i + mp];
    }
  }
}

std::vector<double> arima_inverse_transform_parameters(std::vector<double> coef,
                                                       std::vector<int> &arma) {

  // arma contains the arma structure - 3 leading numbers,
  // number of p parameters (arma[0])
  // number of q parameters (arma[1])
  // number of seasonal p parameters (arma[2])
  // for a function this small, it makes little sense to copy over
  std::vector<double> temp(arma[0]);
  if (arma[0] > 0) {
    for (int i = 0; i < arma[0]; i++) {
      temp[i] = coef[i];
    }
    temp = parameter_transform(temp);
    for (int i = 0; i < arma[0]; i++) {
      coef[i] = temp[i];
    }
  }

  if (arma[2] > 0) {
    temp.resize(arma[2]);
    for (int i = (arma[0] + arma[1]); i < coef.size(); i++) {
      temp[i] = coef[i];
    }
    temp = parameter_transform(temp);
    for (int i = (arma[0] + arma[1]); i < coef.size(); i++) {
      coef[i] = temp[i];
    }
  }
  return coef;
}

std::vector<double> arima_grad_transform(std::vector<double> &coef,
                                         std::vector<int> &arma,
                                         double eps = 1e-3) {
  int mp = arma[0], mq = arma[1], msp = arma[2], n = coef.size();
  std::vector<double> w1_temp(mp), w2_temp(mp), w3_temp(mp);

  std::vector<double> result(n * n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      result[i + j * n] = (i == j);
    }
  }
  if (mp > 0) {
    for (int i = 0; i < mp; i++) {
      w1_temp[i] = coef[i];
    }
    w2_temp = parameter_transform(w1_temp);
    for (int i = 0; i < mp; i++) {
      w1_temp[i] += eps;
      w3_temp = parameter_transform(w1_temp);
      for (int j = 0; j < mp; j++) {
        result[i + j * n] = (w3_temp[j] - w2_temp[j]) / eps;
      }
      w1_temp[i] -= eps;
    }
  }
  if (msp > 0) {
    int v = mp + mq;
    w1_temp.resize(msp);
    w2_temp.resize(msp);
    w3_temp.resize(msp);
    for (int i = 0; i < msp; i++) {
      w1_temp[i] = coef[i + v];
    }
    w2_temp = parameter_transform(w1_temp);
    for (int i = 0; i < msp; i++) {
      w1_temp[i] += eps;
      w1_temp = parameter_transform(w3_temp);
      for (int j = 0; j < msp; j++) {
        result[i + v + (j + v) * n] = (w3_temp[j] - w2_temp[j]) / eps;
      }
      w1_temp[i] -= eps;
    }
  }
  // result is basically a packed matrix
  return result;
}

#endif
