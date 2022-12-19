#ifndef ARIMA_UTILS_2
#define ARIMA_UTILS_2

#include "cmath"
#include "utils/utils.h"
#include "utils/poly.h"
#include "structural_model.h"
#include "third_party/eigen.h"

// defines the arima structure
struct arima_kind{
  arima_kind(){};
  arima_kind( int p, int d, int q, int P, int D, int Q, int s_period ){
    this->arma_p = p;
    this->diff_d = d;
    this->arma_q = q;
    this->sarma_P = P;
    this->seas_diff_D = D;
    this->sarma_Q = Q;
    this->s_period = s_period;
  }
  arima_kind( std::vector<int> arima_def ) {
    this->arma_p = arima_def[0];
    this->diff_d = arima_def[1];
    this->arma_q = arima_def[2];
    this->sarma_P = arima_def[3];
    this->seas_diff_D = arima_def[4];
    this->sarma_Q = arima_def[5];
    this->s_period = arima_def[6];
  }
  const int p() const {
    return this->arma_p;
  }
  const int d() const {
    return this->diff_d;
  }
  const int q() const {
    return this->arma_q;
  }
  const int P() const {
    return this->sarma_P;
  }
  const int D() const {
    return this->seas_diff_D;
  }
  const int Q() const {
    return this->sarma_Q;
  }
  const int period() const {
    return this->s_period;
  }
private:
  int arma_p, diff_d, arma_q, sarma_P, seas_diff_D, sarma_Q, s_period;
};

// originally partrans
std::vector<double> parameter_transform( std::vector<double> & coef ) {
  auto p = coef.size();
  int j, k;
  // we solve this by using a vector rather than an array
  std::vector<double> working_par(p);
  std::vector<double> new_par(p);
  /* Step one: map (-Inf, Inf) to (-1, 1) via tanh
   The parameters are now the pacf phi_{kk} */
  for(j = 0; j < p; j++) {
    new_par[j] = working_par[j] = tanh(coef[j]);
  }
  /* Step two: run the Durbin-Levinson recursions to find phi_{j.},
   j = 2, ..., p and phi_{p.} are the autoregression coefficients */
  for(j = 1; j < p; j++) {
    for(k = 0; k < j; k++) {
      // I believe allocating a is not necessary
      new_par[k] -= working_par[j] * working_par[j - k - 1];
    }
    for(k = 0; k < j; k++) {
      working_par[k] = new_par[k];
    }
  }
  // why not just return work directly
  return new_par;
}
// new partrans
void parameter_transform2( std::vector<double> & coef, int start, int end ) {
  auto p = coef.size();
  int j, k;
  // we solve this by using a vector rather than an array
  // std::vector<double> working_par(p);
  std::vector<double> new_par(p);
  /* Step one: map (-Inf, Inf) to (-1, 1) via tanh
   The parameters are now the pacf phi_{kk} */
  for(j = start; j < end; j++) {
    coef[j] = new_par[j] = tanh(coef[j]);

    // new_par[j] = working_par[j] = tanh(coef[j]);
  }
  /* Step two: run the Durbin-Levinson recursions to find phi_{j.},
   j = 2, ..., p and phi_{p.} are the autoregression coefficients */
  for(j = start; j < end; j++) {
    for(k = 0; k < j; k++) {
      // I believe allocating a is not necessary
      new_par[k] -= coef[j] * coef[j - k - 1];
    }
    for(k = 0; k < j; k++) {
      coef[k] = new_par[k];
    }
  }
  // why not just return work directly
  // return new_par;
}



// this just directly modifies coef
void arima_transform_parameters( std::vector<double> &coef,
                                 const arima_kind &arma,
                                 bool transform = true)
{
  // the coefficients are all 'packed in' inside coef - so we have
  // different types of coefficients. this tells us basically how many
  // of each type there are
  int mp = arma.p(), mq = arma.q(), msp = arma.P(), msq = arma.Q(), ns = arma.period();
  int p = mp + ns * msp;
  int q = mq + ns * msq;

  std::vector<double> phi(p);
  std::vector<double> theta(q);
  int n = mp + mq + msp + msq;
  std::vector<double> params(n);
  int i, j, v;
  for (i = 0; i < coef.size(); i++) {
    params[i] = coef[i];
  }

  if (transform) {
    std::vector<double> temp(mp);
    if (mp > 0) {
      for(i = 0; i < mp; i++) {
        temp[i] = params[i];
      }
      temp = parameter_transform(temp);
      for(i = 0; i < temp.size(); i++) {
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
      for( i = v; i < msp; i++ ) {
        temp[i-v] = coef[i];
      }
      // overwrite
      temp = parameter_transform(temp);
      // write back to parameters
      for( i = v; i < msp; i++ ) {
        params[i] = temp[i-v];
      }
    }
  }
  if (ns > 0) {
    /* expand out seasonal ARMA models */
    for (i = 0; i < mp; i++) {
      phi[i] = params[i];
    }
    for (i = 0; i < mq; i++) {
      theta[i] = params[i + mp];
    }
    for (i = mp; i < p; i++) {
      phi[i] = 0.0;
    }
    for (i = mq; i < q; i++) {
      theta[i] = 0.0;
    }
    for (j = 0; j < msp; j++) {
      phi[(j + 1) * ns - 1] += params[j + mp + mq];
      for (i = 0; i < mp; i++) {
        phi[(j + 1) * ns + i] -= params[i] * params[j + mp + mq];
      }
    }
    for (j = 0; j < msq; j++) {
      theta[(j + 1) * ns - 1] += params[j + mp + mq + msp];
      for (i = 0; i < mq; i++) {
        theta[(j + 1) * ns + i] += params[i + mp] *
          params[j + mp + mq + msp];
      }
    }
  } else {
    for(i = 0; i < mp; i++) {
      phi[i] = params[i];
    }
    for(i = 0; i < mq; i++) {
      theta[i] = params[i + mp];
    }
  }
  // the output is written back to coef
  // ##########################################################
  // this is probabably a memory BUG for seasonal models
  for( i = 0; i < p; i++) {
    coef[i] = phi[i];
  }
  for( i = p; i < phi.size() + theta.size(); i++) {
    coef[i] = theta[i-p];
  }
}

void arima_transform_parameters2( std::vector<double> &coef,
                                  const arima_kind &arma,
                                  bool transform = true)
{
  // the coefficients are all 'packed in' inside coef - so we have
  // different types of coefficients. this tells us basically how many
  // of each type there are
  int mp = arma.p(), mq = arma.q(), msp = arma.P(), msq = arma.Q(), ns = arma.period();
  int n = mp + mq + msp + msq;

  int i, j, v;
  if (transform) {
    if (mp > 0) {
      parameter_transform2(coef, 0, mp);
    }
    v = mp + mq;
    if (msp > 0) {
      parameter_transform2(coef, v, msp);
    }
  }
  if ((msp + msq) > 0) {
    int p = mp + ns * msp;
    int q = mq + ns * msq;
    std::vector<double> params(p+q);
    int coef_old_size = coef.size();
    coef.resize(p+q);
    for( int i=coef_old_size; i < coef.size(); i++) {
      coef[i] = 0;
    }
    // /* expand out seasonal ARMA models */
    for (i = 0; i < coef.size(); i++) {
      params[i] = coef[i];
    }
    for (j = 0; j < msp; j++) {
      coef[(j + 1) * ns - 1] += params[j + mp + mq];
      for (i = 0; i < mp; i++) {
        coef[(j + 1) * ns + i] -= params[i] * params[j + mp + mq];
      }
    }
    for (j = p; j < p + msq; j++) {
      coef[(j + 1) * ns - 1] += params[j + mp + mq + msp];
      for (i = 0; i < mq; i++) {
        coef[(j + 1) * ns + i] += params[i + mp] *
          params[j + mp + mq + msp];
      }
    }
  }
}




#endif
