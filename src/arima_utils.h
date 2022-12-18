#ifndef ARIMA_UTILS
#define ARIMA_UTILS

#include "cmath"
#include "utils/utils.h"
#include "utils/poly.h"
#include "structural_model.h"
#include "third_party/eigen.h"

// defines the arima structure
struct arima_kind{
  arima_kind(){};
  arima_kind( int p, int d, int q, int P, int D, int Q, int s_period ){
    this->arma_p= p;
    this->diff_d = d;
    this->arma_q = q;
    this->sarma_P = P;
    this->seas_diff_D = D;
    this->sarma_Q = Q;
    this->s_period = s_period;
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

// invert the parameter transformation
std::vector<double> inv_parameter_transform( std::vector<double> & phi ) {
  auto p = phi.size();
  std::vector<double> new_pars(p);
  std::vector<double> work(p);
  for(int j = 0; j < p; j++){
    // assigning work[j] = might be redundant - it gets reassigned anyways and
    // only holds intermediaries, but does not do anything for the result
    new_pars[j] = phi[j];
  }
  /* Run the Durbin-Levinson recursions backwards
   to find the PACF phi_{j.} from the autoregression coefficients */
  for(int j = p - 1; j > 0; j--) {
    // this allocation is redundant, we can just take reference to new_pars[j]
    // a = new_pars[j];
    for(int k = 0; k < j; k++) {
      work[k]  = (new_pars[k] + new_pars[j] * new_pars[j - k - 1]) / (1 - new_pars[j] * new_pars[j]);
    }
    for(int k = 0; k < j; k++) {
      new_pars[k] = work[k];
    }
  }
  // revert the tanh transform from earlier
  for(int j = 0; j < p; j++) {
    new_pars[j] = atanh(new_pars[j]);
  }
  return new_pars;
}

void inv_parameter_transform2( std::vector<double> & phi ) {

  auto p = phi.size();
  std::vector<double> temp(p*2);
  // std::vector<double> new_pars(p);
  // std::vector<double> work(p);
  for(int j = 0; j < p; j++){
    // assigning work[j] = might be redundant - it gets reassigned anyways and
    // only holds intermediaries, but does not do anything for the result
    temp[j] = phi[j];
  }
  /* Run the Durbin-Levinson recursions backwards
   to find the PACF phi_{j.} from the autoregression coefficients */
  for(int j = p - 1; j > 0; j--) {
    // this allocation is redundant, we can just take reference to new_pars[j]
    // a = new_pars[j];
    for(int k = 0; k < j; k++) {
      temp[p+k]  = (temp[k] + temp[j] * temp[j - k - 1]) / (1 - temp[j] * temp[j]);
    }
    for(int k = 0; k < j; k++) {
      temp[k] = temp[p+k];
    }
  }
  // revert the tanh transform from earlier
  for(int j = 0; j < p; j++) {
    phi[j] = atanh(temp[j]);
  }
}

bool ar_check( std::vector<double> & ar_coef, double tol = 0.0000001 ) {
  // check if all inverse coefficients are zero - in that case we return
  int p = 0;
  for( int i =0; i < ar_coef.size(); i++) {
    if( abs(-ar_coef[i]) > tol  ) {
      p = i;
    }
  }
  if( !p ) {
    return true;
  }
  // the polyroot part, e.g.
  // all(Mod(polyroot(c(1, -ar[1L:p]))) > 1)
  // goes here instead of the false
  // then the abs(/inverse/roots) must be all greater than > 1 for
  // this to lie within the unit circle
  std::vector<double> coef(p+1);
  coef[0] = 1;
  for( int i = 1; i< coef.size(); i++ ) {
    coef[i] = -ar_coef[i-1];
  }
  auto res = polyroot(coef);
  // check coefficients
  for( int i = 0; i < res.size(); i++ ) {
    if( abs(res[i]) <= 1 ) {
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

std::complex<double> operator /( double left, std::complex<double> right ) {
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
      (a*c)/(pow(c,2) + pow(d,2)),
      // complex part
      (-a*d)/(pow(c,2) + pow(d,2))
  );
  return res;
}

std::vector<double> ma_invert( std::vector<double> & ma_coef, double tol = 0.0000001 ) {
  auto q = ma_coef.size();
  // analogous to ar_check >> find latest/largest order nonzero coef
  int q0 = 0;
  for( int i =0; i < ma_coef.size(); i++) {
    if( abs(ma_coef[i]) > tol  ) {
      q0 = i;
    }
  }
  // if all coefs zero
  if(!q0) {
    return ma_coef;
  }
  // otherwise find roots
  std::vector<double> coef(q0+1);
  coef[0] = 1;
  for( int i = 1; i< coef.size(); i++ ) {
    coef[i] = ma_coef[i-1];
  }
  auto roots = polyroot(coef);
  // find roots that are smaller than 1 - these are not inverse roots,
  // so they should be smaller, rather than larger than 1
  std::vector<bool> indices(roots.size());
  int any = 0;
  bool temp = false;
  for( int i = 0; i < roots.size(); i++ ) {
    temp = abs(roots[i]) < 1;
    indices[i] = temp;
    // this simplifies the all(!ind) part
    any += temp;
  }
  // if all are > 1, we do not need to change them
  // as this implies they are valid inverse roots already - and there is no
  // need to invert
  if(!any) {
    return ma_coef;
  }
  // otherwise if there is only 1 valid coefficient
  if( q0 == 1 ) {
    coef.resize(q);
    // only invert the very first coefficient and otherwise set other coefficients
    // to zero
    coef[0] = 1/ma_coef[0];
    for( int i = 1; i < coef.size(); i++) {
      coef[i] = 0;
    }
    return coef;
  }
  // otherwise invert the roots which need to be inverted
  for( int i = 0; i < roots.size(); i++ ) {
    // I made this work for complex numbers :)
    if( indices[i] ) {
      roots[i] = 1/roots[i];
    }
  }
  std::vector<std::complex<double>> result(roots.size()+1);
  // set first new coef == 1
  result[0] = 1;
  std::vector<std::complex<double>> inv_temp( roots.size()+1);
  for( int i = 0; i < roots.size(); i++ ){
    // take the root
    inv_temp = result;
    for( int k = 1; k <= i+1; k++){
      result[k] = inv_temp[k] - (inv_temp[k-1]/roots[i]);
    }
  }
  coef.resize(ma_coef.size());
  for( int i = 1; i < result.size(); i++ ) {
    coef[i] = std::real(result[i]);
  }
  // if necessary pad with 0s
  for( int i = result.size(); i < q; i++ ) {
    coef[i] = 0;
  }
  return coef;
}

std::vector<std::complex<double>> invert_ma_coef_from_roots( std::vector<std::complex<double>> roots ) {
  std::vector<std::complex<double>> result(roots.size()+1);
  // set first new coef == 1
  result[0] = 1;
  std::vector<std::complex<double>> temp( roots.size()+1);
  for( int i = 0; i < roots.size(); i++ ){
    // take the root
    temp = result;
    for( int k = 1; k <= i+1; k++){
      result[k] = temp[k] - (temp[k-1]/roots[i]);
    }
  }
  return result;
}

// this just directly modifies coef
void arima_transform_parameters( std::vector<double> &coef,
                                 std::vector<int> &arma,
                                 bool transform = true)
{
  // the coefficients are all 'packed in' inside coef - so we have
  // different types of coefficients. this tells us basically how many
  // of each type there are
  int mp = arma[0], mq = arma[1], msp = arma[2], msq = arma[3], ns = arma[4];
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
  for( i = 0; i < p; i++) {
    coef[i] = phi[i];
  }
  for( i = p; i < phi.size() + theta.size(); i++) {
    coef[i] = theta[i-p];
  }
}

// this just directly modifies coef
void arima_transform_parameters( Eigen::VectorXd &coef,
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

void arima_transform_parameters( structural_model<double> model,
                                 arima_kind kind,
                                 bool transform = true)
{
  // the coefficients are all 'packed in' inside coef - so we have
  // different types of coefficients. this tells us basically how many
  // of each type there are
  int mp = kind.p(), mq = kind.q(), msp = kind.P(), msq = kind.Q(), ns = kind.period();
  int p = mp + ns * msp;
  int q = mq + ns * msq;

  int n = mp + mq + msp + msq;
  std::vector<double> params(n);
  int i, j, v;
  for (i = 0; i < mp; i++) {
    params[i] = model.phi[i];
  }
  for( i = mp; i < mq; i++) {
    params[i] = model.theta[i];
  }
  for (i = mq; i < msp; i++) {
    params[i] = model.phi[i];
  }
  for( i = msp; i < msq; i++) {
    params[i] = model.theta[i];
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
        temp[i-v] = params[i];
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
        model.theta[(j + 1) * ns + i] += params[i + mp] *
          params[j + mp + mq + msp];
      }
    }
  } else {
    for(i = 0; i < mp; i++) {
      model.phi[i] = params[i];
    }
    for(i = 0; i < mq; i++) {
      model.theta[i] = params[i + mp];
    }
  }
}


// this just directly modifies coef
void arima_transform_parameters2( std::vector<double> &coef,
                                  arima_kind &kind,
                                  bool transform = true)
{
  // the coefficients are all 'packed in' inside coef - so we have
  // different types of coefficients. this tells us basically how many
  // of each type there are
  int mp = kind.p(), mq = kind.q(), msp = kind.P(), msq = kind.Q(), ns = kind.period();
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
}


std::vector<double> arima_inverse_transform_parameters( std::vector<double> coef,
                                                        std::vector<int> &arma ) {

  // arma contains the arma structure - 3 leading numbers,
  // number of p parameters (arma[0])
  // number of q parameters (arma[1])
  // number of seasonal p parameters (arma[2])
  // for a function this small, it makes little sense to copy over
  std::vector<double> temp(arma[0]);
  if (arma[0] > 0) {
    for( int i = 0; i < arma[0]; i++ ) {
      temp[i] = coef[i];
    }
    temp = parameter_transform(temp);
    for( int i = 0; i < arma[0]; i++ ) {
      coef[i] = temp[i];
    }
  }

  if( arma[2] > 0 ) {
    temp.resize(arma[2]);
    for( int i = (arma[0] + arma[1]); i <coef.size(); i++ ) {
      temp[i] = coef[i];
    }
    temp = parameter_transform(temp);
    for( int i = (arma[0] + arma[1]); i < coef.size(); i++ ) {
      coef[i] = temp[i];
    }
  }
  return coef;
}

std::vector<double> arima_grad_transform( std::vector<double> &coef,
                                          std::vector<int> &arma,
                                          double eps = 1e-3) {
  int mp = arma[0], mq = arma[1], msp = arma[2], n = coef.size();
  std::vector<double> w1_temp(mp), w2_temp(mp), w3_temp(mp);

  std::vector<double> result(n*n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      result[i + j*n] = (i == j);
    }
  }
  if(mp > 0) {
    for (int i = 0; i < mp; i++) {
      w1_temp[i] = coef[i];
    }
    w2_temp = parameter_transform(w1_temp);
    for (int i = 0; i < mp; i++) {
      w1_temp[i] += eps;
      w3_temp = parameter_transform(w1_temp);
      for (int j = 0; j < mp; j++) {
        result[i + j*n] = (w3_temp[j] - w2_temp[j])/eps;
      }
      w1_temp[i] -= eps;
    }
  }
  if(msp > 0) {
    int v = mp + mq;
    w1_temp.resize(msp);
    w2_temp.resize(msp);
    w3_temp.resize(msp);
    for (int i = 0; i < msp; i++) {
      w1_temp[i] = coef[i + v];
    }
    w2_temp = parameter_transform(w1_temp);
    for(int i = 0; i < msp; i++) {
      w1_temp[i] += eps;
      w1_temp = parameter_transform(w3_temp);
      for(int j = 0; j < msp; j++) {
        result[i + v + (j+v)*n] = (w3_temp[j] - w2_temp[j])/eps;
      }
      w1_temp[i] -= eps;
    }
  }
  // result is basically a packed matrix
  return result;
}

#endif
