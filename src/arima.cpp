#include <Rcpp.h>
using namespace Rcpp;
#include "cmath"
#include "utils.h"
#include "poly.h"
// included mainly for isnan()
#include <math.h>

// originally partrans
std::vector<double> parameter_transform( std::vector<double> & coef ) {
  auto p = coef.size();
  int j, k;
  // if(p > 100) error(_("can only transform 100 pars in arima0"));
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
// static void invpartrans(int p, double *phi, double *new)
  // int j, k;
  auto p = phi.size();
  // double a;//, work[100];
  std::vector<double> new_pars(p);
  std::vector<double> work(p);
  // if(p > 100) error(_("can only transform 100 pars in arima0"));
  int j, k;
  for(j = 0; j < p; j++){
    // assigning work[j] = might be redundant - it gets reassigned anyways and
    // only holds intermediaries, but does not do anything for the result
    new_pars[j] = phi[j];
  }
  /* Run the Durbin-Levinson recursions backwards
   to find the PACF phi_{j.} from the autoregression coefficients */
  for(j = p - 1; j > 0; j--) {
    // this allocation is redundant, we can just take reference to new_pars[j]
    // a = new_pars[j];
    for(k = 0; k < j; k++) {
      work[k]  = (new_pars[k] + new_pars[j] * new_pars[j - k - 1]) / (1 - new_pars[j] * new_pars[j]);
    }
    for(k = 0; k < j; k++) {
      // this is a bit ugly - perhaps we can just std::move the whole thing afterwards?
      new_pars[k] = work[k];
    }
  }
  // revert the tanh transform from earlier
  for(j = 0;j < p; j++) {
    new_pars[j] = atanh(new_pars[j]);
  }
  return new_pars;
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

// this should probably be void f() and write the result back to coef
std::vector<double> arima_transform_parameters( std::vector<double> coef,
                                                std::vector<int> arma,
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
  // the output is a 2 element list
  std::vector<double> result(phi.size()+theta.size());
  for( i = 0; i < p; i++) {
    result[i] = phi[i];
  }
  for( i = p; i < phi.size() + theta.size(); i++) {
    result[i] = theta[i-p];
  }
  return result;
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


/* arma is p, q, sp, sq, ns, d, sd
 * Note that this function is very similar to the one that follows it
 * the main difference is in what they return -
 */
double arima_css_ssq( std::vector<double> & y,
                      std::vector<int> & arma,
                      std::vector<double> & phi,
                      std::vector<double> & theta,
                      int n_cond )
{
  // SEXP res, sResid = R_NilValue;
  double ssq = 0.0, tmp = 0.;
  int n = y.size(), p = phi.size(), q = theta.size();
  int ns, nu = 0;
  // Rboolean useResid = asLogical(giveResid);
  std::vector<double> w(n);
  // w = (double *) R_alloc(n, sizeof(double));
  for (int l = 0; l < n; l++) {
    w[l] = y[l];
  }
  // regular differencing, as far as I can tell :)
  for (int i = 0; i < arma[5]; i++) {
    for (int l = n - 1; l > 0; l--) {
      w[l] -= w[l - 1];
    }
  }
  ns = arma[4];
  // seasonal differencing, as far as I can tell :)
  for (int i = 0; i < arma[6]; i++) {
    for (int l = n - 1; l >= ns; l--) {
      w[l] -= w[l - ns];
    }
  }
  // prepare the residuals
  std::vector<double> resid(n);
  for (int l = 0; l < n_cond; l++) {
    resid[l] = 0;
  }

  for (int l = n_cond; l < n; l++) {
    tmp = w[l];
    for (int j = 0; j < p; j++) {
      tmp -= phi[j] * w[l - j - 1];
    }
    for (int j = 0; j < min(l - n_cond, q); j++) {
      tmp -= theta[j] * resid[l - j - 1];
    }
    resid[l] = tmp;
    if (!isnan(tmp)) {
      nu++;
      ssq += tmp * tmp;
    }
  }
  return ssq/nu;
}

std::vector<double> arima_css_resid( std::vector<double> & y,
                                     std::vector<int> & arma,
                                     std::vector<double> & phi,
                                     std::vector<double> & theta,
                                     int n_cond )
{
  // SEXP res, sResid = R_NilValue;
  double ssq = 0.0, tmp = 0.;
  int n = y.size(), p = phi.size(), q = theta.size();
  int ns, nu = 0;
  // Rboolean useResid = asLogical(giveResid);
  std::vector<double> w(n);
  // w = (double *) R_alloc(n, sizeof(double));
  for (int l = 0; l < n; l++) {
    w[l] = y[l];
  }
  // regular differencing, as far as I can tell :)
  for (int i = 0; i < arma[5]; i++) {
    for (int l = n - 1; l > 0; l--) {
      w[l] -= w[l - 1];
    }
  }
  ns = arma[4];
  // seasonal differencing, as far as I can tell :)
  for (int i = 0; i < arma[6]; i++) {
    for (int l = n - 1; l >= ns; l--) {
      w[l] -= w[l - ns];
    }
  }
  // prepare the residuals
  std::vector<double> resid(n);
  for (int l = 0; l < n_cond; l++) {
    resid[l] = 0;
  }

  for (int l = n_cond; l < n; l++) {
    tmp = w[l];
    for (int j = 0; j < p; j++) {
      tmp -= phi[j] * w[l - j - 1];
    }
    for (int j = 0; j < min(l - n_cond, q); j++) {
      tmp -= theta[j] * resid[l - j - 1];
    }
    resid[l] = tmp;
    if (!isnan(tmp)) {
      nu++;
      ssq += tmp * tmp;
    }
  }
  return resid;
}

// potentially use something that represents the model??
// pass a reference to residuals - we will work on them
// directly
// std::vector<double> & residuals,
// we will likewise modify these items
// these are the structural model matrices
// [[Rcpp::export]]
std::vector<double> arima_likelihood( std::vector<double> & y,
                                      std::vector<double> & phi,
                                      std::vector<double> & theta,
                                      std::vector<double> & delta,
                                      std::vector<double> & a,
                                      std::vector<double> & P,
                                      std::vector<double> & Pn ) {
  // define integers needed for further processing - these are mostly used
  // for indexing and offsetting
  int n = y.size(), rd = a.size(), p = phi.size(),
    q = theta.size(), d = delta.size(), r = rd - d,
    nu = 0;

  // define data structures needed for computation intermediaries
  std::vector<double> anew(rd);
  std::vector<double> M(rd);
  std::vector<double> Pnew = Pn;
  // this is only needed if we have any deltas
  std::vector<double> mm(0);
  if(d > 0) {
    mm.resize(rd*rd);
  }

  double tmp, vi, resid, gain, sumlog = 0, ssq = 0;
  for (int l = 0; l < n; l++) {
    for (int i = 0; i < r; i++) {
      tmp = (i < r - 1) ? a[i + 1] : 0.0;
      if (i < p) {
        tmp += phi[i] * a[0];
      }
      anew[i] = tmp;
    }
    if (d > 0) {
      for (int i = r + 1; i < rd; i++) {
        anew[i] = a[i - 1];
      }
      tmp = a[0];
      for (int i = 0; i < d; i++) {
        tmp += delta[i] * a[r + i];
      }
      anew[r] = tmp;
    }
    // only if we are past the first observation
    if (l > 0) {
      // if we have any thetas
      if (d == 0) {
        for (int i = 0; i < r; i++) {
          vi = 0.0;
          // presumably leading coefficient
          if (i == 0) {
            vi = 1.0;
          } else if (i - 1 < q) {
            vi = theta[i - 1];
          }
          for (int j = 0; j < r; j++) {
            tmp = 0.0;
            if (j == 0) {
              tmp = vi;
            } else if(j - 1 < q) {
              tmp = vi * theta[j - 1];
            }
            if (i < p && j < p) {
              tmp += phi[i] * phi[j] * P[0];
            }
            if (i < r - 1 && j < r - 1) {
              tmp += P[i + 1 + r * (j + 1)];
            }
            if (i < p && j < r - 1) {
              tmp += phi[i] * P[j + 1];
            }
            if (j < p && i < r - 1) {
              tmp += phi[j] * P[i + 1];
            }
            // update new P matrix with appropriate entry
            Pnew[i + r * j] = tmp;
          }
        }
    } else {
      /* mm = TP */
      for (int i = 0; i < r; i++) {
        for (int j = 0; j < rd; j++) {
          tmp = 0.0;
          if (i < p) {
            tmp += phi[i] * P[rd * j];
          }
          if (i < r - 1) {
            tmp += P[i + 1 + rd * j];
          }
          mm[i + rd * j] = tmp;
        }
      }
      for (int j = 0; j < rd; j++) {
        tmp = P[rd * j];
        for (int k = 0; k < d; k++) {
          tmp += delta[k] * P[r + k + rd * j];
        }
        mm[r + rd * j] = tmp;
      }
      for (int i = 1; i < d; i++) {
        for (int j = 0; j < rd; j++) {
          mm[r + i + rd * j] = P[r + i - 1 + rd * j];
        }
      }
      /* Pnew = mmT' */
      for (int i = 0; i < r; i++) {
        for (int j = 0; j < rd; j++) {
          tmp = 0.0;
          if (i < p) {
            tmp += phi[i] * mm[j];
          }
          if (i < r - 1) {
            tmp += mm[rd * (i + 1) + j];
          }
          Pnew[j + rd * i] = tmp;
        }
      }
      for (int j = 0; j < rd; j++) {
        tmp = mm[j];
        for (int k = 0; k < d; k++) {
          tmp += delta[k] * mm[rd * (r + k) + j];
        }
        Pnew[rd * r + j] = tmp;
      }
      for (int i = 1; i < d; i++) {
        for (int j = 0; j < rd; j++) {
          Pnew[rd * (r + i) + j] = mm[rd * (r + i - 1) + j];
        }
      }
      /* Pnew <- Pnew + (1 theta) %o% (1 theta) */
      for (int i = 0; i <= q; i++) {
        vi = (i == 0) ? 1. : theta[i - 1];
        for (int j = 0; j <= q; j++) {
          Pnew[i + rd * j] += vi * ((j == 0) ? 1. : theta[j - 1]);
          }
        }
      }
    }
    if (!isnan(y[l])) {
      resid = y[l] - anew[0];
      for (int i = 0; i < d; i++) {
        resid -= delta[i] * anew[r + i];
      }
      for (int i = 0; i < rd; i++) {
        tmp = Pnew[i];
        for (int j = 0; j < d; j++) {
          tmp += Pnew[i + (r + j) * rd] * delta[j];
        }
        M[i] = tmp;
      }
      gain = M[0];
      for (int j = 0; j < d; j++) {
        gain += delta[j] * M[r + j];
      }
      // if gain is reasonable, update nu, residual sum of squares and
      // sum of log gain
      if(gain < 1e4) {
        nu++;
        ssq += resid * resid / gain;
        sumlog += log(gain);
      }
      // you would normally update the residuals here
      // also, you get to update a, and P - this should change them by
      // reference (so that you do not have to return them)
      for (int i = 0; i < rd; i++) {
        a[i] = anew[i] + M[i] * resid / gain;
      }
      for (int i = 0; i < rd; i++) {
        for (int j = 0; j < rd; j++) {
          P[i + j * rd] = Pnew[i + j * rd] - M[i] * M[j] / gain;
        }
      }
    } else {
      for (int i = 0; i < rd; i++) {
        a[i] = anew[i];
      }
      for (int i = 0; i < rd * rd; i++) {
        P[i] = Pnew[i];
      }
      // if you were updating residuals, here you would put in an 'NA' or NaN
    }
  }
  // finally, return
  std::vector<double> res{ssq, sumlog, (double) nu};
  return res;
}

// [[Rcpp::export]]
std::vector<double> make_delta( int n_diff,
                                int seas_period = 1,
                                int n_seas_diff = 0 ) {
  int diff_size = n_diff+1;
  std::vector<double> a(diff_size + (n_seas_diff * seas_period));
  a[0] = 1;
  std::vector<double> temp(diff_size + (n_seas_diff * seas_period));
  for( int k = 0; k < n_diff; k++) {
    for (int i = 0; i <= k ; i++) {
      // the array extend is always 2, hence we can always just do these two operations
      // first this is temp[i+0] += a[i] * 1;
      temp[i] += a[i]; // * 1
      // and this is the other operation
      // a[i] * -1 == -= a[i];
      temp[i+1] -= a[i];
    }
    // move all of the elements of temp to a - but temp has constant size,
    // so we can just use k+2
    for( int i = 0; i < k+2; i++ ) {
      a[i] = std::move(temp[i]);
    }
    std::fill(temp.begin(), temp.end(), 0);
  }
  // seasonal differences:
  for( int k = 0; k < n_seas_diff; k++) {
    for (int i = 0; i < diff_size + (k*seas_period); i++) {
      /* we know that this only operates on the first and last element of
       * the vector - it adds a[i] * 1 to the first element and adds
       * a[i] * -1 to the last - which is effectively adding and subtracting
       * a[i] at various points, i.e.
       * temp[i+0] += a[i] * 1; */
      temp[i] += a[i];
      // and temp[i+seas_period] += a[i] * -1;
      temp[i + seas_period] -= a[i];
      }
    for( int i = 0; i < temp.size(); i++ ) {
      a[i] = std::move(temp[i]);
    }
    std::fill(temp.begin(), temp.end(), 0);
  }
  // remove leading coefficient and flip signs
  pop_front(a);
  for( unsigned long long i = 0; i < a.size(); i++ ) {
    a[i] = -a[i];
  }
  return a;
}


