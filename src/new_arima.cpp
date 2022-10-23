#include <Rcpp.h>
using namespace Rcpp;
#include "cmath"
#include "utils.h"
#include "poly.h"

// originally partrans
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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
// [[Rcpp::export]]
std::vector<double> ts_conv(std::vector<double> & a, std::vector<double> & b)
{
  std::vector<double> result(a.size() + b. size() - 1);
  for (int i = 0; i < a.size(); i++) {
    for (int j = 0; j < b.size(); j++) {
      result[i + j] += a[i] * b[j];
    }
  }
  return result;
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
  // try to reuse a bit of memory :)
  // coef.resize(res.size());
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
  // perform the operation (ac+bd)/(c2 + d2) + (bc-ad)/(c2 + d2)i.
  // left = a+bi
  const double &a = left;
  const double &c = std::real(right);
  const double &d = std::imag(right);
  // the formula simplifies to (ac)/(c^2 + d^2) + (-ad)/(c^2 + d^2)i
  std::complex<double> res(
      // real part
      (left * c)/(pow(c,2) + pow(d,2)),
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

// SEXP ARIMA_transPars(SEXP sin, SEXP sarma, SEXP strans)
// {
//   int *arma = INTEGER(sarma), trans = asLogical(strans);
//   int mp = arma[0], mq = arma[1], msp = arma[2], msq = arma[3],
//                                                            ns = arma[4], i, j, p = mp + ns * msp, q = mq + ns * msq, v;
//   double *in = REAL(sin), *params = REAL(sin), *phi, *theta;
//   SEXP res, sPhi, sTheta;
//
//   PROTECT(res = allocVector(VECSXP, 2));
//   SET_VECTOR_ELT(res, 0, sPhi = allocVector(REALSXP, p));
//   SET_VECTOR_ELT(res, 1, sTheta = allocVector(REALSXP, q));
//   phi = REAL(sPhi);
//   theta = REAL(sTheta);
//   if (trans) {
//     int n = mp + mq + msp + msq;
//
//     params = (double *) R_alloc(n, sizeof(double));
//     for (i = 0; i < n; i++) params[i] = in[i];
//     if (mp > 0) partrans(mp, in, params);
//     v = mp + mq;
//     if (msp > 0) partrans(msp, in + v, params + v);
//   }
//   if (ns > 0) {
//     /* expand out seasonal ARMA models */
//     for (i = 0; i < mp; i++) phi[i] = params[i];
//     for (i = 0; i < mq; i++) theta[i] = params[i + mp];
//     for (i = mp; i < p; i++) phi[i] = 0.0;
//     for (i = mq; i < q; i++) theta[i] = 0.0;
//     for (j = 0; j < msp; j++) {
//       phi[(j + 1) * ns - 1] += params[j + mp + mq];
//       for (i = 0; i < mp; i++)
//         phi[(j + 1) * ns + i] -= params[i] * params[j + mp + mq];
//     }
//     for (j = 0; j < msq; j++) {
//       theta[(j + 1) * ns - 1] += params[j + mp + mq + msp];
//       for (i = 0; i < mq; i++)
//         theta[(j + 1) * ns + i] += params[i + mp] *
//           params[j + mp + mq + msp];
//     }
//   } else {
//     for (i = 0; i < mp; i++) phi[i] = params[i];
//     for (i = 0; i < mq; i++) theta[i] = params[i + mp];
//   }
//   UNPROTECT(1);
//   return res;
// }

// SEXP ARIMA_Invtrans(SEXP in, SEXP sarma)
// {
//   int *arma = INTEGER(sarma), mp = arma[0], mq = arma[1], msp = arma[2],
//                                                                     i, v, n = LENGTH(in);
//   SEXP y = allocVector(REALSXP, n);
//   double *raw = REAL(in), *new = REAL(y);
//
//   for(i = 0; i < n; i++) new[i] = raw[i];
//   if (mp > 0) invpartrans(mp, raw, new);
//   v = mp + mq;
//   if (msp > 0) invpartrans(msp, raw + v, new + v);
//   return y;
// }
//
// #define eps 1e-3
// SEXP ARIMA_Gradtrans(SEXP in, SEXP sarma)
// {
//   int *arma = INTEGER(sarma), mp = arma[0], mq = arma[1], msp = arma[2],
//                                                                     n = LENGTH(in);
//   SEXP y = allocMatrix(REALSXP, n, n);
//   double *raw = REAL(in), *A = REAL(y), w1[100], w2[100], w3[100];
//
//   for (int i = 0; i < n; i++)
//     for (int j = 0; j < n; j++)
//       A[i + j*n] = (i == j);
//   if(mp > 0) {
//     for (int i = 0; i < mp; i++) w1[i] = raw[i];
//     partrans(mp, w1, w2);
//     for (int i = 0; i < mp; i++) {
//       w1[i] += eps;
//       partrans(mp, w1, w3);
//       for (int j = 0; j < mp; j++) A[i + j*n] = (w3[j] - w2[j])/eps;
//       w1[i] -= eps;
//     }
//   }
//   if(msp > 0) {
//     int v = mp + mq;
//     for (int i = 0; i < msp; i++) w1[i] = raw[i + v];
//     partrans(msp, w1, w2);
//     for(int i = 0; i < msp; i++) {
//       w1[i] += eps;
//       partrans(msp, w1, w3);
//       for(int j = 0; j < msp; j++)
//         A[i + v + (j+v)*n] = (w3[j] - w2[j])/eps;
//       w1[i] -= eps;
//     }
//   }
//   return y;
// }


