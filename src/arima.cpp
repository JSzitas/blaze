#include <Rcpp.h>
using namespace Rcpp;
#include "cmath"
#include "utils.h"
#include "poly.h"
// included mainly for isnan()
#include <math.h>,

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

// [[Rcpp::export]]
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

//TODO: test
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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



// SEXP
//   ARIMA_CSS(SEXP sy, SEXP sarma, SEXP sPhi, SEXP sTheta,
//             SEXP sncond, SEXP giveResid)
//   {
//     SEXP res, sResid = R_NilValue;
//     double ssq = 0.0, *y = REAL(sy), tmp;
//     double *phi = REAL(sPhi), *theta = REAL(sTheta), *w, *resid;
//     int n = LENGTH(sy), *arma = INTEGER(sarma), p = LENGTH(sPhi),
//       q = LENGTH(sTheta), ncond = asInteger(sncond);
//     int ns, nu = 0;
//     Rboolean useResid = asLogical(giveResid);
//
//     w = (double *) R_alloc(n, sizeof(double));
//     for (int l = 0; l < n; l++) w[l] = y[l];
//     for (int i = 0; i < arma[5]; i++)
//       for (int l = n - 1; l > 0; l--) w[l] -= w[l - 1];
//     ns = arma[4];
//     for (int i = 0; i < arma[6]; i++)
//       for (int l = n - 1; l >= ns; l--) w[l] -= w[l - ns];
//
//     PROTECT(sResid = allocVector(REALSXP, n));
//     resid = REAL(sResid);
//     if (useResid) for (int l = 0; l < ncond; l++) resid[l] = 0;
//
//     for (int l = ncond; l < n; l++) {
//       tmp = w[l];
//       for (int j = 0; j < p; j++) tmp -= phi[j] * w[l - j - 1];
//       for (int j = 0; j < min(l - ncond, q); j++)
//         tmp -= theta[j] * resid[l - j - 1];
//       resid[l] = tmp;
//       if (!ISNAN(tmp)) {
//         nu++;
//         ssq += tmp * tmp;
//       }
//     }
//     if (useResid) {
//       PROTECT(res = allocVector(VECSXP, 2));
//       SET_VECTOR_ELT(res, 0, ScalarReal(ssq / (double) (nu)));
//       SET_VECTOR_ELT(res, 1, sResid);
//       UNPROTECT(2);
//       return res;
//     } else {
//       UNPROTECT(1);
//       return ScalarReal(ssq / (double) (nu));
//     }
//   }


// SEXP
//   ARIMA_Like(SEXP sy, SEXP mod, SEXP sUP, SEXP giveResid)
//   {
//     SEXP sPhi = getListElement(mod, "phi"),
//       sTheta = getListElement(mod, "theta"),
//       sDelta = getListElement(mod, "Delta"),
//       sa = getListElement(mod, "a"),
//       sP = getListElement(mod, "P"),
//       sPn = getListElement(mod, "Pn");
//
//     if (TYPEOF(sPhi) != REALSXP || TYPEOF(sTheta) != REALSXP ||
//         TYPEOF(sDelta) != REALSXP || TYPEOF(sa) != REALSXP ||
//         TYPEOF(sP) != REALSXP || TYPEOF(sPn) != REALSXP)
//       error(_("invalid argument type"));
//
//     SEXP res, nres, sResid = R_NilValue;
//     int n = LENGTH(sy), rd = LENGTH(sa), p = LENGTH(sPhi),
//       q = LENGTH(sTheta), d = LENGTH(sDelta), r = rd - d;
//     double *y = REAL(sy), *a = REAL(sa), *P = REAL(sP), *Pnew = REAL(sPn);
//     double *phi = REAL(sPhi), *theta = REAL(sTheta), *delta = REAL(sDelta);
//     double sumlog = 0.0, ssq = 0, *anew, *mm = NULL, *M;
//     int nu = 0;
//     Rboolean useResid = asLogical(giveResid);
//     double *rsResid = NULL /* -Wall */;
//
//     anew = (double *) R_alloc(rd, sizeof(double));
//     M = (double *) R_alloc(rd, sizeof(double));
//     if (d > 0) mm = (double *) R_alloc(rd * rd, sizeof(double));
//
//     if (useResid) {
//       PROTECT(sResid = allocVector(REALSXP, n));
//       rsResid = REAL(sResid);
//     }
//
//     for (int l = 0; l < n; l++) {
//       for (int i = 0; i < r; i++) {
//         double tmp = (i < r - 1) ? a[i + 1] : 0.0;
//         if (i < p) tmp += phi[i] * a[0];
//         anew[i] = tmp;
//       }
//       if (d > 0) {
//         for (int i = r + 1; i < rd; i++) anew[i] = a[i - 1];
//         double tmp = a[0];
//         for (int i = 0; i < d; i++) tmp += delta[i] * a[r + i];
//         anew[r] = tmp;
//       }
//       if (l > asInteger(sUP)) {
//         if (d == 0) {
//           for (int i = 0; i < r; i++) {
//             double vi = 0.0;
//             if (i == 0) vi = 1.0; else if (i - 1 < q) vi = theta[i - 1];
//             for (int j = 0; j < r; j++) {
//               double tmp = 0.0;
//               if (j == 0) tmp = vi; else if (j - 1 < q) tmp = vi * theta[j - 1];
//               if (i < p && j < p) tmp += phi[i] * phi[j] * P[0];
//               if (i < r - 1 && j < r - 1) tmp += P[i + 1 + r * (j + 1)];
//               if (i < p && j < r - 1) tmp += phi[i] * P[j + 1];
//               if (j < p && i < r - 1) tmp += phi[j] * P[i + 1];
//               Pnew[i + r * j] = tmp;
//             }
//           }
//         } else {
//           /* mm = TP */
//           for (int i = 0; i < r; i++)
//             for (int j = 0; j < rd; j++) {
//               double tmp = 0.0;
//               if (i < p) tmp += phi[i] * P[rd * j];
//               if (i < r - 1) tmp += P[i + 1 + rd * j];
//               mm[i + rd * j] = tmp;
//             }
//             for (int j = 0; j < rd; j++) {
//               double tmp = P[rd * j];
//               for (int k = 0; k < d; k++)
//                 tmp += delta[k] * P[r + k + rd * j];
//               mm[r + rd * j] = tmp;
//             }
//             for (int i = 1; i < d; i++)
//               for (int j = 0; j < rd; j++)
//                 mm[r + i + rd * j] = P[r + i - 1 + rd * j];
//
//           /* Pnew = mmT' */
//           for (int i = 0; i < r; i++)
//             for (int j = 0; j < rd; j++) {
//               double tmp = 0.0;
//               if (i < p) tmp += phi[i] * mm[j];
//               if (i < r - 1) tmp += mm[rd * (i + 1) + j];
//               Pnew[j + rd * i] = tmp;
//             }
//             for (int j = 0; j < rd; j++) {
//               double tmp = mm[j];
//               for (int k = 0; k < d; k++)
//                 tmp += delta[k] * mm[rd * (r + k) + j];
//               Pnew[rd * r + j] = tmp;
//             }
//             for (int i = 1; i < d; i++)
//               for (int j = 0; j < rd; j++)
//                 Pnew[rd * (r + i) + j] = mm[rd * (r + i - 1) + j];
//           /* Pnew <- Pnew + (1 theta) %o% (1 theta) */
//           for (int i = 0; i <= q; i++) {
//             double vi = (i == 0) ? 1. : theta[i - 1];
//             for (int j = 0; j <= q; j++)
//               Pnew[i + rd * j] += vi * ((j == 0) ? 1. : theta[j - 1]);
//           }
//         }
//       }
//       if (!ISNAN(y[l])) {
//         double resid = y[l] - anew[0];
//         for (int i = 0; i < d; i++)
//           resid -= delta[i] * anew[r + i];
//
//         for (int i = 0; i < rd; i++) {
//           double tmp = Pnew[i];
//           for (int j = 0; j < d; j++)
//             tmp += Pnew[i + (r + j) * rd] * delta[j];
//           M[i] = tmp;
//         }
//
//         double gain = M[0];
//         for (int j = 0; j < d; j++) gain += delta[j] * M[r + j];
//         if(gain < 1e4) {
//           nu++;
//           ssq += resid * resid / gain;
//           sumlog += log(gain);
//         }
//         if (useResid) rsResid[l] = resid / sqrt(gain);
//         for (int i = 0; i < rd; i++)
//           a[i] = anew[i] + M[i] * resid / gain;
//         for (int i = 0; i < rd; i++)
//           for (int j = 0; j < rd; j++)
//             P[i + j * rd] = Pnew[i + j * rd] - M[i] * M[j] / gain;
//       } else {
//         for (int i = 0; i < rd; i++) a[i] = anew[i];
//         for (int i = 0; i < rd * rd; i++) P[i] = Pnew[i];
//         if (useResid) rsResid[l] = NA_REAL;
//       }
//     }
//
//     if (useResid) {
//       PROTECT(res = allocVector(VECSXP, 3));
//       SET_VECTOR_ELT(res, 0, nres = allocVector(REALSXP, 3));
//       REAL(nres)[0] = ssq;
//       REAL(nres)[1] = sumlog;
//       REAL(nres)[2] = (double) nu;
//       SET_VECTOR_ELT(res, 1, sResid);
//       UNPROTECT(2);
//       return res;
//     } else {
//       nres = allocVector(REALSXP, 3);
//       REAL(nres)[0] = ssq;
//       REAL(nres)[1] = sumlog;
//       REAL(nres)[2] = (double) nu;
//       return nres;
//     }
//   }
