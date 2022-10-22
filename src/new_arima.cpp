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


  // int na, nb, nab;
  // SEXP ab;
  // double *ra, *rb, *rab;

  // PROTECT(a = coerceVector(a, REALSXP));
  // PROTECT(b = coerceVector(b, REALSXP));
  // na = LENGTH(a);
  // nb = LENGTH(b);
  // nab = na + nb - 1;
  // PROTECT(ab = allocVector(REALSXP, nab));
  // ra = REAL(a);
  // rb = REAL(b);
  // rab = REAL(ab);
  std::vector<double> result(a.size() + b. size() - 1);
  // for (int i = 0; i < a.size() + b. size() - 1; i++) {
  //   rab[i] = 0.0;
  // }
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
  auto res = polyroot(coef);
  // find roots that are smaller than 1 - these are not inverse roots,
  // so they should be smaller, rather than larger than 1
  std::vector<bool> indices(res.size());
  int any = 0;
  bool temp = false;
  for( int i = 0; i < res.size(); i++ ) {
    temp = abs(res[i]) < 1;
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
  for( int i = 0; i < res.size(); i++ ) {
    // I made this work for complex numbers :)
    if( indices[i] ) {
      res[i] = 1/res[i];
    }
  }
  double x = 1;
  for( auto &root:res ) {
    // somewhat ugly
    // x <- c(x, 0) - c(0,x)/r
  }
  // finally take the real parts
  // c(Re(x[-1L]), rep.int(0, q - q0))
  // and return
}

