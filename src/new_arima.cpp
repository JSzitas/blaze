#include <Rcpp.h>
using namespace Rcpp;
#include "cmath"
#include "utils.h"

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
  return false;
}



maInvert <- function(ma) {
  q <- length(ma)
  q0 <- max(which(c(1, ma) != 0)) - 1L
  if (!q0)
    return(ma)
    roots <- polyroot(c(1, ma[1L:q0]))
    ind <- Mod(roots) < 1
  if (all(!ind))
    return(ma)
    if (q0 == 1)
      return(c(1/ma[1L], rep.int(0, q - q0)))
      roots[ind] <- 1/roots[ind]
    x <- 1
    for (r in roots) x <- c(x, 0) - c(0, x)/r
      c(Re(x[-1L]), rep.int(0, q - q0))
}






