#include <Rcpp.h>
using namespace Rcpp;

#include "eigen.hpp"
#include "optim.hpp"


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


inline double arima_css_wrap(const Eigen::VectorXd& vals_inp, Eigen::VectorXd* grad_out, void* opt_data)
  {
    double obj_val = vals_inp.dot(vals_inp);

    if (grad_out) {
      *grad_out = 2.0*vals_inp;
    }

    return obj_val;
  }



// [[Rcpp::export]]
void test_optim(std::vector<double> x) {

  Eigen::VectorXd y = Eigen::VectorXd::Ones(test_dim);

  return x * 2;
}
