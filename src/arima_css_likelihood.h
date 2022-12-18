#ifndef ARIMA_CSS_LIKELIHOOD
#define ARIMA_CSS_LIKELIHOOD

#include "arima_utils.h"
#include "structural_model.h"

/* arma is p, q, sp, sq, ns, d, sd
 * Note that this function is very similar to the one that follows it
 * the main difference is in what they return -
 */
double arima_css_ssq( std::vector<double> & y,
                      structural_model<double> &model,
                      arima_kind &kind,
                      int n_cond )
{
  double ssq = 0.0, tmp = 0.;
  int n = y.size(), p = model.phi.size(), q = model.theta.size();
  int ns, nu = 0;
  std::vector<double> w(n);
  // w = (double *) R_alloc(n, sizeof(double));
  for (int l = 0; l < n; l++) {
    w[l] = y[l];
  }
  // regular differencing, as far as I can tell :)
  for (int i = 0; i < kind.d(); i++) {
    for (int l = n - 1; l > 0; l--) {
      w[l] -= w[l - 1];
    }
  }
  ns = kind.period();
  // seasonal differencing, as far as I can tell :)
  for (int i = 0; i < kind.D(); i++) {
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
      tmp -= model.phi[j] * w[l - j - 1];
    }
    for (int j = 0; j < min(l - n_cond, q); j++) {
      tmp -= model.theta[j] * resid[l - j - 1];
    }
    resid[l] = tmp;
    if (!isnan(tmp)) {
      nu++;
      ssq += tmp * tmp;
    }
  }
  return ssq/nu;
}

double arima_css_ssq( std::vector<double> & y,
                      const Eigen::VectorXd & pars,
                      arima_kind &kind,
                      int n_cond )
{
  double ssq = 0.0, tmp = 0.;
  int n = y.size(), p = kind.p() + kind.P(), q = kind.q() + kind.Q();
  int ns, nu = 0;
  std::vector<double> w(n);
  // w = (double *) R_alloc(n, sizeof(double));
  for (int l = 0; l < n; l++) {
    w[l] = y[l];
  }
  // regular differencing, as far as I can tell :)
  for (int i = 0; i < kind.d(); i++) {
    for (int l = n - 1; l > 0; l--) {
      w[l] -= w[l - 1];
    }
  }
  ns = kind.period();
  // seasonal differencing, as far as I can tell :)
  for (int i = 0; i < kind.D(); i++) {
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
      tmp -= pars[j] * w[l - j - 1];
    }
    for (int j = p; j < min(l - n_cond, q); j++) {
      tmp -= pars[j] * resid[l - j - 1];
    }
    resid[l] = tmp;
    if (!isnan(tmp)) {
      nu++;
      ssq += tmp * tmp;
    }
  }
  return ssq/nu;
}

// const std::vector<double> & y,

double arima_css_ssq( const Eigen::VectorXd & y,
                      const Eigen::VectorXd & pars,
                      const arima_kind &kind,
                      const int n_cond )
{
  double ssq = 0.0, tmp = 0.0;
  int n = y.size(), p = kind.p() + kind.P(), q = kind.q() + kind.Q();
  int ns, nu = 0;
  // std::cout << "ncond: " << n_cond << std::endl;
  std::vector<double> w(n);
  // w = (double *) R_alloc(n, sizeof(double));
  for (int l = 0; l < n; l++) {
    w[l] = y[l];
  }
  // regular differencing, as far as I can tell :)
  for (int i = 0; i < kind.d(); i++) {
    for (int l = n - 1; l > 0; l--) {
      w[l] -= w[l - 1];
    }
  }
  ns = kind.period();
  // seasonal differencing, as far as I can tell :)
  for (int i = 0; i < kind.D(); i++) {
    for (int l = n - 1; l >= ns; l--) {
      w[l] -= w[l - ns];
    }
  }
  // prepare the residuals
  std::vector<double> resid(n);
  // for (int l = 0; l < n_cond; l++) {
  //   resid[l] = 0;
  // }
  // print_vector(w);

  int ma_offset;
  for (int l = n_cond; l < n; l++) {
    ma_offset = min(l - n_cond, q);
    tmp = w[l];
    // std::cout << "tmp: " << tmp << std::endl;
    for (int j = 0; j < p; j++) {
      tmp -= pars[j] * w[l - j - 1];
    }
    // to offset that this is all in one vector, we need to
    // start at p and go to p + q
    for (int j = 0; j < ma_offset; j++) {
      tmp -= pars[p+j] * resid[l - j - 1];
    }
    resid[l] = tmp;
    if (!isnan(tmp)) {
      nu++;
      ssq += tmp * tmp;
      // std::cout << " tmp post: " << tmp << " ssq " << ssq << std::endl;
    }
  }
  // std::cout << "SSQ: " << ssq << " nu: " << nu << std::endl;
  return ssq/nu;
}

std::vector<double> arima_css_resid( std::vector<double> & y,
                                     structural_model<double> &model,
                                     arima_kind &kind,
                                     int n_cond )
{
  double ssq = 0.0, tmp = 0.;
  int n = y.size(), p = model.phi.size(), q = model.theta.size();
  int ns, nu = 0;
  std::vector<double> w(n);
  for (int l = 0; l < n; l++) {
    w[l] = y[l];
  }
  // regular differencing
  for (int i = 0; i < kind.d(); i++) {
    for (int l = n - 1; l > 0; l--) {
      w[l] -= w[l - 1];
    }
  }
  ns = kind.period();
  // seasonal differencing
  for (int i = 0; i < kind.D(); i++) {
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
      tmp -= model.phi[j] * w[l - j - 1];
    }
    for (int j = 0; j < min(l - n_cond, q); j++) {
      tmp -= model.theta[j] * resid[l - j - 1];
    }
    resid[l] = tmp;
    if (!isnan(tmp)) {
      nu++;
      ssq += tmp * tmp;
    }
  }
  return resid;
}

std::vector<double> arima_likelihood( std::vector<double> & y,
                                      structural_model<double> &model ) {
                                        // define integers needed for further processing - these are mostly used
                                        // for indexing and offsetting
                                        int n = y.size(), rd = model.a.size(), p = model.phi.size(),
                                          q = model.theta.size(), d = model.delta.size(), r = rd - d,
                                          nu = 0;

  // define data structures needed for computation intermediaries
  std::vector<double> anew(rd);
  std::vector<double> M(rd);
  std::vector<double> Pnew = model.Pn;
  // this is only needed if we have any deltas
  std::vector<double> mm(0);
  if(d > 0) {
    mm.resize(rd*rd);
  }
  double tmp, vi, resid, gain, sumlog = 0, ssq = 0;
  for (int l = 0; l < n; l++) {
    for (int i = 0; i < r; i++) {
      tmp = (i < r - 1) ? model.a[i + 1] : 0.0;
      if (i < p) {
        tmp += model.phi[i] * model.a[0];
      }
      anew[i] = tmp;
    }
    if (d > 0) {
      for (int i = r + 1; i < rd; i++) {
        anew[i] = model.a[i - 1];
      }
      tmp = model.a[0];
      for (int i = 0; i < d; i++) {
        tmp += model.delta[i] * model.a[r + i];
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
            vi = model.theta[i - 1];
          }
          for (int j = 0; j < r; j++) {
            tmp = 0.0;
            if (j == 0) {
              tmp = vi;
            } else if(j - 1 < q) {
              tmp = vi * model.theta[j - 1];
            }
            if (i < p && j < p) {
              tmp += model.phi[i] * model.phi[j] * model.P[0];
            }
            if (i < r - 1 && j < r - 1) {
              tmp += model.P[i + 1 + r * (j + 1)];
            }
            if (i < p && j < r - 1) {
              tmp += model.phi[i] * model.P[j + 1];
            }
            if (j < p && i < r - 1) {
              tmp += model.phi[j] * model.P[i + 1];
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
              tmp += model.phi[i] * model.P[rd * j];
            }
            if (i < r - 1) {
              tmp += model.P[i + 1 + rd * j];
            }
            mm[i + rd * j] = tmp;
          }
        }
        for (int j = 0; j < rd; j++) {
          tmp = model.P[rd * j];
          for (int k = 0; k < d; k++) {
            tmp += model.delta[k] * model.P[r + k + rd * j];
          }
          mm[r + rd * j] = tmp;
        }
        for (int i = 1; i < d; i++) {
          for (int j = 0; j < rd; j++) {
            mm[r + i + rd * j] = model.P[r + i - 1 + rd * j];
          }
        }
        /* Pnew = mmT' */
        for (int i = 0; i < r; i++) {
          for (int j = 0; j < rd; j++) {
            tmp = 0.0;
            if (i < p) {
              tmp += model.phi[i] * mm[j];
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
            tmp += model.delta[k] * mm[rd * (r + k) + j];
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
          vi = (i == 0) ? 1. : model.theta[i - 1];
          for (int j = 0; j <= q; j++) {
            Pnew[i + rd * j] += vi * ((j == 0) ? 1. : model.theta[j - 1]);
          }
        }
      }
    }
    if (!isnan(y[l])) {
      resid = y[l] - anew[0];
      for (int i = 0; i < d; i++) {
        resid -= model.delta[i] * anew[r + i];
      }
      for (int i = 0; i < rd; i++) {
        tmp = Pnew[i];
        for (int j = 0; j < d; j++) {
          tmp += Pnew[i + (r + j) * rd] * model.delta[j];
        }
        M[i] = tmp;
      }
      gain = M[0];
      for (int j = 0; j < d; j++) {
        gain += model.delta[j] * M[r + j];
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
        model.a[i] = anew[i] + M[i] * resid / gain;
      }
      for (int i = 0; i < rd; i++) {
        for (int j = 0; j < rd; j++) {
          model.P[i + j * rd] = Pnew[i + j * rd] - M[i] * M[j] / gain;
        }
      }
    } else {
      for (int i = 0; i < rd; i++) {
        model.a[i] = anew[i];
      }
      for (int i = 0; i < rd * rd; i++) {
        model.P[i] = Pnew[i];
      }
      // if you were updating residuals, here you would put in an 'NA' or NaN
    }
  }
  // finally, return
  std::vector<double> res{ssq, sumlog, (double) nu};
  return res;
}

#endif
