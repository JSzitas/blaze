#ifndef STATE_SPACE
#define STATE_SPACE

#include "arima/structures/arima_kind.h"
#include "arima/structures/structural_model.h"
#include "utils/utils.h"
#include "vector"
#include "arima/solvers/initializers.h"
#include "arima/utils/delta.h"

template <typename U = double> struct forecast_result {
  forecast_result<U>( std::vector<U> &forecasts,
                      std::vector<U> &std_errs ) :
  forecast(std::move(forecasts)), std_err(std::move(std_errs)){}
  forecast_result<U>( std::vector<U> &&forecasts,
                      std::vector<U> &&std_errs ) :
  forecast(std::move(forecasts)), std_err(std::move(std_errs)){}
  forecast_result<U>(size_t h) {
    this->forecast = std::vector<U>(h);
    this->std_err = std::vector<U>(h);
  }
  void add(size_t i, U fcst, U se) {
    this->forecast[i] = fcst;
    this->std_err[i] = se;
  }
  std::vector<U> forecast;
  std::vector<U> std_err;
};

/* Forecasts based on state space representation of ARIMA via
 * the kalman filter.
 * TODO: add an option to update model.a and model.P while this runs
 */
template < typename U = double>
forecast_result<U> kalman_forecast(const size_t n_ahead,
                                   structural_model<U> &model) {
  size_t p = model.a.size();
  std::vector<double> anew = model.a;
  std::vector<double> a = model.a;
  std::vector<double> Pnew(p * p);
  std::vector<double> mm(p * p);
  std::vector<double> P = model.P;

  std::vector<double> forecasts(n_ahead);
  std::vector<double> standard_errors(n_ahead);
  double fc = 0.0, tmp = 0.0;
  for (size_t l = 0; l < n_ahead; l++) {
    fc = 0.0;
    for (size_t i = 0; i < p; i++) {
      tmp = 0.0;
      for (size_t k = 0; k < p; k++) {
        tmp += model.T[i + p * k] * a[k];
      }
      anew[i] = tmp;
      fc += tmp * model.Z[i];
    }
    for (size_t i = 0; i < p; i++) {
      a[i] = anew[i];
    }
    forecasts[l] = fc;
    for (size_t i = 0; i < p; i++) {
      for (size_t j = 0; j < p; j++) {
        tmp = 0.0;
        for (size_t k = 0; k < p; k++) {
          tmp += model.T[i + p * k] * P[k + p * j];
        }
        mm[i + p * j] = tmp;
      }
    }
    for (size_t i = 0; i < p; i++) {
      for (size_t j = 0; j < p; j++) {
        tmp = model.V[i + p * j];
        for (size_t k = 0; k < p; k++) {
          tmp += mm[i + p * k] * model.T[j + p * k];
        }
        Pnew[i + p * j] = tmp;
      }
    }
    tmp = model.h;
    for (size_t i = 0; i < p; i++) {
      for (size_t j = 0; j < p; j++) {
        tmp += model.Z[i] * model.Z[j] * Pnew[i + j * p];
        P[i + j * p] = Pnew[i + j * p];
      }
    }
    standard_errors[l] = tmp;
  }
  return forecast_result<U>( forecasts, standard_errors );
}

/* originally an R function - this creates the arima model in state space
 * representation from some values of phi, theta and delta, where phi
 * and theta are embedded in a single vector(std::vector/Eigen::Vector)
 */
template <typename C, typename U = double>
structural_model<U> make_arima( const C &coef,
                                const arima_kind &kind,
                                const U kappa = 1000000,
                                const SSinit state_init = Gardner,
                                const U tol = 1e-9) {
  std::vector<U> delta = make_delta<U>(kind.d(), kind.period(), kind.D());
  const size_t p = kind.p() + (kind.P() * kind.period());
  const size_t q = kind.q() + (kind.Q() * kind.period());
  const size_t r = max(p, q + 1), d = delta.size(), rd = r + d;

  std::vector<U> phi(p);
  std::vector<U> theta(q + max(r - 1 - q, 0));
  // copz out elements of coef into phi and theta
  for( size_t i = 0; i < p; i++) {
    phi[i] = coef[i];
  }
  for( size_t i = p; i < p + q; i++) {
    theta[i-p] = coef[i];
  }
  size_t i, j;
  std::vector<U> Z(rd);
  Z[0] = 1;
  for (i = 1; i < r - 1; i++) {
    Z[i] = 0;
  }
  j = 0;
  for (i = r; i < rd; i++) {
    Z[i] = delta[j];
    j++;
  }
  std::vector<U> T(rd * rd);
  if (p > 0) {
    for (i = 0; i < p; i++) {
      T[i] = phi[i];
    }
  }
  if (r > 1L) {
    /* set '2nd diagonal' elements to 1. since this is hard to understand,
     * here are some examples
     * for a matrix with *rd* == 5, transform it:
     *
     *   (input)          (rd == r)         (rd > r)
     *   x x x x x   =>   x 1 x x x  //  => x 1 x x x
     *   x x x x x   =>   x x 1 x x  //  => x x 1 x x
     *   x x x x x   =>   x x x 1 x  //  => x x x 1 x
     *   x x x x x   =>   x x x x 1  //  => x x x x x
     *   x x x x x   =>   x x x x x  //  => x x x x x
     */
    for (i = 0; i < r - 1; i++) {
      T[(rd * (i + 1)) + i] = 1;
    }
  }
  if (d > 0) {
    // replace row r+1 in R (or row r for us)
    // with whatever is in Z
    for (j = 0; j < rd; j++) {
      T[(j * rd) + r] = Z[j];
    }
    // if there are more than 1 differences d
    if (d > 1) {
      /* start replacing at r and continue until you get to r + d (-1?)
       * replace r + 2:d with 1 - this is similar as above, but it accounts
       * for the first difference differently(that is taken care of in the above
       * code, in the j < rd loop). here we are taking care of the other
       * differences, so if we have 3 differences, we will only be doing 2
       * replacements. these happen after the for a matrix with *rd* == 5, with
       * *d* == 3, transform it: (input)          (d == 3) x x x x x   =>   x x
       * x x x x x x x x   =>   x x x x x x x x x x   =>   x x x x x x x x x x
       * =>   x x 1 x x x x x x x   =>   x x x 1 x
       */
      for (i = r; i < rd - 1; i++) {
        T[((rd + 1) * i) + 1] = 1;
      }
    }
  }
  // this is R <- c(1, theta, rep.int(0, d))
  // we can skip the d part as vectors are 0 initialized.
  std::vector<U> R(1 + theta.size() + d);
  R[0] = 1;
  for (i = 1; i < theta.size() + 1; i++) {
    R[i] = theta[i - 1];
  }
  std::vector<U> V(R.size() * R.size());
  // here we do an outer product, ie: V <- R %o% R
  size_t mat_p = 0;
  for (i = 0; i < R.size(); i++) {
    for (j = 0; j < R.size(); j++) {
      V[mat_p] = R[i] * R[j];
      mat_p++;
    }
  }
  U h = 0;
  std::vector<U> a(rd);
  std::vector<U> P(rd * rd);
  std::vector<U> Pn(rd * rd);
  if (r > 1) {
    // for storing initialization results
    std::vector<U> temp(r * r);
    switch (state_init) {
    case Gardner:
      temp = std::move(get_Q0(phi, theta));
      break;
    case Rossignol:
      temp = std::move(get_Q0_rossignol(phi, theta));
      break;
    };
    /* update a block of first r rows and columns i.e. if we have a 5x5 Pn
     * matrix, and r == 3, then we update the highlighted parts: (input)
     * (updated) x x x x x   =>    y y y|x x x x x x x   =>    y y y|x x x x x x
     * x   =>    y y y|x x
     *                     _____
     *   x x x x x   =>    x x x x x
     *   x x x x x   =>    x x x x x
     */
    mat_p = 0;
    for (j = 0; j < r; j++) {
      for (i = 0; i < r; i++) {
        Pn[(j * rd) + i] = std::move(temp[mat_p]);
        mat_p++;
      }
    }
  } else {
    Pn[0] = (p > 0) * (1 / (1 - pow(phi[0], 2))) + (p == 0);
  }
  if (d > 0L) {
    /* update diagonal elements which come after the coefficients -
     * diagonal entries between r and rd - with kappa
     */
    for (i = r; i < rd; i++) {
      for (j = r; j < rd; j++) {
        // to only update diagonal elements check that we are on the diagonal
        // otherwise we have a zero - as intended
        Pn[(j * rd) + i] = (i == j) * kappa;
      }
    }
  }
  structural_model<U> res(phi, theta, delta, Z, a, P, T, V, h, Pn);
  return res;
}

template <typename C, typename U = double>
structural_model<U> make_arima( const std::vector<U> &phi,
                                const std::vector<U> &theta,
                                const arima_kind &kind,
                                const U kappa = 1000000,
                                const SSinit state_init = Gardner,
                                const U tol = 1e-9) {
  std::vector<U> delta = make_delta<U>(kind.d(), kind.period(), kind.D());
  const size_t p = kind.p() + (kind.P() * kind.period());
  const size_t q = kind.q() + (kind.Q() * kind.period());
  const size_t r = max(p, q + 1), d = delta.size(), rd = r + d;

  size_t i, j;
  std::vector<U> Z(rd);
  Z[0] = 1;
  for (i = 1; i < r - 1; i++) {
    Z[i] = 0;
  }
  j = 0;
  for (i = r; i < rd; i++) {
    Z[i] = delta[j];
    j++;
  }
  std::vector<U> T(rd * rd);
  if (p > 0) {
    for (i = 0; i < p; i++) {
      T[i] = phi[i];
    }
  }
  if (r > 1L) {
    /* set '2nd diagonal' elements to 1. since this is hard to understand,
     * here are some examples
     * for a matrix with *rd* == 5, transform it:
     *
     *   (input)          (rd == r)         (rd > r)
     *   x x x x x   =>   x 1 x x x  //  => x 1 x x x
     *   x x x x x   =>   x x 1 x x  //  => x x 1 x x
     *   x x x x x   =>   x x x 1 x  //  => x x x 1 x
     *   x x x x x   =>   x x x x 1  //  => x x x x x
     *   x x x x x   =>   x x x x x  //  => x x x x x
     */
    for (i = 0; i < r - 1; i++) {
      T[(rd * (i + 1)) + i] = 1;
    }
  }
  if (d > 0) {
    // replace row r+1 in R (or row r for us)
    // with whatever is in Z
    for (j = 0; j < rd; j++) {
      T[(j * rd) + r] = Z[j];
    }
    // if there are more than 1 differences d
    if (d > 1) {
      /* start replacing at r and continue until you get to r + d (-1?)
       * replace r + 2:d with 1 - this is similar as above, but it accounts
       * for the first difference differently(that is taken care of in the above
       * code, in the j < rd loop). here we are taking care of the other
       * differences, so if we have 3 differences, we will only be doing 2
       * replacements. these happen after the for a matrix with *rd* == 5, with
       * *d* == 3, transform it: (input)          (d == 3) x x x x x   =>   x x
       * x x x x x x x x   =>   x x x x x x x x x x   =>   x x x x x x x x x x
       * =>   x x 1 x x x x x x x   =>   x x x 1 x
       */
      for (i = r; i < rd - 1; i++) {
        T[((rd + 1) * i) + 1] = 1;
      }
    }
  }
  // this is R <- c(1, theta, rep.int(0, d))
  // we can skip the d part as vectors are 0 initialized.
  std::vector<U> R(1 + theta.size() + d);
  R[0] = 1;
  for (i = 1; i < theta.size() + 1; i++) {
    R[i] = theta[i - 1];
  }
  std::vector<U> V(R.size() * R.size());
  // here we do an outer product, ie: V <- R %o% R
  size_t mat_p = 0;
  for (i = 0; i < R.size(); i++) {
    for (j = 0; j < R.size(); j++) {
      V[mat_p] = R[i] * R[j];
      mat_p++;
    }
  }
  U h = 0;
  std::vector<U> a(rd);
  std::vector<U> P(rd * rd);
  std::vector<U> Pn(rd * rd);
  if (r > 1) {
    // for storing initialization results
    std::vector<U> temp(r * r);
    switch (state_init) {
    case Gardner:
      temp = std::move(get_Q0(phi, theta));
      break;
    case Rossignol:
      temp = std::move(get_Q0_rossignol(phi, theta));
      break;
    };
    /* update a block of first r rows and columns i.e. if we have a 5x5 Pn
     * matrix, and r == 3, then we update the highlighted parts: (input)
     * (updated) x x x x x   =>    y y y|x x x x x x x   =>    y y y|x x x x x x
     * x   =>    y y y|x x
     *                     _____
     *   x x x x x   =>    x x x x x
     *   x x x x x   =>    x x x x x
     */
    mat_p = 0;
    for (j = 0; j < r; j++) {
      for (i = 0; i < r; i++) {
        Pn[(j * rd) + i] = std::move(temp[mat_p]);
        mat_p++;
      }
    }
  } else {
    Pn[0] = (p > 0) * (1 / (1 - pow(phi[0], 2))) + (p == 0);
  }
  if (d > 0L) {
    /* update diagonal elements which come after the coefficients -
     * diagonal entries between r and rd - with kappa
     */
    for (i = r; i < rd; i++) {
      for (j = r; j < rd; j++) {
        // to only update diagonal elements check that we are on the diagonal
        // otherwise we have a zero - as intended
        Pn[(j * rd) + i] = (i == j) * kappa;
      }
    }
  }
  structural_model<U> res(phi, theta, delta, Z, a, P, T, V, h, Pn);
  return res;
}


template <typename U = double>
void update_arima(structural_model<U> &model, std::vector<U> &phi,
                  std::vector<U> &theta, SSinit state_init = Gardner) {
  const size_t p = phi.size(), q = theta.size(), r = max(p, q + 1), rd = model.Z.size();

  model.phi = phi;
  model.theta = theta;

  if (p > 0) {
    for (size_t i = 0; i < p; i++) {
      model.T[i] = phi[i];
    }
  }

  if (r > 1) {
    // for storing initialization results
    std::vector<U> temp(r * r);
    switch (state_init) {
    case Gardner:
      temp = std::move(get_Q0(phi, theta));
      break;
    case Rossignol:
      temp = std::move(get_Q0_rossignol(phi, theta));
      break;
    };
    size_t mat_p = 0;
    /* update a block of first r rows and columns i.e. if we have a 5x5 Pn
     * matrix, and r == 3, then we update the highlighted parts:
     *   (input)           (updated)
     *   x x x x x   =>    y y y|x x
     *   x x x x x   =>    y y y|x x
     *   x x x x x   =>    y y y|x x
     *                     _____
     *   x x x x x   =>    x x x x x
     *   x x x x x   =>    x x x x x
     */
    for (size_t j = 0; j < r; j++) {
      for (size_t i = 0; i < r; i++) {
        model.Pn[(j * rd) + i] = std::move(temp[mat_p]);
        mat_p++;
      }
    }
  } else {
    model.Pn[0] = (p > 0) * (1 / (1 - pow(phi[0], 2))) + (p == 0);
  }
  // set a to all zero:
  std::fill(model.a.begin(), model.a.end(), 0);
}

template <class T, typename U = double>
void update_arima(structural_model<U> &model,
                  T &coef,
                  const arima_kind kind,
                  SSinit state_init = Gardner) {
  const size_t p = kind.p() + kind.period() * kind.P(),
    q = kind.q() + kind.period() * kind.Q(),
    r = max(p, q + 1), rd = model.Z.size();

  // copy out elements of coef into phi and theta
  for( size_t i = 0; i < p; i++) {
    model.phi[i] = coef[i];
  }
  for( size_t i = p; i < p + q; i++) {
    model.theta[i-p] = coef[i];
  }

  if (p > 0) {
    for (size_t i = 0; i < p; i++) {
      model.T[i] = model.phi[i];
    }
  }

  if (r > 1) {
    // for storing initialization results
    std::vector<U> temp(r * r);
    switch (state_init) {
    case Gardner:
      temp = std::move(get_Q0(model.phi, model.theta));
      break;
    case Rossignol:
      temp = std::move(get_Q0_rossignol(model.phi, model.theta));
      break;
    };
    size_t mat_p = 0;
    /* update a block of first r rows and columns i.e. if we have a 5x5 Pn
     * matrix, and r == 3, then we update the highlighted parts:
     *   (input)           (updated)
     *   x x x x x   =>    y y y|x x
     *   x x x x x   =>    y y y|x x
     *   x x x x x   =>    y y y|x x
     *                     _____
     *   x x x x x   =>    x x x x x
     *   x x x x x   =>    x x x x x
     */
    for (size_t j = 0; j < r; j++) {
      for (size_t i = 0; i < r; i++) {
        model.Pn[(j * rd) + i] = std::move(temp[mat_p]);
        mat_p++;
      }
    }
  } else {
    model.Pn[0] = (p > 0) * (1 / (1 - pow(model.phi[0], 2))) + (p == 0);
  }
  // set a to all zero:
  std::fill(model.a.begin(), model.a.end(), 0);
}

template <typename U = double> void recompute_v(structural_model<U> & model) {
  // reexpand V for an already existing structural model

  size_t i, j;
  std::vector<U> R(1 + model.theta.size() + model.delta.size());
  R[0] = 1;
  for (i = 1; i < model.theta.size() + 1; i++) {
    R[i] = model.theta[i - 1];
  }
  std::vector<U> V(R.size() * R.size());
  // here we do an outer product, ie: V <- R %o% R
  size_t mat_p = 0;
  for (i = 0; i < R.size(); i++) {
    for (j = 0; j < R.size(); j++) {
      V[mat_p] = R[i] * R[j];
      mat_p++;
    }
  }
  model.V = V;
}

#endif
