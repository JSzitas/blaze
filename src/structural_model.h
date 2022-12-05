#ifndef STRUCT_MODEL
#define STRUCT_MODEL

#include "vector"
#include "utils/utils.h"

enum SSinit {
  Gardner = 1,
  Rossignol = 2
};

// Kalman filtering model structure
template <typename U=double> struct structural_model {
  structural_model<U>( std::vector<U> phi,
                            std::vector<U> theta,
                            std::vector<U> delta,
                            std::vector<U> Z,
                            std::vector<U> a,
                            std::vector<U> P,
                            std::vector<U> T,
                            std::vector<U> V,
                            U h,
                            std::vector<U> Pn) :
  phi(std::move(phi)), theta(std::move(theta)), delta(std::move(delta)),
  Z(std::move(Z)), a(std::move(a)), P(std::move(P)),
  T(std::move(T)), V(std::move(V)), h(std::move(h)), Pn(std::move(Pn)){};

  structural_model<U>( std::vector<U>&& phi,
                            std::vector<U>&& theta,
                            std::vector<U>&& delta,
                            std::vector<U>&& Z,
                            std::vector<U>&& a,
                            std::vector<U>&& P,
                            std::vector<U>&& T,
                            std::vector<U>&& V,
                            U&& h,
                            std::vector<U>&& Pn) :
    phi(std::move(phi)), theta(std::move(theta)), delta(std::move(delta)),
    Z(std::move(Z)), a(std::move(a)), P(std::move(P)),
    T(std::move(T)), V(std::move(V)), h(std::move(h)), Pn(std::move(Pn)){};
  // private:
  std::vector<U> phi;
  std::vector<U> theta;
  std::vector<U> delta;
  std::vector<U> Z;
  std::vector<U> a;
  std::vector<U> P;
  std::vector<U> T;
  std::vector<U> V;
  U h;
  std::vector<U> Pn;
};

template <typename U=double>struct forecast_result {
  forecast_result<U>(int h) {
    this->forecast = std::vector<U>(h);
    this->se = std::vector<U>(h);
  }
  void add(int i, U forecast, U se) {
    this->forecast[i] = forecast;
    this->se[i] = se;
  }
  std::vector<U> forecast;
  std::vector<U> se;
};

/* Forecasts based on state space representation of ARIMA via
 * the kalman filter.
 */
template <typename U=double> forecast_result<U> kalman_forecast( int n_ahead,
                                                                 structural_model<U> &model,
                                                                 bool update = false ){
  int p = model.a.size();
  std::vector<double> anew(p);
  std::vector<double> Pnew(p*p);
  std::vector<double> mm(p*p);

  forecast_result<U> res(n_ahead);

  double fc, tmp;
  for (int l = 0; l < n_ahead; l++) {
    fc = 0.0;
    for (int i = 0; i < p; i++) {
      tmp = 0.0;
      for (int k = 0; k < p; k++) {
        tmp += model.T[i + p * k] * model.a[k];
      }
      anew[i] = tmp;
      fc += tmp * model.Z[i];
    }
    for (int i = 0; i < p; i++) {
      model.a[i] = update * anew[i] - (1-update) * model.a[i];
    }
    for (int i = 0; i < p; i++) {
      for (int j = 0; j < p; j++) {
        tmp = 0.0;
        for (int k = 0; k < p; k++) {
          tmp += model.T[i + p * k] * model.P[k + p * j];
        }
        mm[i + p * j] = tmp;
      }
    }
    for (int i = 0; i < p; i++) {
      for (int j = 0; j < p; j++) {
        tmp = model.V[i + p * j];
        for (int k = 0; k < p; k++) {
          tmp += mm[i + p * k] * model.T[j + p * k];
        }
        Pnew[i + p * j] = tmp;
      }
      tmp = model.h;
    }
    for (int i = 0; i < p; i++) {
      for (int j = 0; j < p; j++) {
        tmp += model.Z[i] * model.Z[j] * Pnew[i + j * p];
        model.P[i + j * p] = update * Pnew[i + j * p] +
          (1-update) * model.P[i + j * p];
      }
      res.add(l, fc, tmp );
    }
  }
  return res;
}

/* originally an R function - this creates the arima model in state space representation
 * from some values of phi, theta and delta
 */
template <typename U=double> structural_model<U> make_arima( std::vector<U> phi,
                                                             std::vector<U> theta,
                                                             std::vector<U> delta,
                                                             U kappa = 1000000,
                                                             SSinit state_init = Gardner,
                                                             U tol = 1e-9 ) {
  int p = phi.size(), q = theta.size(), r = max(p, q+1), d = delta.size(), rd = r + d;
  int i, j;
  std::vector<U> Z(rd);
  Z[0] = 1;
  for(i = 1; i < r-1; i++) {
    Z[i] = 0;
  }
  j = 0;
  for(i = r; i < rd; i++) {
    Z[i] = delta[j];
    j++;
  }
  std::vector<U> T(rd*rd);
  if (p > 0) {
    for(i=0; i < p; i++) {
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
    for( i = 0; i < r-1; i++ ) {
      T[(rd * (i+1)) + i] = 1;
    }
  }
  if (d > 0) {
    // replace row r+1 in R (or row r for us)
    // with whatever is in Z
    for(j = 0; j < rd; j++) {
      T[(j*rd) + r] = Z[j];
    }
    // if there are more than 1 differences d
    if (d > 1) {
      /* start replacing at r and continue until you get to r + d (-1?)
       * replace r + 2:d with 1 - this is similar as above, but it accounts
       * for the first difference differently(that is taken care of in the above code,
       * in the j < rd loop). here we are taking care of the other differences,
       * so if we have 3 differences, we will only be doing 2 replacements.
       * these happen after the
       * for a matrix with *rd* == 5, with *d* == 3, transform it:
       *   (input)          (d == 3)
       *   x x x x x   =>   x x x x x
       *   x x x x x   =>   x x x x x
       *   x x x x x   =>   x x x x x
       *   x x x x x   =>   x x 1 x x
       *   x x x x x   =>   x x x 1 x
       */
      for( i = r; i < rd-1; i++ ) {
        T[((rd+1)*i)+1] = 1;
      }
    }
  }
  if (q < r-1) {
    theta.resize(theta.size() + r-1-q);
  }
  // this is R <- c(1, theta, rep.int(0, d))
  // we can skip the d part as vectors are 0 initialized.
  std::vector<U> R(1+theta.size()+d);
  R[0] = 1;
  for( i=1; i < theta.size()+1; i++ ) {
    R[i] = theta[i-1];
  }
  std::vector<U> V(R.size() * R.size());
  // here we do an outer product, ie: V <- R %o% R
  int mat_p = 0;
  for(i = 0; i < R.size(); i++) {
    for(j=0; j < R.size(); j++) {
      V[mat_p] = R[i] * R[j];
      mat_p++;
    }
  }
  U h = 0;
  std::vector<U> a(rd);
  std::vector<U> P(rd*rd);
  std::vector<U> Pn(rd*rd);
  if (r > 1) {
    // for storing initialization results
    std::vector<U> temp(r*r);
    switch(state_init) {
    case Gardner:
      temp = std::move(get_Q0(phi, theta));
      break;
    case Rossignol:
      temp = std::move(get_Q0_rossignol(phi, theta));
      break;
    };
    mat_p = 0;
    /* update a block of first r rows and columns i.e. if we have a 5x5 Pn matrix,
     * and r == 3, then we update the highlighted parts:
     *   (input)           (updated)
     *   x x x x x   =>    y y y|x x
     *   x x x x x   =>    y y y|x x
     *   x x x x x   =>    y y y|x x
     *                     _____
     *   x x x x x   =>    x x x x x
     *   x x x x x   =>    x x x x x
     */
    for(j=0; j < r; j++) {
      for(i=0; i < r;i++) {
        Pn[(j*rd)+i] = std::move(temp[mat_p]);
        mat_p++;
      }
    }
  } else {
    Pn[0] = (p >0) *(1/(1-pow(phi[0],2))) + (p == 0);
  }
  if (d > 0L) {
    /* update diagonal elements which come after the coefficients -
     * diagonal entries between r and rd - with kappa
     */
    for(i = r; i < rd; i++) {
      for(j=r; j < rd; j++) {
        // to only update diagonal elements check that we are on the diagonal
        // otherwise we have a zero - as intended
        Pn[(j*rd) + i] = (i == j) * kappa;
      }
    }
  }
  structural_model<U> res(phi, theta, delta, Z, a, P, T, V, h, Pn);
  return res;
}


template<typename U=double> void update_arima(structural_model<U> &model,
                                              std::vector<U> &phi,
                                              std::vector<U> &theta,
                                              SSinit state_init = Gardner) {
  int p = phi.size(), q = theta.size(), r = max(p, q+1), rd = model.Z.size();
  model.phi = std::move(phi);
  model.theta = std::move(theta);

  if(p > 0) {
    for(int i=0; i < p; i++) {
      model.T[i] = phi[i];
    }
  }

  if (r > 1) {
    // for storing initialization results
    std::vector<U> temp(r*r);
    switch(state_init) {
    case Gardner:
      temp = std::move(get_Q0(phi, theta));
      break;
    case Rossignol:
      temp = std::move(get_Q0_rossignol(phi, theta));
      break;
    };
    int mat_p = 0;
    /* update a block of first r rows and columns i.e. if we have a 5x5 Pn matrix,
     * and r == 3, then we update the highlighted parts:
     *   (input)           (updated)
     *   x x x x x   =>    y y y|x x
     *   x x x x x   =>    y y y|x x
     *   x x x x x   =>    y y y|x x
     *                     _____
     *   x x x x x   =>    x x x x x
     *   x x x x x   =>    x x x x x
     */
    for(int j=0; j < r; j++) {
      for(int i=0; i < r;i++) {
        model.Pn[(j*rd)+i] = std::move(temp[mat_p]);
        mat_p++;
      }
    }
  } else {
    model.Pn[0] = (p > 0) *(1/(1-pow(phi[0],2))) + (p == 0);
  }
  // set a to all zero:
  std::fill(model.a.begin(), model.a.end(), 0);
}

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


#endif
