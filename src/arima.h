#ifndef ARIMA_WRAPPER
#define ARIMA_WRAPPER

#include "initializers.h"
#include "utils.h"

/* Arima class and other structures for modelling.
 */

template <typename U=double> class Arima {
  Arima<U>(){};
  void fit(){};
  void predict(){};



  std::vector<U> y;
  std::vector<int> order;
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

enum SSinit {
  Gardner = 1,
  Rossignol = 2
};


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

#endif
