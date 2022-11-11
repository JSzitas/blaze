#ifndef ARIMA_WRAPPER
#define ARIMA_WRAPPER

#include "initializers.h"

/* Arima class and other structures for modelling.
 */

// Kalman filtering model structure
template <typename U=double> struct structural_model {
//   structural_model() {
//
//   }
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
                                                             U kappa = 100000,
                                                             SSinit state_init = Gardner,
                                                             U tol = 1e-9 ) {
  int p = phi.size(), q = theta.size(), r = max(p, q+1), d = delta.size(), rd = r + d;
        // p <- length(phi)
        // q <- length(theta)
        // r <- max(p, q + 1L)
        // d <- length(Delta)
        // rd <- r + d
  int i, j;
  std::vector<U> Z(rd);
  Z[0] = 1;
  for(i = 1; i < r-1; i++) {
    Z[i] = 0;
  }
  j = 0;
  for(; i < rd; i++) {
    Z[i] = delta[j];
    j++;
  }
  std::vector<U> T(rd*rd);
  if (p > 0) {
    T[1L:p, 1L] <- phi
  }
  if (r > 1L) {
    ind <- 2:r
    T[cbind(ind - 1L, ind)] <- 1
  }
  if (d > 0L) {
    T[r + 1L, ] <- Z
    if (d > 1L) {
      ind <- r + 2:d
      T[cbind(ind, ind - 1)] <- 1
    }
  }
  if (q < r - 1L) {
    theta <- c(theta, rep.int(0, r - 1L - q))
  }
  R <- c(1, theta, rep.int(0, d))
  V <- R %o% R
  h <- 0
  a <- rep(0, rd)
  Pn <- P <- matrix(0, rd, rd)
  if (r > 1L) {
    Pn[1L:r, 1L:r] <- switch(match.arg(SSinit), Gardner1980 = .Call(C_getQ0,
                                       phi, theta), Rossignol2011 = .Call(C_getQ0bis, phi,
                                       theta, tol), stop("invalid 'SSinit'"))
  } else {
    Pn[1L, 1L] <- if (p > 0)
      1/(1 - phi^2)
      else 1
  }
  if (d > 0L) {
    Pn[cbind(r + 1L:d, r + 1L:d)] <- kappa
  }
  list(phi = phi, theta = theta, Delta = Delta, Z = Z, a = a,
       P = P, T = T, V = V, h = h, Pn = Pn)
}























#endif
