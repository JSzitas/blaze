#ifndef STL
#define STL

#include "third_party/eigen.h"

template <typename T, typename scalar_t=double> T div_vec(
  const T &x, const scalar_t divisor) {
  T result(x.size());
  size_t p = 0;
  for( auto &val: x ) {
    result[p++] = val/divisor;
  }
  return result;
}

template <typename T> void tricube(T& x) {
  for( auto &val:x ) {
    val = pow(1 - pow(val,3), 3);
  }
}

template <typename T, typename scalar_t=double> T sqr_dist(
  const T & x, const scalar_t value) {
  T result(x.size());
  size_t p = 0;
  for( auto &val:x ) {
    result[p] = make_pair(pow(val - value,2), p);
    p++;
  }
  return result;
}


template <typename scalar_t = double> std::vector<scalar_t> loess(
  const std::vector<scalar_t> & y,
  const scalar_t alpha = 0.3,
  const size_t degree = 2) {

  using EigMat = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;
  using EigVec = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;

  const auto span = (size_t)(alpha * y.size());
  std::vector<scalar_t> result(y.size());
  std::vector<std::pair<scalar_t, size_t>> distances(y.size());
  EigVec W(span), nearest(span), Y(span);
  // preallocate projection matrix
  EigMat P(span, degree + 1);

  for( size_t i = 0; i < y.size(); i++ ) {
    distances = sqr_dist(y, y[i]);
    std::sort(distances.begin(), distances.end());
    for( size_t j = 0; j < span; j++ ) nearest[i] = distances[i];
    // worst point
    scalar_t worst = std::get<0>(distances[y.size()-1]);
    // update weights
    for( size_t i = 0; i < span; i++ ){
      W[i] = std::get<0>(distances[i]);
    }
    tricube(W);
    // create polynomial
    for(size_t d = 0; d < degree+1; d++) {
      P.col(d) = pow(nearest, d);
    }
    auto V = P.t() * W * P;
    auto Y =



    // W = np.diag(w)
    // A = np.vander(Nx, N=1+deg)

    V = np.matmul(np.matmul(A.T, W), A)
    Y = np.matmul(np.matmul(A.T, W), Ny)
    Q, R = qr(V)
    p = solve_triangular(R, np.matmul(Q.T, Y))
    result[i] = np.polyval(p, val)
  }
  return result;
}

// def loess(X, y, alpha, deg, all_x = True, num_points = 100):
//   '''
// Parameters
//   ----------
//     X : numpy array 1D
//   Explanatory variable.
// y : numpy array 1D
//   Response varible.
// alpha : double
//   Proportion of the samples to include in local regression.
// deg : int
//   Degree of the polynomial to fit. Option 1 or 2 only.
// all_x : boolean, optional
//   Include all x points as target. The default is True.
// num_points : int, optional
//   Number of points to include if all_x is false. The default is 100.
//
//   Returns
//     -------
//       y_hat : numpy array 1D
//     Y estimations at each focal point.
//   x_space : numpy array 1D
//     X range used to calculate each estimation of y.
//
//   '''

  // assert (deg == 1) or (deg == 2), "Deg has to be 1 or 2"
  // assert (alpha > 0) and (alpha <=1), "Alpha has to be between 0 and 1"
  // assert len(X) == len(y), "Length of X and y are different"



      n = len(X)
      span = int(np.ceil(alpha * n))


      y_hat = np.zeros(len(X_domain))
      x_space = np.zeros_like(X_domain)

      for i, val in enumerate(X_domain):

        distance = abs(X - val)
        sorted_dist = np.sort(distance)
        ind = np.argsort(distance)

        Nx = X[ind[:span]]
      Ny = y[ind[:span]]

      delx0 = sorted_dist[span-1]

      u = distance[ind[:span]] / delx0
        w = (1 - u**3)**3

      W = np.diag(w)
        A = np.vander(Nx, N=1+deg)

        V = np.matmul(np.matmul(A.T, W), A)
        Y = np.matmul(np.matmul(A.T, W), Ny)
        Q, R = qr(V)
        p = solve_triangular(R, np.matmul(Q.T, Y))
        y_hat[i] = np.polyval(p, val)
        x_space[i] = val

        return y_hat, x_space



#endif
