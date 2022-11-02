#include <Rcpp.h>
using namespace Rcpp;
#include "utils.h"
#include <Eigen/Dense>


// map 2 vectors to Eigen matrices and call solve
// std::vector<double> solve_mat_vec( std::vector<double> &mat,
//                                    std::vector<double> &vec ) {
//   const int n = vec.size();
//   double *mat_ptr = &mat[n*n];
//   double *vec_ptr = &vec_ptr[n];
//   Eigen::Map<Eigen::MatrixXd> new_mat(map_ptr, n, n);
//   Eigen::Map<Eigen::VectorXd> new_vec(vec);
// }


/*
 * Reimplementation of the Matwey V. Kornilov's implementation of algorithm by
 * Dr. Raphael Rossignol, See https://bugs.r-project.org/show_bug.cgi?id=14682
 * for details (of algorithm).
 * Avoids R specific data structures (i.e. SEXPs and related types) in favor
 * of standard C++
 */
std::vector<double> get_Q0_rossignol(std::vector<double> & phi_coef,
                                     std::vector<double> & theta_coef,
                                     double tol)
{
  const int p = phi_coef.size(), q = theta_coef.size();
  int i,j, r = max(p, q + 1);
  // in the original, you create a pointer to the result, which is this,
  // and is a matrix (P),
  // but it makes more sense to me to model it as a std::vector -
  // since working with R matrices is not better than that.
  std::vector<double> P(r*r);
  /* Final result is block product
   *   Q0 = A1 SX A1^T + A1 SXZ A2^T + (A1 SXZ A2^T)^T + A2 A2^T ,
   * where A1 [i,j] = phi[i+j],
   *       A2 [i,j] = ttheta[i+j],  and SX, SXZ are defined below */

  // Initialize and allocate theta
  std::vector<double> theta( q+1 );
  theta[0] = 1.0;
  for( i = 1; i < q+1; i++ ) {
    theta[i] = theta_coef[i-1];
  }

  if( p > 0 ) {
    int r2 = max(p + q, p + 1);
    // initialize phi
    std::vector<double> phi(p+1);
    // fill with correct values (i.e. a leading 1)
    phi[0] = 1.0;
    for (i = 1; i < p + 1; i++) {
      phi[i] = -phi_coef[i - 1];
    }
    /* Compute the autocovariance function of U, the AR part of X */
    /* Gam := C1 + C2 ; initialize */
    // gam is a matrix
    std::vector<double> gam( r2 * r2 );
    std::vector<double> g(r2);
    /* C1[E] */
    for (j = 0; j < r2; ++j) {
      for (i = j; i < r2 && i - j < p + 1; ++i) {
        gam[j*r2 + i] += phi[i-j];
      }
    }
    /* C2[E] */
    for (i = 0; i < r2; ++i) {
      for (j = 1; j < r2 && i + j < p + 1; ++j) {
        gam[j*r2 + i] += phi[i+j];
      }
    }
    /* Initialize g = (1 0 0 .... 0) */
    // This can be replaced in modern C++ by just calling the constructor and writing
    // 1 to the first position
    g[0] = 1.;
    for (i = 1; i < r2; ++i) {
      g[i] = 0.;
    }
    /* rU = solve(Gam, g)  -> solve.default() -> .Internal(La_solve, .,)
     * --> fiddling with R-objects -> C and then F77_CALL(.) of dgesv, dlange, dgecon
     * FIXME: call these directly here, possibly even use 'info' instead of error(.)
     * e.g., in case of exact singularity.
     * ^ The above was a call for help in the original code to replace a
     * convoluted call to solve.default with something saner. I will try
     * to comply using Eigen, and calling Eigen's solvers directly.
     */

    SEXP callS = PROTECT(lang4(install("solve.default"), gam, g, tol)),
    su = PROTECT(eval(callS, R_BaseEnv));
    double *u = REAL(su);
    /* SX = A SU A^T */
    /* A[i,j]  = theta[j-i] */
    /* SU[i,j] = u[abs(i-j)] */
    /* Q0 += ( A1 SX A1^T == A1 A SU A^T A1^T) */
    // (relying on good compiler optimization here:)
    for (i = 0; i < r; ++i) {
      for (j = i; j < r; ++j) {
        for (int k = 0; i + k < p; ++k) {
          for (int L = k; L - k < q + 1; ++L) {
            for (int m = 0; j + m < p; ++m) {
              for (int n = m; n - m < q + 1; ++n) {
                P[r*i + j] += phi[i + k] * phi[j + m] *
                  theta[L - k] * theta[n - m] * u[abs(L - n)];
              }
            }
          }
        }
      }
    }
    /* Compute correlation matrix between X and Z */
    /* forwardsolve(C1, g) */
    /* C[i,j] = tphi[i-j] */
    /* g[i] = _ttheta(i) */
    std::vector<double> rrz(q);
    if(q > 0) {
      for (i = 0; i < q; ++i) {
        rrz[i] = theta[i];
        for (j = max(0, i - p); j < i; ++j) {
          rrz[i] -= rrz[j] * phi[i-j];
        }
      }
    }
    int k, L;
    /* Q0 += A1 SXZ A2^T + (A1 SXZ A2^T)^T */
    /* SXZ[i,j] = rrz[j-i-1], j > 0 */
    for (i = 0; i < r; ++i)
      for (j = i; j < r; ++j) {
        for (k = 0; i + k < p; ++k) {
          for (L = k+1; j + L < q + 1; ++L) {
            P[r*i + j] += phi[i + k] * theta[j + L] * rrz[L - k - 1];
          }
        }
        for (k = 0; j + k < p; ++k) {
          for (L = k+1; i + L < q + 1; ++L) {
            P[r*i + j] += phi[j + k] * theta[i + L] * rrz[L - k - 1];
          }
        }
      }
  } // end if(p > 0)

  /* Q0 += A2 A2^T */
  for (i = 0; i < r; ++i) {
    for (j = i; j < r; ++j) {
      for (int k = 0; j + k < q + 1; ++k) {
        P[r*i + j] += theta[i + k] * theta[j + k];
      }
    }
  }
  /* Symmetrize result */
  for (i = 0; i < r; ++i) {
    for (j = i+1; j < r; ++j) {
      P[r*j + i] = P[r*i + j];
    }
  }
  return P;
}
