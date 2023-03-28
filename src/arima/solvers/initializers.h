#ifndef INITIALIZERS
#define INITIALIZERS

#include "third_party/eigen.h"
#include "utils/utils.h"
#include "arima/structures/ss_init.h"

// map 2 vectors to Eigen matrices and call solve
template <typename scalar_t = float>
std::vector<scalar_t> solve_mat_vec(std::vector<scalar_t> &mat,
                                    std::vector<scalar_t> &vec) {
  const size_t n = vec.size();

  Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic> new_mat =
      Eigen::Map<Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>>(
          mat.data(), n, n);
  Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> new_vec =
      Eigen::Map<Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>>(vec.data(), n, 1);
  Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> res =
      new_mat.completeOrthogonalDecomposition().solve(new_vec);

  std::vector<scalar_t> result(res.data(), res.data() + res.rows() * res.cols());
  return result;
}

/* A mild reimplementation of the Matwey V. Kornilov's implementation of
 * algorithm by Dr. Raphael Rossignol, See
 * https://bugs.r-project.org/show_bug.cgi?id=14682 for details (of algorithm).
 * Avoids R specific data structures (i.e. SEXPs and related types) in favor
 * of standard C++, uses Eigen (rather than a call to the solver used by R),
 * to solve a system of linear equations, where the Eigen solver might be faster
 * (ie completeOrthogonalDecomposition **should** be a better approach than the
 * regular solver)
 */
template <typename T, typename scalar_t> std::vector<scalar_t>
get_Q0_rossignol(
  const T &phi_coef,
  const T &theta_coef) {
  const size_t p = phi_coef.size(), q = theta_coef.size();
  size_t i, j, r = max(p, q + 1);
  // in the original, you create a pointer to the result, which is this,
  // and is a matrix (P),
  // but it makes more sense to me to model it as a std::vector -
  // since working with R matrices is not better than that.
  std::vector<scalar_t> P(r * r);
  /* Final result is block product
   *   Q0 = A1 SX A1^T + A1 SXZ A2^T + (A1 SXZ A2^T)^T + A2 A2^T ,
   * where A1 [i,j] = phi[i+j],
   *       A2 [i,j] = ttheta[i+j],  and SX, SXZ are defined below */

  // Initialize and allocate theta
  std::vector<scalar_t> theta(q + 1);
  theta[0] = 1.0;
  for (i = 1; i < q + 1; i++) {
    theta[i] = theta_coef[i - 1];
  }

  if (p > 0) {
    size_t r2 = max(p + q, p + 1);
    // initialize phi
    std::vector<scalar_t> phi(p + 1);
    // fill with correct values (i.e. a leading 1)
    phi[0] = 1.0;
    for (i = 1; i < p + 1; i++) {
      phi[i] = -phi_coef[i - 1];
    }
    /* Compute the autocovariance function of U, the AR part of X */
    /* Gam := C1 + C2 ; initialize */
    // gam is a matrix
    std::vector<scalar_t> gam(r2 * r2);
    std::vector<scalar_t> g(r2);
    /* C1[E] */
    for (j = 0; j < r2; ++j) {
      for (i = j; i < r2 && i - j < p + 1; ++i) {
        gam[j * r2 + i] += phi[i - j];
      }
    }
    /* C2[E] */
    for (i = 0; i < r2; ++i) {
      for (j = 1; j < r2 && i + j < p + 1; ++j) {
        gam[j * r2 + i] += phi[i + j];
      }
    }
    /* Initialize g = (1 0 0 .... 0) */
    // This can be replaced in modern C++ by just calling the constructor and
    // writing 1 to the first position
    g[0] = 1.;
    for (i = 1; i < r2; ++i) {
      g[i] = 0.;
    }
    // Solve the system of linear equations - using Eigen
    std::vector<scalar_t> u = solve_mat_vec(gam, g);
    /* SX = A SU A^T */
    /* A[i,j]  = theta[j-i] */
    /* SU[i,j] = u[abs(i-j)] */
    /* Q0 += ( A1 SX A1^T == A1 A SU A^T A1^T) */
    // (relying on good compiler optimization here:)
    for (i = 0; i < r; ++i) {
      for (j = i; j < r; ++j) {
        for (size_t k = 0; i + k < p; ++k) {
          for (size_t L = k; L - k < q + 1; ++L) {
            for (size_t m = 0; j + m < p; ++m) {
              for (size_t n = m; n - m < q + 1; ++n) {
                P[r * i + j] += phi[i + k] * phi[j + m] * theta[L - k] *
                                theta[n - m] * u[abs(L - n)];
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
    std::vector<scalar_t> rrz(q);
    if (q > 0) {
      for (i = 0; i < q; ++i) {
        rrz[i] = theta[i];
        for (j = max(0, i - p); j < i; ++j) {
          rrz[i] -= rrz[j] * phi[i - j];
        }
      }
    }
    size_t k, L;
    /* Q0 += A1 SXZ A2^T + (A1 SXZ A2^T)^T */
    /* SXZ[i,j] = rrz[j-i-1], j > 0 */
    for (i = 0; i < r; ++i)
      for (j = i; j < r; ++j) {
        for (k = 0; i + k < p; ++k) {
          for (L = k + 1; j + L < q + 1; ++L) {
            P[r * i + j] += phi[i + k] * theta[j + L] * rrz[L - k - 1];
          }
        }
        for (k = 0; j + k < p; ++k) {
          for (L = k + 1; i + L < q + 1; ++L) {
            P[r * i + j] += phi[j + k] * theta[i + L] * rrz[L - k - 1];
          }
        }
      }
  } // end if(p > 0)

  /* Q0 += A2 A2^T */
  for (i = 0; i < r; ++i) {
    for (j = i; j < r; ++j) {
      for (size_t k = 0; j + k < q + 1; ++k) {
        P[r * i + j] += theta[i + k] * theta[j + k];
      }
    }
  }
  /* Symmetrize result */
  for (i = 0; i < r; ++i) {
    for (j = i + 1; j < r; ++j) {
      P[r * j + i] = P[r * i + j];
    }
  }
  return P;
}

/* based on code from AS154 */
template <typename scalar_t> inline static void inclu2(
    size_t np,
    std::vector<scalar_t> &xnext,
    std::vector<scalar_t> &xrow, scalar_t ynext,
    std::vector<scalar_t> &d,
    std::vector<scalar_t> &rbar,
    std::vector<scalar_t> &thetab) {
  scalar_t cbar, sbar, di, xi, xk, rbthis, dpi;
  size_t i, k, ithisr;
  /*   This subroutine updates d, rbar, thetab by the inclusion
   of xnext and ynext. */
  for (i = 0; i < np; i++) {
    xrow[i] = xnext[i];
  }

  for (ithisr = 0, i = 0; i < np; i++) {
    if (xrow[i] != 0.0) {
      xi = xrow[i];
      di = d[i];
      dpi = di + xi * xi;
      d[i] = dpi;
      cbar = di / dpi;
      sbar = xi / dpi;
      for (k = i + 1; k < np; k++) {
        xk = xrow[k];
        rbthis = rbar[ithisr];
        xrow[k] = xk - xi * rbthis;
        rbar[ithisr++] = cbar * rbthis + sbar * xk;
      }
      xk = ynext;
      ynext = xk - xi * thetab[i];
      thetab[i] = cbar * thetab[i] + sbar * xk;
      if (di == 0.0) {
        return;
      }
    } else {
      ithisr = ithisr + np - i - 1;
    }
  }
}

template <typename T, typename scalar_t>
std::vector<scalar_t> get_Q0_impl(
    const T &phi_coef,
    const T &theta_coef,
    std::vector<scalar_t> &xnext,
    std::vector<scalar_t> &xrow,
    std::vector<scalar_t> &rbar,
    std::vector<scalar_t> &thetab,
    std::vector<scalar_t> &V,
    std::vector<scalar_t> &P) {

  const size_t p = phi_coef.size(), q = theta_coef.size();
  size_t i, j, r = max(p, q + 1);
  size_t np = r * (r + 1) / 2, nrbar = np * (np - 1) / 2, npr, npr1;
  size_t indi, indj, indn, ithisr, ind, ind1, ind2, im, jm;

  scalar_t vj, vi, bi, ynext, phii, phij;
  for (ind = 0, j = 0; j < r; j++) {
    vj = 0.0;
    if (j == 0) {
      vj = 1.0;
    } else if (j - 1 < q) {
      vj = theta_coef[j - 1];
    }
    for (i = j; i < r; i++) {
      vi = 0.0;
      if (i == 0) {
        vi = 1.0;
      } else if (i - 1 < q) {
        vi = theta_coef[i - 1];
      }
      V[ind++] = vi * vj;
    }
  }
  // result vector
  if (r == 1) {
    if (p == 0) {
      P[0] = 1.0;
    } else {
      P[0] = 1.0 / (1.0 - phi_coef[0] * phi_coef[0]);
    }
    return P;
  }
  if (p > 0) {
    /*      The set of equations s * vec(P0) = vec(v) is solved for
     vec(P0).  s is generated row by row in the array xnext.  The
     order of elements in P is changed, so as to bring more leading
     zeros into the rows of s. */
    for (i = 0; i < nrbar; i++) {
      rbar[i] = 0.0;
    }
    for (i = 0; i < np; i++) {
      P[i] = 0.0;
      thetab[i] = 0.0;
      xnext[i] = 0.0;
    }
    ind = 0;
    ind1 = -1;
    npr = np - r;
    npr1 = npr + 1;
    indj = npr;
    ind2 = npr - 1;
    for (j = 0; j < r; j++) {
      phij = (j < p) ? phi_coef[j] : 0.0;
      xnext[indj++] = 0.0;
      indi = npr1 + j;
      for (i = j; i < r; i++) {
        ynext = V[ind++];
        phii = (i < p) ? phi_coef[i] : 0.0;
        if (j != r - 1) {
          xnext[indj] = -phii;
          if (i != r - 1) {
            xnext[indi] -= phij;
            xnext[++ind1] = -1.0;
          }
        }
        xnext[npr] = -phii * phij;
        if (++ind2 >= np) {
          ind2 = 0;
        }
        xnext[ind2] += 1.0;
        inclu2(np, xnext, xrow, ynext, P, rbar, thetab);
        xnext[ind2] = 0.0;
        if (i != r - 1) {
          xnext[indi++] = 0.0;
          xnext[ind1] = 0.0;
        }
      }
    }
    ithisr = nrbar - 1;
    im = np - 1;
    for (i = 0; i < np; i++) {
      bi = thetab[im];
      for (jm = np - 1, j = 0; j < i; j++) {
        bi -= rbar[ithisr--] * P[jm--];
      }
      P[im--] = bi;
    }
    /*  now re-order p. */
    ind = npr;
    for (i = 0; i < r; i++) {
      xnext[i] = P[ind++];
    }
    ind = np - 1;
    ind1 = npr - 1;
    for (i = 0; i < npr; i++) {
      P[ind--] = P[ind1--];
    }
    for (i = 0; i < r; i++) {
      P[i] = xnext[i];
    }
  } else {
    /* P0 is obtained by backsubstitution for a moving average process. */
    indn = np;
    ind = np;
    for (i = 0; i < r; i++) {
      for (j = 0; j <= i; j++) {
        --ind;
        P[ind] = V[ind];
        if (j != 0) {
          P[ind] += P[--indn];
        }
      }
    }
  }
  /* now unpack to a full matrix */
  for (i = r - 1, ind = np; i > 0; i--) {
    for (j = r - 1; j >= i; j--) {
      P[r * i + j] = P[--ind];
    }
  }
  for (i = 0; i < r - 1; i++) {
    for (j = i + 1; j < r; j++) {
      P[i + r * j] = P[j + r * i];
    }
  }
  return P;
}

template <typename T, typename scalar_t>
std::vector<scalar_t> get_Q0(
    const T &phi_coef,
    const T &theta_coef) {
  const size_t p = phi_coef.size(), q = theta_coef.size(),
    r = max(p, q + 1), np = r * (r + 1) / 2, nrbar = np * (np - 1) / 2;

  std::vector<scalar_t> xnext(np), xrow(np), rbar(nrbar),
    thetab(np), V(np), P(r * r);
  return get_Q0_impl(phi_coef, theta_coef, xnext, xrow, rbar, thetab, V, P);
}

template <typename T, typename scalar_t>
std::vector<scalar_t> get_Q0(
    const T &phi_coef,
    const T &theta_coef,
    std::vector<scalar_t> &xnext,
    std::vector<scalar_t> &xrow,
    std::vector<scalar_t> &rbar,
    std::vector<scalar_t> &thetab,
    std::vector<scalar_t> &V,
    std::vector<scalar_t> &P) {
  return get_Q0_impl(phi_coef, theta_coef, xnext, xrow, rbar, thetab, V, P);
}

#endif
