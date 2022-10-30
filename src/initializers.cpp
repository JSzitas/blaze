#include <Rcpp.h>
using namespace Rcpp;

// /*
//  Matwey V. Kornilov's implementation of algorithm by
//  Dr. Raphael Rossignol
//  See https://bugs.r-project.org/show_bug.cgi?id=14682 for details.
//  */
SEXP getQ0bis(SEXP sPhi, SEXP sTheta, SEXP sTol)
{
  SEXP res;
  int p = LENGTH(sPhi), q = LENGTH(sTheta);
  double *phi = REAL(sPhi), *theta = REAL(sTheta); // tol = REAL(sTol)[0];

  int i,j, r = max(p, q + 1);

  /* Final result is block product
   *   Q0 = A1 SX A1^T + A1 SXZ A2^T + (A1 SXZ A2^T)^T + A2 A2^T ,
   * where A1 [i,j] = phi[i+j],
   *       A2 [i,j] = ttheta[i+j],  and SX, SXZ are defined below */
  PROTECT(res = allocMatrix(REALSXP, r, r));
  double *P = REAL(res);

  /* Clean P */
  Memzero(P, r*r);

#define _ttheta(j) ttheta[j]
#define _tphi(j) tphi[j]
#define _rrz(j)  rrz [j]

  double *ttheta = (double *) R_alloc(q + 1, sizeof(double));
  /* Init ttheta = c(1, theta) */
  ttheta[0] = 1.;
  for (i = 1; i < q + 1; ++i) ttheta[i] = theta[i - 1];

  if( p > 0 ) {
    int r2 = max(p + q, p + 1);
    SEXP sgam = PROTECT(allocMatrix(REALSXP, r2, r2)),
      sg = PROTECT(allocVector(REALSXP, r2));
    double *gam = REAL(sgam);
    double *g = REAL(sg);
    double *tphi = (double *) R_alloc(p + 1, sizeof(double));
    /* Init tphi = c(1, -phi) */
    tphi[0] = 1.;
    for (i = 1; i < p + 1; ++i) tphi[i] = -phi[i - 1];

    /* Compute the autocovariance function of U, the AR part of X */

    /* Gam := C1 + C2 ; initialize */
    Memzero(gam, r2*r2);

    /* C1[E] */
    for (j = 0; j < r2; ++j)
      for (i = j; i < r2 && i - j < p + 1; ++i)
        gam[j*r2 + i] += _tphi(i-j);

    /* C2[E] */
    for (i = 0; i < r2; ++i)
      for (j = 1; j < r2 && i + j < p + 1; ++j)
        gam[j*r2 + i] += _tphi(i+j);

    /* Initialize g = (1 0 0 .... 0) */
    g[0] = 1.;
    for (i = 1; i < r2; ++i)
      g[i] = 0.;

    /* rU = solve(Gam, g)  -> solve.default() -> .Internal(La_solve, .,)
     * --> fiddling with R-objects -> C and then F77_CALL(.) of dgesv, dlange, dgecon
     * FIXME: call these directly here, possibly even use 'info' instead of error(.)
     * e.g., in case of exact singularity.
     */
    SEXP callS = PROTECT(lang4(install("solve.default"), sgam, sg, sTol)),
      su = PROTECT(eval(callS, R_BaseEnv));
    double *u = REAL(su);
    /* SX = A SU A^T */
    /* A[i,j]  = ttheta[j-i] */
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
                  _ttheta(L - k) * _ttheta(n - m) * u[abs(L - n)];
              }
            }
          }
        }
      }
    }

    UNPROTECT(4);
    /* Compute correlation matrix between X and Z */
    /* forwardsolve(C1, g) */
    /* C[i,j] = tphi[i-j] */
    /* g[i] = _ttheta(i) */
    double *rrz = (double *) R_alloc(q, sizeof(double));
    if(q > 0) {
      for (i = 0; i < q; ++i) {
        rrz[i] = _ttheta(i);
        for (j = max(0, i - p); j < i; ++j)
          rrz[i] -= _rrz(j) * _tphi(i-j);
      }
    }

    /* Q0 += A1 SXZ A2^T + (A1 SXZ A2^T)^T */
    /* SXZ[i,j] = rrz[j-i-1], j > 0 */
    for (i = 0; i < r; ++i)
      for (j = i; j < r; ++j) {
        int k, L;
        for (k = 0; i + k < p; ++k)
          for (L = k+1; j + L < q + 1; ++L)
            P[r*i + j] += phi[i + k] * _ttheta(j + L) * _rrz(L - k - 1);
        for (k = 0; j + k < p; ++k)
          for (L = k+1; i + L < q + 1; ++L)
            P[r*i + j] += phi[j + k] * _ttheta(i + L) * _rrz(L - k - 1);
      }
  } // end if(p > 0)

  /* Q0 += A2 A2^T */
  for (i = 0; i < r; ++i)
    for (j = i; j < r; ++j)
      for (int k = 0; j + k < q + 1; ++k)
        P[r*i + j] += _ttheta(i + k) * _ttheta(j + k);

  /* Symmetrize result */
  for (i = 0; i < r; ++i)
    for (j = i+1; j < r; ++j)
      P[r*j + i] = P[r*i + j];

  UNPROTECT(1);
  return res;
}

// /* based on code from AS154 */
static void  inclu2(size_t np, double *xnext, double *xrow, double ynext,
                    double *d, double *rbar, double *thetab)
{
    double cbar, sbar, di, xi, xk, rbthis, dpi;
    size_t i, k, ithisr;

    /*   This subroutine updates d, rbar, thetab by the inclusion
     of xnext and ynext. */

    for (i = 0; i < np; i++) xrow[i] = xnext[i];

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
        if (di == 0.0) return;
      } else
        ithisr = ithisr + np - i - 1;
    }
  }

SEXP getQ0(SEXP sPhi, SEXP sTheta)
{
  SEXP res;
  int  p = LENGTH(sPhi), q = LENGTH(sTheta);
  double *phi = REAL(sPhi), *theta = REAL(sTheta);

  /* thetab[np], xnext[np], xrow[np].  rbar[rbar] */
  /* NB: nrbar could overflow */
  int r = max(p, q + 1);
  size_t np = r * (r + 1) / 2, nrbar = np * (np - 1) / 2, npr, npr1;
  size_t indi, indj, indn, i, j, ithisr, ind, ind1, ind2, im, jm;


  /* This is the limit using an int index.  We could use
   size_t and get more on a 64-bit system,
   but there seems no practical need. */
  // if(r > 350) error(_("maximum supported lag is 350"));
  double *xnext, *xrow, *rbar, *thetab, *V;
  xnext = (double *) R_alloc(np, sizeof(double));
  xrow = (double *) R_alloc(np, sizeof(double));
  rbar = (double *) R_alloc(nrbar, sizeof(double));
  thetab = (double *) R_alloc(np, sizeof(double));
  V = (double *) R_alloc(np, sizeof(double));
  for (ind = 0, j = 0; j < r; j++) {
    double vj = 0.0;
    if (j == 0) vj = 1.0; else if (j - 1 < q) vj = theta[j - 1];
    for (i = j; i < r; i++) {
      double vi = 0.0;
      if (i == 0) vi = 1.0; else if (i - 1 < q) vi = theta[i - 1];
      V[ind++] = vi * vj;
    }
  }

  PROTECT(res = allocMatrix(REALSXP, r, r));
  double *P = REAL(res);

  if (r == 1) {
    if (p == 0) P[0] = 1.0; // PR#16419
    else P[0] = 1.0 / (1.0 - phi[0] * phi[0]);
    UNPROTECT(1);
    return res;
  }
  if (p > 0) {
    /*      The set of equations s * vec(P0) = vec(v) is solved for
     vec(P0).  s is generated row by row in the array xnext.  The
     order of elements in P is changed, so as to bring more leading
     zeros into the rows of s. */

    for (i = 0; i < nrbar; i++) rbar[i] = 0.0;
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
      double phij = (j < p) ? phi[j] : 0.0;
      xnext[indj++] = 0.0;
      indi = npr1 + j;
      for (i = j; i < r; i++) {
        double ynext = V[ind++];
        double phii = (i < p) ? phi[i] : 0.0;
        if (j != r - 1) {
          xnext[indj] = -phii;
          if (i != r - 1) {
            xnext[indi] -= phij;
            xnext[++ind1] = -1.0;
          }
        }
        xnext[npr] = -phii * phij;
        if (++ind2 >= np) ind2 = 0;
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
      double bi = thetab[im];
      for (jm = np - 1, j = 0; j < i; j++)
        bi -= rbar[ithisr--] * P[jm--];
      P[im--] = bi;
    }

    /*        now re-order p. */

    ind = npr;
    for (i = 0; i < r; i++) xnext[i] = P[ind++];
    ind = np - 1;
    ind1 = npr - 1;
    for (i = 0; i < npr; i++) P[ind--] = P[ind1--];
    for (i = 0; i < r; i++) P[i] = xnext[i];
  } else {

    /* P0 is obtained by backsubstitution for a moving average process. */

    indn = np;
    ind = np;
    for (i = 0; i < r; i++)
      for (j = 0; j <= i; j++) {
        --ind;
        P[ind] = V[ind];
        if (j != 0) P[ind] += P[--indn];
      }
  }
  /* now unpack to a full matrix */
  for (i = r - 1, ind = np; i > 0; i--)
    for (j = r - 1; j >= i; j--)
      P[r * i + j] = P[--ind];
  for (i = 0; i < r - 1; i++)
    for (j = i + 1; j < r; j++)
      P[i + r * j] = P[j + r * i];
  UNPROTECT(1);
  return res;
}

