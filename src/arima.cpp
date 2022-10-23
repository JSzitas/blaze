#include <Rcpp.h>
using namespace Rcpp;
#include <cmath>


// A large part of the original C code base is almost valid standalone C (save some
// protect/unprotect calls, SEXPS and other stuff like that), which is close to
// being almost valid C++... so :)


//
//
// SEXP
//   ARIMA_Like(SEXP sy, SEXP mod, SEXP sUP, SEXP giveResid)
//   {
//     SEXP sPhi = getListElement(mod, "phi"),
//       sTheta = getListElement(mod, "theta"),
//       sDelta = getListElement(mod, "Delta"),
//       sa = getListElement(mod, "a"),
//       sP = getListElement(mod, "P"),
//       sPn = getListElement(mod, "Pn");
//
//     if (TYPEOF(sPhi) != REALSXP || TYPEOF(sTheta) != REALSXP ||
//         TYPEOF(sDelta) != REALSXP || TYPEOF(sa) != REALSXP ||
//         TYPEOF(sP) != REALSXP || TYPEOF(sPn) != REALSXP)
//       error(_("invalid argument type"));
//
//     SEXP res, nres, sResid = R_NilValue;
//     int n = LENGTH(sy), rd = LENGTH(sa), p = LENGTH(sPhi),
//       q = LENGTH(sTheta), d = LENGTH(sDelta), r = rd - d;
//     double *y = REAL(sy), *a = REAL(sa), *P = REAL(sP), *Pnew = REAL(sPn);
//     double *phi = REAL(sPhi), *theta = REAL(sTheta), *delta = REAL(sDelta);
//     double sumlog = 0.0, ssq = 0, *anew, *mm = NULL, *M;
//     int nu = 0;
//     Rboolean useResid = asLogical(giveResid);
//     double *rsResid = NULL /* -Wall */;
//
//     anew = (double *) R_alloc(rd, sizeof(double));
//     M = (double *) R_alloc(rd, sizeof(double));
//     if (d > 0) mm = (double *) R_alloc(rd * rd, sizeof(double));
//
//     if (useResid) {
//       PROTECT(sResid = allocVector(REALSXP, n));
//       rsResid = REAL(sResid);
//     }
//
//     for (int l = 0; l < n; l++) {
//       for (int i = 0; i < r; i++) {
//         double tmp = (i < r - 1) ? a[i + 1] : 0.0;
//         if (i < p) tmp += phi[i] * a[0];
//         anew[i] = tmp;
//       }
//       if (d > 0) {
//         for (int i = r + 1; i < rd; i++) anew[i] = a[i - 1];
//         double tmp = a[0];
//         for (int i = 0; i < d; i++) tmp += delta[i] * a[r + i];
//         anew[r] = tmp;
//       }
//       if (l > asInteger(sUP)) {
//         if (d == 0) {
//           for (int i = 0; i < r; i++) {
//             double vi = 0.0;
//             if (i == 0) vi = 1.0; else if (i - 1 < q) vi = theta[i - 1];
//             for (int j = 0; j < r; j++) {
//               double tmp = 0.0;
//               if (j == 0) tmp = vi; else if (j - 1 < q) tmp = vi * theta[j - 1];
//               if (i < p && j < p) tmp += phi[i] * phi[j] * P[0];
//               if (i < r - 1 && j < r - 1) tmp += P[i + 1 + r * (j + 1)];
//               if (i < p && j < r - 1) tmp += phi[i] * P[j + 1];
//               if (j < p && i < r - 1) tmp += phi[j] * P[i + 1];
//               Pnew[i + r * j] = tmp;
//             }
//           }
//         } else {
//           /* mm = TP */
//           for (int i = 0; i < r; i++)
//             for (int j = 0; j < rd; j++) {
//               double tmp = 0.0;
//               if (i < p) tmp += phi[i] * P[rd * j];
//               if (i < r - 1) tmp += P[i + 1 + rd * j];
//               mm[i + rd * j] = tmp;
//             }
//             for (int j = 0; j < rd; j++) {
//               double tmp = P[rd * j];
//               for (int k = 0; k < d; k++)
//                 tmp += delta[k] * P[r + k + rd * j];
//               mm[r + rd * j] = tmp;
//             }
//             for (int i = 1; i < d; i++)
//               for (int j = 0; j < rd; j++)
//                 mm[r + i + rd * j] = P[r + i - 1 + rd * j];
//
//           /* Pnew = mmT' */
//           for (int i = 0; i < r; i++)
//             for (int j = 0; j < rd; j++) {
//               double tmp = 0.0;
//               if (i < p) tmp += phi[i] * mm[j];
//               if (i < r - 1) tmp += mm[rd * (i + 1) + j];
//               Pnew[j + rd * i] = tmp;
//             }
//             for (int j = 0; j < rd; j++) {
//               double tmp = mm[j];
//               for (int k = 0; k < d; k++)
//                 tmp += delta[k] * mm[rd * (r + k) + j];
//               Pnew[rd * r + j] = tmp;
//             }
//             for (int i = 1; i < d; i++)
//               for (int j = 0; j < rd; j++)
//                 Pnew[rd * (r + i) + j] = mm[rd * (r + i - 1) + j];
//           /* Pnew <- Pnew + (1 theta) %o% (1 theta) */
//           for (int i = 0; i <= q; i++) {
//             double vi = (i == 0) ? 1. : theta[i - 1];
//             for (int j = 0; j <= q; j++)
//               Pnew[i + rd * j] += vi * ((j == 0) ? 1. : theta[j - 1]);
//           }
//         }
//       }
//       if (!ISNAN(y[l])) {
//         double resid = y[l] - anew[0];
//         for (int i = 0; i < d; i++)
//           resid -= delta[i] * anew[r + i];
//
//         for (int i = 0; i < rd; i++) {
//           double tmp = Pnew[i];
//           for (int j = 0; j < d; j++)
//             tmp += Pnew[i + (r + j) * rd] * delta[j];
//           M[i] = tmp;
//         }
//
//         double gain = M[0];
//         for (int j = 0; j < d; j++) gain += delta[j] * M[r + j];
//         if(gain < 1e4) {
//           nu++;
//           ssq += resid * resid / gain;
//           sumlog += log(gain);
//         }
//         if (useResid) rsResid[l] = resid / sqrt(gain);
//         for (int i = 0; i < rd; i++)
//           a[i] = anew[i] + M[i] * resid / gain;
//         for (int i = 0; i < rd; i++)
//           for (int j = 0; j < rd; j++)
//             P[i + j * rd] = Pnew[i + j * rd] - M[i] * M[j] / gain;
//       } else {
//         for (int i = 0; i < rd; i++) a[i] = anew[i];
//         for (int i = 0; i < rd * rd; i++) P[i] = Pnew[i];
//         if (useResid) rsResid[l] = NA_REAL;
//       }
//     }
//
//     if (useResid) {
//       PROTECT(res = allocVector(VECSXP, 3));
//       SET_VECTOR_ELT(res, 0, nres = allocVector(REALSXP, 3));
//       REAL(nres)[0] = ssq;
//       REAL(nres)[1] = sumlog;
//       REAL(nres)[2] = (double) nu;
//       SET_VECTOR_ELT(res, 1, sResid);
//       UNPROTECT(2);
//       return res;
//     } else {
//       nres = allocVector(REALSXP, 3);
//       REAL(nres)[0] = ssq;
//       REAL(nres)[1] = sumlog;
//       REAL(nres)[2] = (double) nu;
//       return nres;
//     }
//   }
//
// /* do differencing here */
// /* arma is p, q, sp, sq, ns, d, sd */
// SEXP
//   ARIMA_CSS(SEXP sy, SEXP sarma, SEXP sPhi, SEXP sTheta,
//             SEXP sncond, SEXP giveResid)
//   {
//     SEXP res, sResid = R_NilValue;
//     double ssq = 0.0, *y = REAL(sy), tmp;
//     double *phi = REAL(sPhi), *theta = REAL(sTheta), *w, *resid;
//     int n = LENGTH(sy), *arma = INTEGER(sarma), p = LENGTH(sPhi),
//       q = LENGTH(sTheta), ncond = asInteger(sncond);
//     int ns, nu = 0;
//     Rboolean useResid = asLogical(giveResid);
//
//     w = (double *) R_alloc(n, sizeof(double));
//     for (int l = 0; l < n; l++) w[l] = y[l];
//     for (int i = 0; i < arma[5]; i++)
//       for (int l = n - 1; l > 0; l--) w[l] -= w[l - 1];
//     ns = arma[4];
//     for (int i = 0; i < arma[6]; i++)
//       for (int l = n - 1; l >= ns; l--) w[l] -= w[l - ns];
//
//     PROTECT(sResid = allocVector(REALSXP, n));
//     resid = REAL(sResid);
//     if (useResid) for (int l = 0; l < ncond; l++) resid[l] = 0;
//
//     for (int l = ncond; l < n; l++) {
//       tmp = w[l];
//       for (int j = 0; j < p; j++) tmp -= phi[j] * w[l - j - 1];
//       for (int j = 0; j < min(l - ncond, q); j++)
//         tmp -= theta[j] * resid[l - j - 1];
//       resid[l] = tmp;
//       if (!ISNAN(tmp)) {
//         nu++;
//         ssq += tmp * tmp;
//       }
//     }
//     if (useResid) {
//       PROTECT(res = allocVector(VECSXP, 2));
//       SET_VECTOR_ELT(res, 0, ScalarReal(ssq / (double) (nu)));
//       SET_VECTOR_ELT(res, 1, sResid);
//       UNPROTECT(2);
//       return res;
//     } else {
//       UNPROTECT(1);
//       return ScalarReal(ssq / (double) (nu));
//     }
//   }
//
//
// /*
//  Matwey V. Kornilov's implementation of algorithm by
//  Dr. Raphael Rossignol
//  See https://bugs.r-project.org/show_bug.cgi?id=14682 for details.
//  */
// SEXP getQ0bis(SEXP sPhi, SEXP sTheta, SEXP sTol)
// {
//   SEXP res;
//   int p = LENGTH(sPhi), q = LENGTH(sTheta);
//   double *phi = REAL(sPhi), *theta = REAL(sTheta); // tol = REAL(sTol)[0];
//
//   int i,j, r = max(p, q + 1);
//
//   /* Final result is block product
//    *   Q0 = A1 SX A1^T + A1 SXZ A2^T + (A1 SXZ A2^T)^T + A2 A2^T ,
//    * where A1 [i,j] = phi[i+j],
//    *       A2 [i,j] = ttheta[i+j],  and SX, SXZ are defined below */
//   PROTECT(res = allocMatrix(REALSXP, r, r));
//   double *P = REAL(res);
//
//   /* Clean P */
//   Memzero(P, r*r);
//
// #ifdef DEBUG_Q0bis
// #define _ttheta(j) chk_V(ttheta, "ttheta", j, q+1)// was  r
// #define _tphi(j)   chk_V(tphi,   "tphi",   j, p+1)
// #define _rrz(j)    chk_V(rrz,    "rrz",    j, q)
// #else
// #define _ttheta(j) ttheta[j]
// #define _tphi(j) tphi[j]
// #define _rrz(j)  rrz [j]
// #endif
//
//   double *ttheta = (double *) R_alloc(q + 1, sizeof(double));
//   /* Init ttheta = c(1, theta) */
//   ttheta[0] = 1.;
//   for (i = 1; i < q + 1; ++i) ttheta[i] = theta[i - 1];
//
//   if( p > 0 ) {
//     int r2 = max(p + q, p + 1);
//     SEXP sgam = PROTECT(allocMatrix(REALSXP, r2, r2)),
//       sg = PROTECT(allocVector(REALSXP, r2));
//     double *gam = REAL(sgam);
//     double *g = REAL(sg);
//     double *tphi = (double *) R_alloc(p + 1, sizeof(double));
//     /* Init tphi = c(1, -phi) */
//     tphi[0] = 1.;
//     for (i = 1; i < p + 1; ++i) tphi[i] = -phi[i - 1];
//
//     /* Compute the autocovariance function of U, the AR part of X */
//
//     /* Gam := C1 + C2 ; initialize */
//     Memzero(gam, r2*r2);
//
//     /* C1[E] */
//     for (j = 0; j < r2; ++j)
//       for (i = j; i < r2 && i - j < p + 1; ++i)
//         gam[j*r2 + i] += _tphi(i-j);
//
//     /* C2[E] */
//     for (i = 0; i < r2; ++i)
//       for (j = 1; j < r2 && i + j < p + 1; ++j)
//         gam[j*r2 + i] += _tphi(i+j);
//
//     /* Initialize g = (1 0 0 .... 0) */
//     g[0] = 1.;
//     for (i = 1; i < r2; ++i)
//       g[i] = 0.;
//
//     /* rU = solve(Gam, g)  -> solve.default() -> .Internal(La_solve, .,)
//      * --> fiddling with R-objects -> C and then F77_CALL(.) of dgesv, dlange, dgecon
//      * FIXME: call these directly here, possibly even use 'info' instead of error(.)
//      * e.g., in case of exact singularity.
//      */
//     SEXP callS = PROTECT(lang4(install("solve.default"), sgam, sg, sTol)),
//       su = PROTECT(eval(callS, R_BaseEnv));
//     double *u = REAL(su);
//     /* SX = A SU A^T */
//     /* A[i,j]  = ttheta[j-i] */
//     /* SU[i,j] = u[abs(i-j)] */
//     /* Q0 += ( A1 SX A1^T == A1 A SU A^T A1^T) */
//     // (relying on good compiler optimization here:)
//     for (i = 0; i < r; ++i)
//       for (j = i; j < r; ++j)
//         for (int k = 0; i + k < p; ++k)
//           for (int L = k; L - k < q + 1; ++L)
//             for (int m = 0; j + m < p; ++m)
//               for (int n = m; n - m < q + 1; ++n)
//                 P[r*i + j] += phi[i + k] * phi[j + m] *
//                   _ttheta(L - k) * _ttheta(n - m) * u[abs(L - n)];
//     UNPROTECT(4);
//     /* Compute correlation matrix between X and Z */
//     /* forwardsolve(C1, g) */
//     /* C[i,j] = tphi[i-j] */
//     /* g[i] = _ttheta(i) */
//     double *rrz = (double *) R_alloc(q, sizeof(double));
//     if(q > 0) {
//       for (i = 0; i < q; ++i) {
//         rrz[i] = _ttheta(i);
//         for (j = max(0, i - p); j < i; ++j)
//           rrz[i] -= _rrz(j) * _tphi(i-j);
//       }
//     }
//
//     /* Q0 += A1 SXZ A2^T + (A1 SXZ A2^T)^T */
//     /* SXZ[i,j] = rrz[j-i-1], j > 0 */
//     for (i = 0; i < r; ++i)
//       for (j = i; j < r; ++j) {
//         int k, L;
//         for (k = 0; i + k < p; ++k)
//           for (L = k+1; j + L < q + 1; ++L)
//             P[r*i + j] += phi[i + k] * _ttheta(j + L) * _rrz(L - k - 1);
//         for (k = 0; j + k < p; ++k)
//           for (L = k+1; i + L < q + 1; ++L)
//             P[r*i + j] += phi[j + k] * _ttheta(i + L) * _rrz(L - k - 1);
//       }
//   } // end if(p > 0)
//
//   /* Q0 += A2 A2^T */
//   for (i = 0; i < r; ++i)
//     for (j = i; j < r; ++j)
//       for (int k = 0; j + k < q + 1; ++k)
//         P[r*i + j] += _ttheta(i + k) * _ttheta(j + k);
//
//   /* Symmetrize result */
//   for (i = 0; i < r; ++i)
//     for (j = i+1; j < r; ++j)
//       P[r*j + i] = P[r*i + j];
//
//   UNPROTECT(1);
//   return res;
// }
//
// SEXP getQ0(SEXP sPhi, SEXP sTheta)
// {
//   SEXP res;
//   int  p = LENGTH(sPhi), q = LENGTH(sTheta);
//   double *phi = REAL(sPhi), *theta = REAL(sTheta);
//
//   /* thetab[np], xnext[np], xrow[np].  rbar[rbar] */
//   /* NB: nrbar could overflow */
//   int r = max(p, q + 1);
//   size_t np = r * (r + 1) / 2, nrbar = np * (np - 1) / 2, npr, npr1;
//   size_t indi, indj, indn, i, j, ithisr, ind, ind1, ind2, im, jm;
//
//
//   /* This is the limit using an int index.  We could use
//    size_t and get more on a 64-bit system,
//    but there seems no practical need. */
//   if(r > 350) error(_("maximum supported lag is 350"));
//   double *xnext, *xrow, *rbar, *thetab, *V;
//   xnext = (double *) R_alloc(np, sizeof(double));
//   xrow = (double *) R_alloc(np, sizeof(double));
//   rbar = (double *) R_alloc(nrbar, sizeof(double));
//   thetab = (double *) R_alloc(np, sizeof(double));
//   V = (double *) R_alloc(np, sizeof(double));
//   for (ind = 0, j = 0; j < r; j++) {
//     double vj = 0.0;
//     if (j == 0) vj = 1.0; else if (j - 1 < q) vj = theta[j - 1];
//     for (i = j; i < r; i++) {
//       double vi = 0.0;
//       if (i == 0) vi = 1.0; else if (i - 1 < q) vi = theta[i - 1];
//       V[ind++] = vi * vj;
//     }
//   }
//
//   PROTECT(res = allocMatrix(REALSXP, r, r));
//   double *P = REAL(res);
//
//   if (r == 1) {
//     if (p == 0) P[0] = 1.0; // PR#16419
//     else P[0] = 1.0 / (1.0 - phi[0] * phi[0]);
//     UNPROTECT(1);
//     return res;
//   }
//   if (p > 0) {
//     /*      The set of equations s * vec(P0) = vec(v) is solved for
//      vec(P0).  s is generated row by row in the array xnext.  The
//      order of elements in P is changed, so as to bring more leading
//      zeros into the rows of s. */
//
//     for (i = 0; i < nrbar; i++) rbar[i] = 0.0;
//     for (i = 0; i < np; i++) {
//       P[i] = 0.0;
//       thetab[i] = 0.0;
//       xnext[i] = 0.0;
//     }
//     ind = 0;
//     ind1 = -1;
//     npr = np - r;
//     npr1 = npr + 1;
//     indj = npr;
//     ind2 = npr - 1;
//     for (j = 0; j < r; j++) {
//       double phij = (j < p) ? phi[j] : 0.0;
//       xnext[indj++] = 0.0;
//       indi = npr1 + j;
//       for (i = j; i < r; i++) {
//         double ynext = V[ind++];
//         double phii = (i < p) ? phi[i] : 0.0;
//         if (j != r - 1) {
//           xnext[indj] = -phii;
//           if (i != r - 1) {
//             xnext[indi] -= phij;
//             xnext[++ind1] = -1.0;
//           }
//         }
//         xnext[npr] = -phii * phij;
//         if (++ind2 >= np) ind2 = 0;
//         xnext[ind2] += 1.0;
//         inclu2(np, xnext, xrow, ynext, P, rbar, thetab);
//         xnext[ind2] = 0.0;
//         if (i != r - 1) {
//           xnext[indi++] = 0.0;
//           xnext[ind1] = 0.0;
//         }
//       }
//     }
//
//     ithisr = nrbar - 1;
//     im = np - 1;
//     for (i = 0; i < np; i++) {
//       double bi = thetab[im];
//       for (jm = np - 1, j = 0; j < i; j++)
//         bi -= rbar[ithisr--] * P[jm--];
//       P[im--] = bi;
//     }
//
//     /*        now re-order p. */
//
//     ind = npr;
//     for (i = 0; i < r; i++) xnext[i] = P[ind++];
//     ind = np - 1;
//     ind1 = npr - 1;
//     for (i = 0; i < npr; i++) P[ind--] = P[ind1--];
//     for (i = 0; i < r; i++) P[i] = xnext[i];
//   } else {
//
//     /* P0 is obtained by backsubstitution for a moving average process. */
//
//     indn = np;
//     ind = np;
//     for (i = 0; i < r; i++)
//       for (j = 0; j <= i; j++) {
//         --ind;
//         P[ind] = V[ind];
//         if (j != 0) P[ind] += P[--indn];
//       }
//   }
//   /* now unpack to a full matrix */
//   for (i = r - 1, ind = np; i > 0; i--)
//     for (j = r - 1; j >= i; j--)
//       P[r * i + j] = P[--ind];
//   for (i = 0; i < r - 1; i++)
//     for (j = i + 1; j < r; j++)
//       P[i + r * j] = P[j + r * i];
//   UNPROTECT(1);
//   return res;
// }
