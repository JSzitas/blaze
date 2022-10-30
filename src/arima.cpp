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

