#include <Rcpp.h>
using namespace Rcpp;

// SEXP
//   KalmanLike(SEXP sy, SEXP mod, SEXP sUP, SEXP op, SEXP update)
//   {
//     int lop = asLogical(op);
//     mod = PROTECT(duplicate(mod));
//
//     SEXP sZ = getListElement(mod, "Z"), sa = getListElement(mod, "a"),
//       sP = getListElement(mod, "P"), sT = getListElement(mod, "T"),
//       sV = getListElement(mod, "V"), sh = getListElement(mod, "h"),
//       sPn = getListElement(mod, "Pn");
//
//     if (TYPEOF(sy) != REALSXP || TYPEOF(sZ) != REALSXP ||
//         TYPEOF(sa) != REALSXP || TYPEOF(sP) != REALSXP ||
//         TYPEOF(sPn) != REALSXP ||
//         TYPEOF(sT) != REALSXP || TYPEOF(sV) != REALSXP)
//       error(_("invalid argument type"));
//
//     int n = LENGTH(sy), p = LENGTH(sa);
//     double *y = REAL(sy), *Z = REAL(sZ), *T = REAL(sT), *V = REAL(sV),
//       *P = REAL(sP), *a = REAL(sa), *Pnew = REAL(sPn), h = asReal(sh);
//
//     double *anew = (double *) R_alloc(p, sizeof(double));
//     double *M = (double *) R_alloc(p, sizeof(double));
//     double *mm = (double *) R_alloc(p * p, sizeof(double));
//     // These are only used if(lop), but avoid -Wall trouble
//     SEXP ans = R_NilValue, resid = R_NilValue, states = R_NilValue;
//     if(lop) {
//       PROTECT(ans = allocVector(VECSXP, 3));
//       SET_VECTOR_ELT(ans, 1, resid = allocVector(REALSXP, n));
//       SET_VECTOR_ELT(ans, 2, states = allocMatrix(REALSXP, n, p));
//       SEXP nm = PROTECT(allocVector(STRSXP, 3));
//       SET_STRING_ELT(nm, 0, mkChar("values"));
//       SET_STRING_ELT(nm, 1, mkChar("resid"));
//       SET_STRING_ELT(nm, 2, mkChar("states"));
//       setAttrib(ans, R_NamesSymbol, nm);
//       UNPROTECT(1);
//     }
//
//     double sumlog = 0.0, ssq = 0.0;
//     int nu = 0;
//     for (int l = 0; l < n; l++) {
//       for (int i = 0; i < p; i++) {
//         double tmp = 0.0;
//         for (int k = 0; k < p; k++)
//           tmp += T[i + p * k] * a[k];
//         anew[i] = tmp;
//       }
//       if (l > asInteger(sUP)) {
//         for (int i = 0; i < p; i++)
//           for (int j = 0; j < p; j++) {
//             double tmp = 0.0;
//             for (int k = 0; k < p; k++)
//               tmp += T[i + p * k] * P[k + p * j];
//             mm[i + p * j] = tmp;
//           }
//           for (int i = 0; i < p; i++)
//             for (int j = 0; j < p; j++) {
//               double tmp = V[i + p * j];
//               for (int k = 0; k < p; k++)
//                 tmp += mm[i + p * k] * T[j + p * k];
//               Pnew[i + p * j] = tmp;
//             }
//       }
//       if (!ISNAN(y[l])) {
//         nu++;
//         double *rr = NULL /* -Wall */;
//         if(lop) rr = REAL(resid);
//         double resid0 = y[l];
//         for (int i = 0; i < p; i++)
//           resid0 -= Z[i] * anew[i];
//         double gain = h;
//         for (int i = 0; i < p; i++) {
//           double tmp = 0.0;
//           for (int j = 0; j < p; j++)
//             tmp += Pnew[i + j * p] * Z[j];
//           M[i] = tmp;
//           gain += Z[i] * M[i];
//         }
//         ssq += resid0 * resid0 / gain;
//         if(lop) rr[l] = resid0 / sqrt(gain);
//         sumlog += log(gain);
//         for (int i = 0; i < p; i++)
//           a[i] = anew[i] + M[i] * resid0 / gain;
//         for (int i = 0; i < p; i++)
//           for (int j = 0; j < p; j++)
//             P[i + j * p] = Pnew[i + j * p] - M[i] * M[j] / gain;
//       } else {
//         double *rr = NULL /* -Wall */;
//         if(lop) rr = REAL(resid);
//         for (int i = 0; i < p; i++)
//           a[i] = anew[i];
//         for (int i = 0; i < p * p; i++)
//           P[i] = Pnew[i];
//         if(lop) rr[l] = NA_REAL;
//       }
//       if(lop) {
//         double *rs = REAL(states);
//         for (int j = 0; j < p; j++) rs[l + n*j] = a[j];
//       }
//     }
//
//     SEXP res = PROTECT(allocVector(REALSXP, 2));
//     REAL(res)[0] = ssq/nu; REAL(res)[1] = sumlog/nu;
//     if(lop) {
//       SET_VECTOR_ELT(ans, 0, res);
//       if(asLogical(update)) setAttrib(ans, install("mod"), mod);
//       UNPROTECT(3);
//       return ans;
//     } else {
//       if(asLogical(update)) setAttrib(res, install("mod"), mod);
//       UNPROTECT(2);
//       return res;
//     }
//   }
//
// SEXP
//   KalmanSmooth(SEXP sy, SEXP mod, SEXP sUP)
//   {
//     SEXP sZ = getListElement(mod, "Z"), sa = getListElement(mod, "a"),
//       sP = getListElement(mod, "P"), sT = getListElement(mod, "T"),
//       sV = getListElement(mod, "V"), sh = getListElement(mod, "h"),
//       sPn = getListElement(mod, "Pn");
//
//     if (TYPEOF(sy) != REALSXP || TYPEOF(sZ) != REALSXP ||
//         TYPEOF(sa) != REALSXP || TYPEOF(sP) != REALSXP ||
//         TYPEOF(sT) != REALSXP || TYPEOF(sV) != REALSXP)
//       error(_("invalid argument type"));
//
//     SEXP ssa, ssP, ssPn, res, states = R_NilValue, sN;
//     int n = LENGTH(sy), p = LENGTH(sa);
//     double *y = REAL(sy), *Z = REAL(sZ), *a, *P,
//       *T = REAL(sT), *V = REAL(sV), h = asReal(sh), *Pnew;
//     double *at, *rt, *Pt, *gains, *resids, *Mt, *L, gn, *Nt;
//     Rboolean var = TRUE;
//
//     PROTECT(ssa = duplicate(sa)); a = REAL(ssa);
//     PROTECT(ssP = duplicate(sP)); P = REAL(ssP);
//     PROTECT(ssPn = duplicate(sPn)); Pnew = REAL(ssPn);
//
//     PROTECT(res = allocVector(VECSXP, 2));
//     SEXP nm = PROTECT(allocVector(STRSXP, 2));
//     SET_STRING_ELT(nm, 0, mkChar("smooth"));
//     SET_STRING_ELT(nm, 1, mkChar("var"));
//     setAttrib(res, R_NamesSymbol, nm);
//     UNPROTECT(1);
//     SET_VECTOR_ELT(res, 0, states = allocMatrix(REALSXP, n, p));
//     at = REAL(states);
//     SET_VECTOR_ELT(res, 1, sN = allocVector(REALSXP, n*p*p));
//     Nt = REAL(sN);
//
//     double *anew, *mm, *M;
//     anew = (double *) R_alloc(p, sizeof(double));
//     M = (double *) R_alloc(p, sizeof(double));
//     mm = (double *) R_alloc(p * p, sizeof(double));
//
//     Pt = (double *) R_alloc(n * p * p, sizeof(double));
//     gains = (double *) R_alloc(n, sizeof(double));
//     resids = (double *) R_alloc(n, sizeof(double));
//     Mt = (double *) R_alloc(n * p, sizeof(double));
//     L = (double *) R_alloc(p * p, sizeof(double));
//
//     for (int l = 0; l < n; l++) {
//       for (int i = 0; i < p; i++) {
//         double tmp = 0.0;
//         for (int k = 0; k < p; k++)
//           tmp += T[i + p * k] * a[k];
//         anew[i] = tmp;
//       }
//       if (l > asInteger(sUP)) {
//         for (int i = 0; i < p; i++)
//           for (int j = 0; j < p; j++) {
//             double tmp = 0.0;
//             for (int k = 0; k < p; k++)
//               tmp += T[i + p * k] * P[k + p * j];
//             mm[i + p * j] = tmp;
//           }
//           for (int i = 0; i < p; i++)
//             for (int j = 0; j < p; j++) {
//               double tmp = V[i + p * j];
//               for (int k = 0; k < p; k++)
//                 tmp += mm[i + p * k] * T[j + p * k];
//               Pnew[i + p * j] = tmp;
//             }
//       }
//       for (int i = 0; i < p; i++) at[l + n*i] = anew[i];
//       for (int i = 0; i < p*p; i++) Pt[l + n*i] = Pnew[i];
//       if (!ISNAN(y[l])) {
//         double resid0 = y[l];
//         for (int i = 0; i < p; i++)
//           resid0 -= Z[i] * anew[i];
//         double gain = h;
//         for (int i = 0; i < p; i++) {
//           double tmp = 0.0;
//           for (int j = 0; j < p; j++)
//             tmp += Pnew[i + j * p] * Z[j];
//           Mt[l + n*i] = M[i] = tmp;
//           gain += Z[i] * M[i];
//         }
//         gains[l] = gain;
//         resids[l] = resid0;
//         for (int i = 0; i < p; i++)
//           a[i] = anew[i] + M[i] * resid0 / gain;
//         for (int i = 0; i < p; i++)
//           for (int j = 0; j < p; j++)
//             P[i + j * p] = Pnew[i + j * p] - M[i] * M[j] / gain;
//       } else {
//         for (int i = 0; i < p; i++) {
//           a[i] = anew[i];
//           Mt[l + n * i] = 0.0;
//         }
//         for (int i = 0; i < p * p; i++)
//           P[i] = Pnew[i];
//         gains[l] = NA_REAL;
//         resids[l] = NA_REAL;
//       }
//     }
//
//     /* rt stores r_{t-1} */
//     rt = (double *) R_alloc(n * p, sizeof(double));
//     for (int l = n - 1; l >= 0; l--) {
//       if (!ISNAN(gains[l])) {
//         gn = 1/gains[l];
//         for (int i = 0; i < p; i++)
//           rt[l + n * i] = Z[i] * resids[l] * gn;
//       } else {
//         for (int i = 0; i < p; i++) rt[l + n * i] = 0.0;
//         gn = 0.0;
//       }
//
//       if (var) {
//         for (int i = 0; i < p; i++)
//           for (int j = 0; j < p; j++)
//             Nt[l + n*i + n*p*j] = Z[i] * Z[j] * gn;
//       }
//
//       if (l < n - 1) {
//         /* compute r_{t-1} */
//         for (int i = 0; i < p; i++)
//           for (int j = 0; j < p; j++)
//             mm[i + p * j] = ((i==j) ? 1:0) - Mt[l + n * i] * Z[j] * gn;
//         for (int i = 0; i < p; i++)
//           for (int j = 0; j < p; j++) {
//             double tmp = 0.0;
//             for (int k = 0; k < p; k++)
//               tmp += T[i + p * k] * mm[k + p * j];
//             L[i + p * j] = tmp;
//           }
//           for (int i = 0; i < p; i++) {
//             double tmp = 0.0;
//             for (int j = 0; j < p; j++)
//               tmp += L[j + p * i] * rt[l + 1 + n * j];
//             rt[l + n * i] += tmp;
//           }
//           if(var) { /* compute N_{t-1} */
//         for (int i = 0; i < p; i++)
//           for (int j = 0; j < p; j++) {
//             double tmp = 0.0;
//             for (int k = 0; k < p; k++)
//               tmp += L[k + p * i] * Nt[l + 1 + n*k + n*p*j];
//             mm[i + p * j] = tmp;
//           }
//           for (int i = 0; i < p; i++)
//             for (int j = 0; j < p; j++) {
//               double tmp = 0.0;
//               for (int k = 0; k < p; k++)
//                 tmp += mm[i + p * k] * L[k + p * j];
//               Nt[l + n*i + n*p*j] += tmp;
//             }
//           }
//       }
//
//       for (int i = 0; i < p; i++) {
//         double tmp = 0.0;
//         for (int j = 0; j < p; j++)
//           tmp += Pt[l + n*i + n*p*j] * rt[l + n * j];
//         at[l + n*i] += tmp;
//       }
//     }
//     if (var)
//       for (int l = 0; l < n; l++) {
//         for (int i = 0; i < p; i++)
//           for (int j = 0; j < p; j++) {
//             double tmp = 0.0;
//             for (int k = 0; k < p; k++)
//               tmp += Pt[l + n*i + n*p*k] * Nt[l + n*k + n*p*j];
//             mm[i + p * j] = tmp;
//           }
//           for (int i = 0; i < p; i++)
//             for (int j = 0; j < p; j++) {
//               double tmp = Pt[l + n*i + n*p*j];
//               for (int k = 0; k < p; k++)
//                 tmp -= mm[i + p * k] * Pt[l + n*k + n*p*j];
//               Nt[l + n*i + n*p*j] = tmp;
//             }
//       }
//       UNPROTECT(4);
//     return res;
//   }


// SEXP
//   KalmanFore(SEXP nahead, SEXP mod, SEXP update)
//   {
//     mod = PROTECT(duplicate(mod));
//     SEXP sZ = getListElement(mod, "Z"), sa = getListElement(mod, "a"),
//       sP = getListElement(mod, "P"), sT = getListElement(mod, "T"),
//       sV = getListElement(mod, "V"), sh = getListElement(mod, "h");
//
//     if (TYPEOF(sZ) != REALSXP ||
//         TYPEOF(sa) != REALSXP || TYPEOF(sP) != REALSXP ||
//         TYPEOF(sT) != REALSXP || TYPEOF(sV) != REALSXP)
//       error(_("invalid argument type"));
//
//     int  n = asInteger(nahead), p = LENGTH(sa);
//     double *Z = REAL(sZ), *a = REAL(sa), *P = REAL(sP), *T = REAL(sT),
//       *V = REAL(sV), h = asReal(sh);
//     double *mm, *anew, *Pnew;
//
//     anew = (double *) R_alloc(p, sizeof(double));
//     Pnew = (double *) R_alloc(p * p, sizeof(double));
//     mm = (double *) R_alloc(p * p, sizeof(double));
//     SEXP res, forecasts, se;
//     PROTECT(res = allocVector(VECSXP, 2));
//     SET_VECTOR_ELT(res, 0, forecasts = allocVector(REALSXP, n));
//     SET_VECTOR_ELT(res, 1, se = allocVector(REALSXP, n));
//     {
//       SEXP nm = PROTECT(allocVector(STRSXP, 2));
//       SET_STRING_ELT(nm, 0, mkChar("pred"));
//       SET_STRING_ELT(nm, 1, mkChar("var"));
//       setAttrib(res, R_NamesSymbol, nm);
//       UNPROTECT(1);
//     }
//     for (int l = 0; l < n; l++) {
//       double fc = 0.0;
//       for (int i = 0; i < p; i++) {
//         double tmp = 0.0;
//         for (int k = 0; k < p; k++)
//           tmp += T[i + p * k] * a[k];
//         anew[i] = tmp;
//         fc += tmp * Z[i];
//       }
//       for (int i = 0; i < p; i++)
//         a[i] = anew[i];
//       REAL(forecasts)[l] = fc;
//
//       for (int i = 0; i < p; i++)
//         for (int j = 0; j < p; j++) {
//           double tmp = 0.0;
//           for (int k = 0; k < p; k++)
//             tmp += T[i + p * k] * P[k + p * j];
//           mm[i + p * j] = tmp;
//         }
//         for (int i = 0; i < p; i++)
//           for (int j = 0; j < p; j++) {
//             double tmp = V[i + p * j];
//             for (int k = 0; k < p; k++)
//               tmp += mm[i + p * k] * T[j + p * k];
//             Pnew[i + p * j] = tmp;
//           }
//           double tmp = h;
//       for (int i = 0; i < p; i++)
//         for (int j = 0; j < p; j++) {
//           P[i + j * p] = Pnew[i + j * p];
//           tmp += Z[i] * Z[j] * P[i + j * p];
//         }
//         REAL(se)[l] = tmp;
//     }
//     if(asLogical(update)) setAttrib(res, install("mod"), mod);
//     UNPROTECT(2);
//     return res;
//   }
