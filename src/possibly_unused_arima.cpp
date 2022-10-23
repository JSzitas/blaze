#include <Rcpp.h>
using namespace Rcpp;

// /* based on code from AS154 */
// I dont think this is actually used anywhere
// static void  inclu2(size_t np, double *xnext, double *xrow, double ynext,
//                     double *d, double *rbar, double *thetab)
// {
//     double cbar, sbar, di, xi, xk, rbthis, dpi;
//     size_t i, k, ithisr;
//
//     /*   This subroutine updates d, rbar, thetab by the inclusion
//      of xnext and ynext. */
//
//     for (i = 0; i < np; i++) xrow[i] = xnext[i];
//
//     for (ithisr = 0, i = 0; i < np; i++) {
//       if (xrow[i] != 0.0) {
//         xi = xrow[i];
//         di = d[i];
//         dpi = di + xi * xi;
//         d[i] = dpi;
//         cbar = di / dpi;
//         sbar = xi / dpi;
//         for (k = i + 1; k < np; k++) {
//           xk = xrow[k];
//           rbthis = rbar[ithisr];
//           xrow[k] = xk - xi * rbthis;
//           rbar[ithisr++] = cbar * rbthis + sbar * xk;
//         }
//         xk = ynext;
//         ynext = xk - xi * thetab[i];
//         thetab[i] = cbar * thetab[i] + sbar * xk;
//         if (di == 0.0) return;
//       } else
//         ithisr = ithisr + np - i - 1;
//     }
//   }
