// Copyright 2020, https://github.com/PatWie/CppNumericalSolvers
//
#ifndef INCLUDE_CPPOPTLIB_LINESEARCH_MORE_THUENTE_H_
#define INCLUDE_CPPOPTLIB_LINESEARCH_MORE_THUENTE_H_

#include <algorithm>
#include <cmath>

namespace cppoptlib::solver::linesearch {

template <typename Function, int Ord>
class MoreThuente {
 public:
  using scalar_t = typename Function::scalar_t;
  using vector_t = typename Function::vector_t;
  using state_t = typename Function::state_t;
  /**
   * @brief use MoreThuente Rule for (strong) Wolfe conditiions
   * @details [long description]
   *
   * @param search_direction search direction for next update step
   * @param function handle to problem
   *
   * @return step-width
   */
  static scalar_t Search(const state_t &state, const vector_t &search_direction,
                         Function &function,
                         const scalar_t alpha_init = 1.0) {
    scalar_t alpha = alpha_init;
    vector_t g = state.gradient, s = search_direction.eval(), xx = state.x;

    cvsrch(function, &xx, state.value, &g, &alpha, s);
    return alpha;
  }
  // Note that the result of this gets discarded, so for all intends and purposes 
  // this should be type void - or perhaps directly return alpha
  static void cvsrch(Function &function, vector_t *x, scalar_t f,
                     vector_t *g, scalar_t *stp, const vector_t &s) {
    // we rewrite this from MIN-LAPACK and some MATLAB code
    // TODO:  are the info variables necessary - see return type comment
    // int info = 0;
    // int infoc = 1;
    constexpr scalar_t xtol = 1e-15, ftol = 1e-4, gtol = 1e-2, 
      stpmin = 1e-15, stpmax = 1e15, xtrapf = 4;
    constexpr int maxfev = 20;
    int nfev = 0;

    scalar_t dginit = g->dot(s);
    // There is no descent direction. TODO: handle
    if (dginit >= 0.0) return;
    bool brackt = false, stage1 = true;
    
    scalar_t finit = f, dgtest = ftol * dginit, width = stpmax - stpmin,
      width1 = 2 * width, stx = 0.0, fx = finit, dgx = dginit, sty = 0.0,
      fy = finit, dgy = dginit, stmin, stmax;
    
    vector_t wa = x->eval();
    while (true) {
      // Make sure we stay in the interval when setting min/max-step-width.
      if (brackt) {
        stmin = std::min<scalar_t>(stx, sty);
        stmax = std::max<scalar_t>(stx, sty);
      } else {
        stmin = stx;
        stmax = *stp + xtrapf * (*stp - stx);
      }

      // Force the step to be within the bounds stpmax and stpmin.
      *stp = std::clamp(*stp, stpmin, stpmax);

      // Oops, let us return the last reliable values.
      if ((brackt && ((*stp <= stmin) || (*stp >= stmax))) ||
          (nfev >= maxfev - 1) || (
            brackt && ((stmax - stmin) <= (xtol * stmax)))) {
        *stp = stx;
      }

      // Test new point.
      *x = wa + *stp * s;
      f = function(*x);
      function.Gradient(*x, g);
      nfev++;
      scalar_t dg = g->dot(s);
      scalar_t ftest1 = finit + *stp * dgtest;

      // All possible convergence tests - optimizable since they run in order
      if ((nfev >= maxfev) |
          (brackt & ((*stp <= stmin) | (*stp >= stmax))) |
          (*stp == stpmax) & (f <= ftest1) & (dg <= dgtest) |
          (*stp == stpmin) & ((f > ftest1) | (dg >= dgtest)) |
          (brackt & (stmax - stmin <= xtol * stmax)) |
          ((f <= ftest1) & (fabs(dg) <= gtol * (-dginit)))) return;

      if (stage1 & (f <= ftest1) &
          (dg >= std::min<scalar_t>(ftol, gtol) * dginit)) stage1 = false;

      if (stage1 & (f <= fx) & (f > ftest1)) {
        scalar_t fm = f - *stp * dgtest,
          fxm = fx - stx * dgtest, fym = fy - sty * dgtest,
          dgm = dg - dgtest, dgxm = dgx - dgtest, dgym = dgy - dgtest;

        cstep(stx, fxm, dgxm, sty, fym, dgym, *stp, fm, dgm, brackt, stmin,
              stmax);
        // updates 
        fx = fxm + stx * dgtest;
        fy = fym + sty * dgtest;
        dgx = dgxm + dgtest;
        dgy = dgym + dgtest;
      } else {
        // This is ugly and some variables should be moved to the class scope.
        cstep(stx, fx, dgx, sty, fy, dgy, *stp, f, dg, brackt, stmin, stmax);
      }
      if (brackt) {
        if (fabs(sty - stx) >= 0.66 * width1) {
          *stp = stx + 0.5 * (sty - stx);
        }
        width1 = width;
        width = fabs(sty - stx);
      }
    }
    return;
  }
  // TODO(patwie): cpplint prefers pointers here, but this would make the code
  // unreadable. As these are all changing values a configuration structure
  // would be helpful.
  static void cstep(scalar_t &stx, scalar_t &fx, scalar_t &dx,
                    scalar_t &sty,                                
                    scalar_t &fy, scalar_t &dy, scalar_t &stp,
                    scalar_t &fp, scalar_t &dp,
                    bool &brackt,
                    scalar_t &stpmin, scalar_t &stpmax) {
    bool bound = false;
    // Check the input parameters for errors.
    if ((brackt & ((stp <= std::min<scalar_t>(stx, sty)) |
                   (stp >= std::max<scalar_t>(stx, sty)))) |
        (dx * (stp - stx) >= 0.0) | (stpmax < stpmin)) {
      return;
    }

    scalar_t sgnd = dp * (dx / fabs(dx));
    scalar_t stpf = 0, stpc = 0, stpq = 0;

    if (fp > fx) {
      bound = true;
      scalar_t theta = 3. * (fx - fp) / (stp - stx) + dx + dp;
      scalar_t s = max_abs(theta, dx, dp);
      scalar_t gamma =
          s * sqrt((theta / s) * (theta / s) - (dx / s) * (dp / s));
      if (stp < stx) gamma = -gamma;
      scalar_t p = (gamma - dx) + theta;
      scalar_t q = ((gamma - dx) + gamma) + dp;
      scalar_t r = p / q;
      stpc = stx + r * (stp - stx);
      stpq = stx + ((dx / ((fx - fp) / (stp - stx) + dx)) / 2.) * (stp - stx);
      if (fabs(stpc - stx) < fabs(stpq - stx))
        stpf = stpc;
      else
        stpf = stpc + (stpq - stpc) / 2;
      brackt = true;
    } else if (sgnd < 0.0) {
      scalar_t theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
      scalar_t s = max_abs(theta, dx, dp);
      scalar_t gamma =
          s * sqrt((theta / s) * (theta / s) - (dx / s) * (dp / s));
      if (stp > stx) gamma = -gamma;

      scalar_t p = (gamma - dp) + theta;
      scalar_t q = ((gamma - dp) + gamma) + dx;
      scalar_t r = p / q;
      stpc = stp + r * (stx - stp);
      stpq = stp + (dp / (dp - dx)) * (stx - stp);
      if (fabs(stpc - stp) > fabs(stpq - stp))
        stpf = stpc;
      else
        stpf = stpq;
      brackt = true;
    } else if (fabs(dp) < fabs(dx)) {
      bound = true;
      scalar_t theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
      scalar_t s = max_abs(theta, dx, dp);
      scalar_t gamma = s * sqrt(std::max<scalar_t>(static_cast<scalar_t>(0.),
                                                   (theta / s) * (theta / s) -
                                                       (dx / s) * (dp / s)));
      if (stp > stx) gamma = -gamma;
      scalar_t p = (gamma - dp) + theta;
      scalar_t q = (gamma + (dx - dp)) + gamma;
      scalar_t r = p / q;
      if ((r < 0.0) & (gamma != 0.0)) {
        stpc = stp + r * (stx - stp);
      } else if (stp > stx) {
        stpc = stpmax;
      } else {
        stpc = stpmin;
      }
      stpq = stp + (dp / (dp - dx)) * (stx - stp);
      if (brackt) {
        if (fabs(stp - stpc) < fabs(stp - stpq)) {
          stpf = stpc;
        } else {
          stpf = stpq;
        }
      } else {
        if (fabs(stp - stpc) > fabs(stp - stpq)) {
          stpf = stpc;
        } else {
          stpf = stpq;
        }
      }
    } else {
      if (brackt) {
        scalar_t theta = 3 * (fp - fy) / (sty - stp) + dy + dp;
        scalar_t s = max_abs(theta, dy, dp);
        scalar_t gamma =
            s * sqrt((theta / s) * (theta / s) - (dy / s) * (dp / s));
        if (stp > sty) gamma = -gamma;

        scalar_t p = (gamma - dp) + theta;
        scalar_t q = ((gamma - dp) + gamma) + dy;
        scalar_t r = p / q;
        stpc = stp + r * (sty - stp);
        stpf = stpc;
      } else if (stp > stx) {
        stpf = stpmax;
      } else {
        stpf = stpmin;
      }
    }

    if (fp > fx) {
      sty = stp;
      fy = fp;
      dy = dp;
    } else {
      if (sgnd < 0.0) {
        sty = stx;
        fy = fx;
        dy = dx;
      }
      stx = stp;
      fx = fp;
      dx = dp;
    }
    stp = std::clamp(stpf, stpmin, stpmax);

    if (brackt & bound) {
      if (sty > stx) {
        stp = std::min<scalar_t>(
            stx + static_cast<scalar_t>(0.66) * (sty - stx), stp);
      } else {
        stp = std::max<scalar_t>(
            stx + static_cast<scalar_t>(0.66) * (sty - stx), stp);
      }
    }
    return;
  }

  static scalar_t max_abs(scalar_t x, scalar_t y, scalar_t z) {
    return std::max(std::abs(x), std::max(std::abs(y), std::abs(z)));
  }
};
}  // namespace cppoptlib::solver::linesearch

#endif  // INCLUDE_CPPOPTLIB_LINESEARCH_MORE_THUENTE_H_
