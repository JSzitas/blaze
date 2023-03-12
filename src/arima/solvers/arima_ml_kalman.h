#ifndef ARIMA_ML_KALMAN
#define ARIMA_ML_KALMAN

#include "arima/structures/structural_model.h"
#include "arima/structures/arima_kind.h"
#include "arima/structures/ss_init.h"

#include "arima/utils/transforms.h"
#include "arima/utils/delta.h"

#include "third_party/eigen.h"
#include "utils/utils.h"

#include "arima/solvers/state_space.h"


// main class for doing Kalman filtering of an ARIMA model

template <const SSinit ss_type, const bool seasonal,
          const bool has_xreg, const bool transform,
          const bool update_P=true, typename scalar_t=double> class KalmanARIMA {

using EigVec = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;
using EigMat = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;
using vec = std::vector<scalar_t>;

// these things should be const
vec y;
size_t n, arma_pars, p, q, r, d, rd, r2;
arima_kind kind;
scalar_t kappa;

vec anew, M, mm, gam, g, rrz;
EigVec u;
EigMat xreg;
size_t np, nrbar, npr, npr1;
vec xnext, xrow, rbar, thetab;

EigVec y_temp, new_x;
vec residual, temp_phi, temp_theta;

structural_model<scalar_t> model;
scalar_t sigma2;

public:
  KalmanARIMA<ss_type, seasonal, has_xreg, transform, update_P, scalar_t>(){};
  KalmanARIMA<ss_type, seasonal, has_xreg, transform, update_P, scalar_t>(
      const vec &y,
      const arima_kind &kind,
      EigMat xreg,
      const scalar_t kappa ) : y(y), n(y.size()),
      arma_pars(kind.p() + kind.P() + kind.q() + kind.Q()),
      p(kind.p() + (kind.P() * kind.period())),
      q(kind.q() + (kind.Q() * kind.period())),
      r(max(p,q+1)),
      d((kind.d()+1) +(kind.period() * kind.D())),
      rd(r+d),
      r2(max(p + q, p + 1)),
      kind(kind), kappa(kappa), xreg(xreg) {

    this->anew = vec(rd);
    this->M = vec(rd);
    this->mm = vec( (d > 0) * rd * rd);

    if constexpr( ss_type == SSinit::Rossignol) {
      this->gam = vec(r2 * r2);
      this->g = vec(r2);
      this->rrz = vec(q);
      this->u = EigVec(r2);
    }
    if constexpr( ss_type == SSinit::Gardner) {
      // all of these following things are static
      this->np = this->r * (this->r + 1) / 2;
      this->nrbar = this->np * (this->np - 1) / 2;
      this->npr = this->np - this->r;
      // preallocate expansion vectors
      this->xnext = vec(this->np);
      this->xrow = vec(this->np);
      this->rbar = vec(this->nrbar);
      this->thetab = vec(this->np);
    }
    this->y_temp = EigVec(this->n);
    for (size_t i = 0; i < this->n; i++) this->y_temp(i) = y[i];
    // pre-allocate new_x
    this->new_x = EigVec(kind.p() + (kind.P() * kind.period()) + kind.q() +
      (kind.Q() * kind.period()) + this->xreg.cols());
    // pre-allocate model residuals
    this->residual = vec(this->n);
    // pre-allocate transformation helper vector - this is only necessary
    // for expanding seasonal models
    this->temp_phi = vec(kind.p() + (kind.P() * kind.period()));
    this->temp_theta = vec(kind.q() + (kind.Q() * kind.period()));
    if constexpr(seasonal || transform) {
      // the expansion of arima parameters is only necessary for seasonal models
      arima_transform_parameters<EigVec, seasonal, transform>(
          this->new_x, this->kind, this->temp_phi, this->temp_theta
      );
    }
    // // initialize state space model

    if constexpr( ss_type == SSinit::Gardner ) {
      this->model = make_arima( this->new_x, kind, kappa, SSinit::Gardner);
    }
    if constexpr( ss_type == SSinit::Rossignol) {
      this->model = make_arima( this->new_x, kind, kappa, SSinit::Rossignol);
    }
    this->sigma2 = 0;
  };
  scalar_t operator()(const EigVec &x) {
    for (size_t i = 0; i < x.size(); i++) this->new_x(i) = x(i);
    if constexpr (has_xreg) {
      // refresh y_temp and load it with original y data
      for (size_t i = 0; i < this->n; i++) this->y_temp(i) = this->y[i];
      this->y_temp = this->y_temp - this->xreg * x.tail(x.size() - this->arma_pars);
    }
    /* I figured out that I can basically expand this out altogether for non-seasonal
     * models - the compiler should insert an empty function anyways, but just to
     * make sure that this gets compiled away - we can make sure its a dead branch
     */
    if constexpr(seasonal || transform) {
      // the expansion of arima parameters is only necessary for seasonal models
      arima_transform_parameters<EigVec, seasonal, transform>(
          this->new_x, this->kind, this->temp_phi, this->temp_theta
      );
    }
    // update arima
    this->update_arima(this->new_x);
    // return likelihood - check if this updated model or not (it ideally
    // should not, not here)
    // const std::array<scalar_t,2> res = this->arima_likelihood(this->y_temp, this->model);
    this->sigma2 = 1.0;//res[0];
    return 0.0;
  }
  scalar_t get_sigma() const { return this->sigma2; }
  structural_model<scalar_t> get_structural_model() const { return this->model; }
private:
  // rewrite
  // template <typename C>
  // void make_arima( const C &coef,
  //                  const scalar_t tol = 1e-9) {
  //   vec delta = make_delta(this->kind.d(), kind.D(), kind.period());
  //
  //   vec phi = vec(this->p);
  //   vec theta = vec(this->q + max(this->r - 1 - this->q, 0));
  //   // copz out elements of coef into phi and theta
  //   for( size_t i = 0; i < this->p; i++) phi[i] = coef[i];
  //   for( size_t i = p; i < this->p + this->q; i++) theta[i-p] = coef[i];
  //   size_t i, j;
  //   auto Z = vec(rd);
  //   Z[0] = 1;
  //   for (i = 1; i < r - 1; i++) Z[i] = 0;
  //   j = 0;
  //   for (i = r; i < rd; i++) {
  //     Z[i] = delta[j];
  //     j++;
  //   }
  //   auto T = vec(rd * rd);
  //   if (p > 0) {
  //     for (i = 0; i < p; i++) {
  //       T[i] = phi[i];
  //     }
  //   }
  //   if (r > 1L) {
  //     /* set '2nd diagonal' elements to 1. since this is hard to understand,
  //      * here are some examples
  //      * for a matrix with *rd* == 5, transform it:
  //      *
  //      *   (input)          (rd == r)         (rd > r)
  //      *   x x x x x   =>   x 1 x x x  //  => x 1 x x x
  //      *   x x x x x   =>   x x 1 x x  //  => x x 1 x x
  //      *   x x x x x   =>   x x x 1 x  //  => x x x 1 x
  //      *   x x x x x   =>   x x x x 1  //  => x x x x x
  //      *   x x x x x   =>   x x x x x  //  => x x x x x
  //      */
  //     for (i = 0; i < r - 1; i++) {
  //       T[(rd * (i + 1)) + i] = 1;
  //     }
  //   }
  //   if (d > 0) {
  //     // replace row r+1 in R (or row r for us)
  //     // with whatever is in Z
  //     for (j = 0; j < rd; j++) {
  //       T[(j * rd) + r] = Z[j];
  //     }
  //     // if there are more than 1 differences d
  //     if (d > 1) {
  //       /* start replacing at r and continue until you get to r + d (-1?)
  //        * replace r + 2:d with 1 - this is similar as above, but it accounts
  //        * for the first difference differently(that is taken care of in the above
  //        * code, in the j < rd loop). here we are taking care of the other
  //        * differences, so if we have 3 differences, we will only be doing 2
  //        * replacements. these happen after the for a matrix with *rd* == 5, with
  //        * *d* == 3, transform it: (input)          (d == 3) x x x x x   =>   x x
  //        * x x x x x x x x   =>   x x x x x x x x x x   =>   x x x x x x x x x x
  //        * =>   x x 1 x x x x x x x   =>   x x x 1 x
  //        */
  //       for (i = r; i < rd - 1; i++) {
  //         T[((rd + 1) * i) + 1] = 1;
  //       }
  //     }
  //   }
  //   // this is R <- c(1, theta, rep.int(0, d))
  //   // we can skip the d part as vectors are 0 initialized.
  //   vec R = vec(1 + theta.size() + d);
  //   R[0] = 1;
  //   for (i = 1; i < theta.size() + 1; i++) {
  //     R[i] = theta[i - 1];
  //   }
  //   const auto r_sz = R.size();
  //   auto V = vec(r_sz * r_sz);
  //   // here we do an outer product, ie: V <- R %o% R
  //   size_t mat_p = 0;
  //   for (i = 0; i < r_sz; i++) {
  //     for (j = 0; j < r_sz; j++) {
  //       V[mat_p] = R[i] * R[j];
  //       mat_p++;
  //     }
  //   }
  //   scalar_t h = 0.0;
  //   auto a = vec(rd);
  //   auto P = vec(rd * rd);
  //   auto Pn = vec(rd * rd);
  //   // create model
  //   this->model = structural_model<scalar_t>(
  //     phi, theta, delta, Z, a, P, T, V, Pn, h);
  //   // compute new Pn
  //   if (r > 1) {
  //     if constexpr(update_P) {
  //       if constexpr(ss_type == SSinit::Rossignol) {
  //         get_Q0_rossignol();
  //       }
  //       if constexpr(ss_type == SSinit::Gardner) {
  //         get_Q0();
  //       }
  //     }
  //   } else {
  //     this->model.Pn[0] = (p > 0) * (1 / (1 - pow(this->model.phi[0], 2))) + (p == 0);
  //   }
  //   if (d > 0L) {
  //     /* update diagonal elements which come after the coefficients -
  //      * diagonal entries between r and rd - with kappa
  //      */
  //     for (i = r; i < rd; i++) {
  //       for (j = r; j < rd; j++) {
  //         // to only update diagonal elements check that we are on the diagonal
  //         // otherwise we have a zero - as intended
  //         this->model.Pn[(j * rd) + i] = (i == j) * this->kappa;
  //       }
  //     }
  //   }
  // }
  // map 2 vectors to Eigen matrices and call solve
  void solve_mat_vec(vec &mat, vec &vec) {
    EigMat new_mat = Eigen::Map<EigMat>(mat.data(), r2, r2);
    EigVec new_vec = Eigen::Map<EigVec>(vec.data(), r2, 1);
    this->u = new_mat.completeOrthogonalDecomposition().solve(new_vec);
  }
  void get_Q0_rossignol() {
    /* A mild reimplementation of the Matwey V. Kornilov's implementation of
     * algorithm by Dr. Raphael Rossignol, See
     * https://bugs.r-project.org/show_bug.cgi?id=14682 for details (of algorithm).
     * Avoids R specific data structures (i.e. SEXPs and related types) in favor
     * of standard C++, uses Eigen (rather than a call to the solver used by R),
     * to solve a system of linear equations, where the Eigen solver might be faster
     * (ie completeOrthogonalDecomposition **should** be a better approach than the
     * regular solver)
     */
    size_t i, j;

    // Initialize and allocate theta
    // UNNECESSARY_ALLOC
    vec theta_t(q + 1);
    theta_t[0] = 1.0;
    for (i = 1; i < q + 1; i++) theta_t[i] = model.theta[i - 1];

    if (p > 0) {
      // initialize phi
      // UNNECESSARY_ALLOC
      vec phi_t(p + 1);
      // fill with correct values (i.e. a leading 1)
      phi_t[0] = 1.0;
      for (i = 1; i < p + 1; i++) phi_t[i] = -model.phi[i - 1];
      /* Compute the autocovariance function of U, the AR part of X */
      /* Gam := C1 + C2 ; initialize */
      // gam is a matrix
      /* C1[E] */
      for (j = 0; j < r2; ++j) {
        for (i = j; i < r2 && i - j < p + 1; ++i) {
          this->gam[j * r2 + i] += phi_t[i - j];
        }
      }
      /* C2[E] */
      for (i = 0; i < r2; ++i) {
        for (j = 1; j < r2 && i + j < p + 1; ++j) {
          this->gam[j * r2 + i] += phi_t[i + j];
        }
      }
      // reset g => [1, 0, 0, 0, ..., 0]
      std::fill(this->g.begin(), this->g.end(), 0.0);
      this->g[0] = 1.;
      // Solve the system of linear equations - using Eigen
      this->solve_mat_vec(this->gam, this->g);
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
                  this->model.Pn[r * i + j] += phi_t[i + k] * phi_t[j + m] * theta_t[L - k] *
                    theta_t[n - m] * this->u[abs(L - n)];
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

      // I think this condition is unnecessary
      if (q > 0) {
        for (i = 0; i < q; ++i) {
          this->rrz[i] = theta_t[i];
          for (j = max(0, i - p); j < i; ++j) {
            this->rrz[i] -= this->rrz[j] * phi_t[i - j];
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
              this->model.Pn[r * i + j] += phi_t[i + k] * theta_t[j + L] * this->rrz[L - k - 1];
            }
          }
          for (k = 0; j + k < p; ++k) {
            for (L = k + 1; i + L < q + 1; ++L) {
              this->model.Pn[r * i + j] += phi_t[j + k] * theta_t[i + L] * this->rrz[L - k - 1];
            }
          }
        }
    } // end if(p > 0)
    /* Q0 += A2 A2^T */
    for (i = 0; i < r; ++i) {
      for (j = i; j < r; ++j) {
        for (size_t k = 0; j + k < q + 1; ++k) {
          this->model.Pn[r * i + j] += theta_t[i + k] * theta_t[j + k];
        }
      }
    }
    /* Symmetrize result */
    for (i = 0; i < r; ++i) {
      for (j = i + 1; j < r; ++j) {
        this->model.Pn[r * j + i] = this->model.Pn[r * i + j];
      }
    }
  }
  /* based on code from AS154 */
  inline void inclu2( scalar_t ynext ) {
    /*   This subroutine updates d, rbar, thetab by the inclusion
     of xnext and ynext. */
    for (size_t i = 0; i < this->np; i++) this->xrow[i] = this->xnext[i];

    scalar_t cbar, sbar, di, xi, xk, rbthis, dpi;
    for (size_t ithisr = 0, i = 0; i < this->np; i++) {
      if (this->xrow[i] != 0.0) {
        xi = this->xrow[i];
        // TODO: this inclu_d is probably unnecessary and needs to be replaced by smth
        di = this->model.Pn[i];
        dpi = di + xi * xi;
        this->model.Pn[i] = dpi;
        cbar = di / dpi;
        sbar = xi / dpi;
        for (size_t k = i + 1; k < this->np; k++) {
          xk = this->xrow[k];
          rbthis = this->rbar[ithisr];
          this->xrow[k] = xk - xi * rbthis;
          this->rbar[ithisr++] = cbar * rbthis + sbar * xk;
        }
        xk = ynext;
        ynext = xk - xi * this->thetab[i];
        this->thetab[i] = cbar * this->thetab[i] + sbar * xk;
        // throws uninitialized value error
        if (this->model.Pn[i] == 0.0) {
          return;
        }
      } else {
        ithisr = ithisr + this->np - i - 1;
      }
    }
  }
  void get_Q0() {
    size_t i, j, indi, indj, indn, ithisr, ind, ind1, ind2, im, jm;

    scalar_t vj, vi, bi, ynext, phii, phij;
    for (ind = 0, j = 0; j < r; j++) {
      vj = 0.0;
      if (j == 0) vj = 1.0;
      else if (j - 1 < q) vj = this->model.theta[j - 1];
      for (i = j; i < r; i++) {
        vi = 0.0;
        if (i == 0) vi = 1.0;
        else if (i - 1 < q) vi = this->model.theta[i - 1];
        this->model.V[ind++] = vi * vj;
      }
    }
    if (r == 1) {
      this->model.Pn[0] = (scalar_t)(p == 0) +
        (scalar_t)(p != 0) / (1.0 - this->model.phi[0] * this->model.phi[0]);
      return;
    }
    if (p > 0) {
      /*      The set of equations s * vec(P0) = vec(v) is solved for
       vec(P0).  s is generated row by row in the array xnext.  The
       order of elements in P is changed, so as to bring more leading
       zeros into the rows of s. */
      std::fill(this->rbar.begin(), this->rbar.end(), 0.0);
      std::fill(this->model.Pn.begin(), this->model.Pn.end(), 0.0);
      std::fill(this->thetab.begin(), this->thetab.end(), 0.0);
      std::fill(this->xnext.begin(), this->xnext.end(), 0.0);

      ind = 0;
      ind1 = -1;
      indj = npr;
      ind2 = npr - 1;

      indi = npr + 1;
      for (j = 0; j < r; j++) {
        phij = (j < p) ? this->model.phi[j] : 0.0;
        this->xnext[indj++] = 0.0;
        indi++;
        for (i = j; i < r; i++) {
          ynext = this->model.V[ind++];
          phii = (i < p) ? this->model.phi[i] : 0.0;
          if (j != r - 1) {
            this->xnext[indj] = -phii;
            if (i != r - 1) {
              this->xnext[indi] -= phij;
              this->xnext[++ind1] = -1.0;
            }
          }
          this->xnext[npr] = -phii * phij;
          if (++ind2 >= np) {
            ind2 = 0;
          }
          this->xnext[ind2] += 1.0;
          // update Pn vy inclusion of ynext
          inclu2(ynext);
          this->xnext[ind2] = 0.0;
          if (i != r - 1) {
            // valgrind hates this
            this->xnext[indi++] = 0.0;
            this->xnext[ind1] = 0.0;
          }
        }
      }
      ithisr = nrbar - 1;
      im = this->np - 1;
      for (i = 0; i < this->np; i++) {
        bi = this->thetab[im];
        for (jm = this->np - 1, j = 0; j < i; j++) {
          bi -= this->rbar[ithisr--] * this->model.Pn[jm--];
        }
        this->model.Pn[im--] = bi;
      }
      // re-order p.
      ind = npr;
      for (i = 0; i < r; i++) this->xnext[i] = this->model.Pn[ind++];
      ind = this->np - 1;
      ind1 = npr - 1;
      for (i = 0; i < npr; i++) this->model.Pn[ind--] = this->model.Pn[ind1--];
      for (i = 0; i < r; i++) this->model.Pn[i] = this->xnext[i];
    } else {
      // P0 is obtained by back-substitution for a moving average process.
      indn = this->np;
      ind = this->np;
      for (i = 0; i < r; i++) {
        for (j = 0; j <= i; j++) {
          --ind;
          this->model.Pn[ind] = this->model.V[ind];
          if (j != 0) this->model.Pn[ind] += this->model.Pn[--indn];
        }
      }
    }
    // unpack to a full matrix
    for (i = r - 1, ind = np; i > 0; i--) {
      for (j = r - 1; j >= i; j--) {
        this->model.Pn[r * i + j] = this->model.Pn[--ind];
      }
    }
    // and symmetrize
    for (i = 0; i < r - 1; i++) {
      for (j = i + 1; j < r; j++) {
        this->model.Pn[i + r * j] = this->model.Pn[j + r * i];
      }
    }
  }
  template <class T>
  void update_arima(T &coef) {
    // copy out elements of coef into phi and theta
    for( size_t i = 0; i < p; i++) this->model.phi[i] = coef[i];
    for( size_t i = p; i < p + q; i++) this->model.theta[i-p] = coef[i];
    for (size_t i = 0; i < p; i++) this->model.T[i] = coef[i];
    if constexpr(update_P) {
      if constexpr(ss_type == SSinit::Rossignol) {
        get_Q0_rossignol();
      }
      if constexpr(ss_type == SSinit::Gardner) {
        get_Q0();
      }
    }
    // set a to all zero:
    std::fill(this->model.a.begin(), this->model.a.end(), 0);
  }
  template <typename T> std::array<scalar_t,2> arima_likelihood(
      const T &y,
      structural_model<scalar_t> &model) {

    scalar_t ssq = 0, sumlog = 0;
    size_t nu = this->n;

    for (size_t l = 0; l < n; l++) {
      for (size_t i = 0; i < r; i++) {
        scalar_t tmp = (i < r - 1) ? model.a[i + 1] : 0.0;
        if (i < p) tmp += this->model.phi[i] * this->model.a[0];
        this->anew[i] = tmp;
      }
      // template candidate
      if (d > 0) {
        for (size_t i = r + 1; i < rd; i++) this->anew[i] = this->model.a[i - 1];
        scalar_t tmp = this->model.a[0];
        for (size_t i = 0; i < d; i++) tmp += this->model.delta[i] * this->model.a[r + i];
        this->anew[r] = tmp;
      }
      if(l > 0) {
        // template candidate
        if (d == 0) {
          for (size_t i = 0; i < r; i++) {
            scalar_t vi = 0.0;
            if (i == 0) vi = 1.0; else if (i - 1 < q) vi = this->model.theta[i - 1];
            for (size_t j = 0; j < r; j++) {
              scalar_t tmp = 0.0;
              if (j == 0) tmp = vi; else if (j - 1 < q) tmp = vi * this->model.theta[j - 1];
              if (i < p && j < p) tmp += this->model.phi[i] * this->model.phi[j] * this->model.P[0];
              if (i < r - 1 && j < r - 1) tmp += model.P[i + 1 + r * (j + 1)];
              if (i < p && j < r - 1) tmp += this->model.phi[i] * this->model.P[j + 1];
              if (j < p && i < r - 1) tmp += this->model.phi[j] * this->model.P[i + 1];
              this->model.Pn[i + r * j] = tmp;
            }
          }
        } else {
          /* mm = TP */
          for (size_t i = 0; i < r; i++)
            for (size_t j = 0; j < rd; j++) {
              scalar_t tmp = 0.0;
              if (i < p) tmp += this->model.phi[i] * this->model.P[rd * j];
              if (i < r - 1) tmp += this->model.P[i + 1 + rd * j];
              this->mm[i + rd * j] = tmp;
            }
              for (size_t j = 0; j < rd; j++) {
                scalar_t tmp = model.P[rd * j];
                for (size_t k = 0; k < d; k++)
                  tmp += this->model.delta[k] * this->model.P[r + k + rd * j];
                this->mm[r + rd * j] = tmp;
                }
              for (size_t i = 1; i < d; i++)
                for (size_t j = 0; j < rd; j++)
                  this->mm[r + i + rd * j] = this->model.P[r + i - 1 + rd * j];
            /* Pnew = mmT' */
            for (size_t i = 0; i < r; i++)
              for (size_t j = 0; j < rd; j++) {
                scalar_t tmp = 0.0;
                if (i < p) tmp += this->model.phi[i] * mm[j];
                if (i < r - 1) tmp += this->mm[rd * (i + 1) + j];
                this->model.Pn[j + rd * i] = tmp;
                }
              for (size_t j = 0; j < rd; j++) {
                scalar_t tmp = mm[j];
                for (size_t k = 0; k < d; k++)
                  tmp += model.delta[k] * mm[rd * (r + k) + j];
                this->model.Pn[rd * r + j] = tmp;
                }
              for (size_t i = 1; i < d; i++)
                for (size_t j = 0; j < rd; j++)
                  model.Pn[rd * (r + i) + j] = mm[rd * (r + i - 1) + j];
            /* Pnew <- Pnew + (1 theta) %o% (1 theta) */
            for (size_t i = 0; i <= q; i++) {
              scalar_t vi = (i == 0) ? 1. : this->model.theta[i - 1];
              for (size_t j = 0; j <= q; j++)
                this->model.Pn[i + rd * j] += vi * ((j == 0) ? 1. : this->model.theta[j - 1]);
              }
            }
          }
      if (!isnan(y[l])) {
        scalar_t resid = y[l] - anew[0];
        for (size_t i = 0; i < d; i++)
          resid -= this->model.delta[i] * this->anew[r + i];
        for (size_t i = 0; i < rd; i++) {
          scalar_t tmp = this->model.Pn[i];
          for (size_t j = 0; j < d; j++)
            tmp += this->model.Pn[i + (r + j) * rd] * this->model.delta[j];
          M[i] = tmp;
          }
        scalar_t gain = M[0];
        for (size_t j = 0; j < d; j++) gain += this->model.delta[j] * M[r + j];
        if (gain < 1e4) {
          ssq += resid * resid / gain;
          sumlog += log(gain);
        } else {
          nu--;
        }
        for (size_t i = 0; i < rd; i++)
          this->model.a[i] = this->anew[i] + M[i] * resid / gain;
        for (size_t i = 0; i < rd; i++)
          for (size_t j = 0; j < rd; j++)
            this->model.P[i + j * rd] = this->model.Pn[i + j * rd] - this->M[i] * this->M[j] / gain;
        } else {
          for (size_t i = 0; i < rd; i++) this->model.a[i] = this->anew[i];
          for (size_t i = 0; i < rd * rd; i++) this->model.P[i] = this->model.Pn[i];
        }
     }
    // finally, compute likelihood and return
    const scalar_t s2 = ssq/(scalar_t)nu;
    const scalar_t loglik = 0.5 * (log(s2) + (sumlog/(scalar_t)nu));
    return std::array<scalar_t,2>{s2, loglik};
  }
};

#endif
