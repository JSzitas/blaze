#ifndef ARIMA_ML_KALMAN
#define ARIMA_ML_KALMAN

#include "arima/structures/structural_model.h"
#include "arima/structures/arima_kind.h"
#include "arima/structures/ss_init.h"

#include "arima/utils/transforms.h"
#include "arima/utils/delta.h"

#include "third_party/eigen.h"
#include "utils/utils.h"

#include "arima/solvers/arima_css_likelihood.h"
#include "arima/solvers/state_space.h"


// main class for doing Kalman filtering of an ARIMA model

template <const SSinit ss_type, const bool seasonal,
          const bool has_xreg, const bool transform,
          typename scalar_t=double> class KalmanARIMA {

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
vec xnext, xrow, rbar, thetab, V, P;

EigVec y_temp, new_x;
vec residual, temp_phi, temp_theta;

structural_model<scalar_t> model;
scalar_t sigma2;

public:
  KalmanARIMA<ss_type, seasonal, has_xreg, transform, scalar_t>(){};
  KalmanARIMA<ss_type, seasonal, has_xreg, transform, scalar_t>(
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
      this->V = vec(this->np);
      this->P = vec(r * r);
    }
    this->y_temp = EigVec(this->n);
    for (size_t i = 0; i < this->n; i++) this->y_temp(i) = y[i];
    // pre-allocate new_x
    this->new_x = EigVec::Zero(kind.p() + (kind.P() * kind.period()) + kind.q() +
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
    // initialize state space model
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
      this->y_temp -= this->xreg * x.tail(x.size() - this->arma_pars);
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
    // this->update_arima();
    if constexpr( ss_type == SSinit::Gardner ) {
      update_arima(this->model, this->new_x, this->kind,
                   this->xnext, this->xrow, this->rbar,
                   this->thetab, this->V, this->P,
                   SSinit::Gardner);
    }
    if constexpr( ss_type == SSinit::Rossignol ) {
      update_arima(this->model, this->new_x, this->kind, SSinit::Rossignol);
    }
    // return likelihood - check if this updated model or not (it ideally
    // should not, not here)
    const std::array<scalar_t,2> res = arima_likelihood(
      this->y_temp, this->model, this->anew, this->M, this->mm);
    this->sigma2 = res[0];
    return res[1];
  }
  scalar_t get_sigma() const { return this->sigma2; }
  structural_model<scalar_t> get_structural_model() const { return this->model; }
};

#endif
