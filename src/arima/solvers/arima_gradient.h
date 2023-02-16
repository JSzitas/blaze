#ifndef ARIMA_CSS_GRAD
#define ARIMA_CSS_GRAD

#include "arima/structures/arima_kind.h"

#include "utils/utils.h"

#include <cmath>
#include <math.h>

#include "third_party/eigen.h"

/* Computing the (S)ARIMA(X) gradient can be costly as many operations
 * have to be repeated even when they are potentially unnecessary (e.g.
 * handling of xreg, parameter expansion). To combat this, when computing
 * gradients, we can actually leverage this knowledge, and only run through a
 * subset of the objective function relevant to that specific parameter.
 * To do that, we use this custom class (since this will require some temporary
 * objects).
 */
template <const bool seasonal, const bool has_xreg, typename scalar_t = double> class ArimaLossGradient {
  // convenience typedefs
  using EigVec = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;
  using EigMat = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;

  std::vector<scalar_t> y;
  size_t n, n_cond, accuracy,
  // arima kind size definitions
  p, q, d, D, ns, mp, mq, msp, msq, arma_pars;
  bool has_intercept;
  EigMat xreg;
  EigVec phi, theta, resid, y_temp;
  // finite difference approximation parameters
  size_t inner_steps;
  scalar_t dd_val;
  // taken from cppNumericalSolvers, see:
  // https://github.com/PatWie/CppNumericalSolvers/blob/v2/include/cppoptlib/utils/derivatives.h
  // i.e. the same library we use for BFGS
  // Many thanks - very useful to not have to rethink this computation (and find a decent eps)
  constexpr static scalar_t eps = 2.2204e-6;
  std::array<std::vector<scalar_t>, 4> coeff;
  std::array<std::vector<scalar_t>, 4> coeff2;
  std::array<scalar_t, 4> dd;
public:
  ArimaLossGradient<seasonal, has_xreg, scalar_t>(){};
  ArimaLossGradient<seasonal, has_xreg, scalar_t>(
      const std::vector<scalar_t> &y,
      const arima_kind &kind,
      const bool has_intercept,
      EigMat xreg,
      const size_t n_cond = 0,
      const size_t accuracy = 0) {

    this->y = y; this->n = y.size();
    this->n_cond = n_cond;
    this->accuracy = accuracy;
    // basically just decompose arima_kind and compute some intermediary values
    this->p = kind.p() + (kind.period() * kind.P());
    this->q = kind.q() + (kind.period() * kind.Q());
    this->d = kind.d(); this->D = kind.D(); this->ns = kind.period();
    this->mp = kind.p(); this->mq = kind.q(); this->msp = kind.P();
    this->msq = kind.Q();
    // the size of the coefficient vector that will be passed to us
    // if you ignore xreg
    this->arma_pars = kind.p() + kind.P() + kind.q() + kind.Q();
    this->has_intercept = has_intercept;
    // the xreg fixed matrix
    this->xreg = xreg;
    // helper vectors for holding and updating parts of the model
    // only used for seasonal models
    if constexpr(seasonal) {
      this->phi = EigVec(this->p);
      this->theta = EigVec(this->q);
    }
    // residuals from running the ar part of the model
    this->resid = EigVec(this->n).setZero();
    // temporary y values
    this->y_temp = EigVec(this->n);
    for (size_t i = 0; i < n; i++) this->y_temp(i) = this->y[i];
    // inner gradient approximation steps
    this->inner_steps = 2 * (this->accuracy + 1);

    // only if you have no xreg can you do differencing here
    if constexpr(!has_xreg) this->diff_y_temp();
    this->coeff = { {{1, -1},
                     {1, -8, 8, -1},
                     {-1, 9, -45, 45, -9, 1},
                     {3, -32, 168, -672, 672, -168, 32, -3}}
      };
    this->coeff2 = { {{1, -1},
                      {-2, -1, 1, 2},
                      {-3, -2, -1, 1, 2, 3},
                      {-4, -3, -2, -1, 1, 2, 3, 4}}
      };
    this->dd = {2, 12, 60, 840};
    this->dd_val = dd[this->accuracy] * this->eps;
  }
  scalar_t loss(const EigVec & x) {
    // complete loss computation
    // apply xreg
    this->apply_xreg(x);
    // expand phi and theta
    if constexpr(seasonal) {
      this->expand_phi(x);
      this->expand_theta(x);
    }
    // update ar part of the model
    // this->update_ar(x);
    // update theta part of the model and return loss
    // return this->update_ma_loss(x);
    return this->arma_loss(x);
  }
  EigVec Gradient(const EigVec & x) {
    // get rid of const - the interface is const, and we do not actually
    // modify the underlying parameters at the end, but it is more efficient
    // for us to modify them here

    EigVec grad(x.size());
    grad.setZero();
    EigVec &pars = const_cast<EigVec &>(x);
    // carry out xreg evaluation once at the start since we are iterating from ar coefficients
    this->apply_xreg(pars);
    // expand phi and update arpart of model
    if constexpr(seasonal) {
      this->expand_phi(pars);
      this->update_ar(pars);
    }
    /* run operations based on current iteration block - i.e. it is kind of
     * smart to run them from back to front:
     * MA operations first => only require MA updates
     * AR operations second => require AR AND MA updates
     * xreg => basically the only thing that requires the whole objective function
     * however MA does not change the state of the gradient class whatsoever
     * whereas AR only updates ar_res and xreg requires evaluation of everything.
     * hence we start with MA terms
     */
     for( size_t i = this->mp; i < this->mp + this->mq; i++ ) {
       for (size_t s = 0; s < this->inner_steps; ++s) {
         scalar_t tmp = pars[i];
         // update eps value
         pars[i] += this->coeff2[accuracy][s] * this->eps;
         // expand theta with new parameter value
         if constexpr(seasonal) this->expand_theta(pars);
         grad[i] += this->coeff[accuracy][s] * arma_loss(pars);
         pars[i] = tmp;
       }
     }
     // expand out seasonal MA
     for( size_t i = this->mp + this->mq + this->msp; i < arma_pars; i++) {
       for (size_t s = 0; s < this->inner_steps; ++s) {
         scalar_t tmp = pars[i];
         // update eps value
         pars[i] += this->coeff2[accuracy][s] * this->eps;
         // expand theta with new parameter value
         if constexpr(seasonal) this->expand_theta(pars);
         grad[i] += this->coeff[accuracy][s] * this->arma_loss(pars);
         pars[i] = tmp;
       }
     }
     // now, update AR part of the model
     for( size_t i = 0; i < this->mp; i++ ) {
       for (size_t s = 0; s < this->inner_steps; ++s) {
         scalar_t tmp = pars[i];
         // update eps value
         pars[i] += this->coeff2[accuracy][s] * this->eps;
         // expand both phi and theta with new parameter value
         if constexpr(seasonal) this->expand_phi(pars);
         grad[i] += this->coeff[accuracy][s] * this->arma_loss(pars);
         pars[i] = tmp;
       }
     }
     // procceed to seasonal AR
     for( size_t i = this->mp + this->mq;
          i < this->mp + this->mq + this->msp; i++ ) {
       for (size_t s = 0; s < this->inner_steps; ++s) {
         scalar_t tmp = pars[i];
         // update eps value
         pars[i] += this->coeff2[accuracy][s] * this->eps;
         // expand phi with new parameter value
         // note that since theta does not change, we do not need to expand it
         if constexpr(seasonal) this->expand_phi(pars);
         // update all arma coefficients
         grad[i] += this->coeff[accuracy][s] * this->arma_loss(pars);
         pars[i] = tmp;
       }
     }
     // update both expansions to unchanged coefficients
     if constexpr(seasonal) {
       this->expand_phi(pars);
       this->expand_theta(pars);
     }
     // finally, get to the xreg part of the model
     for( size_t i = this->arma_pars; i < pars.size(); i++ ) {
       for (size_t s = 0; s < this->inner_steps; ++s) {
         scalar_t tmp = pars[i];
         // update eps value
         pars[i] += this->coeff2[accuracy][s] * this->eps;
         // apply xreg and refilter the series - note that there is no need
         // to expand
         this->apply_xreg(pars);
         grad[i] += this->coeff[accuracy][s] * this->arma_loss(pars);
         pars[i] = tmp;
       }
     }
     // rescale gradient
     for( auto &val : grad) {
       val /= dd_val;
     }
     return grad;
  }
  const EigVec get_y_temp(const EigVec & x) {
    // modify y_temp to acount for xreg
    for (size_t i = 0; i < this->n; i++) this->y_temp[i] = this->y[i];
    if constexpr(has_xreg) {
      this->y_temp -= this->xreg * x.tail(x.size() - this->arma_pars);
    }
    return this->y_temp;
  }
  const EigVec get_expanded_coef(const EigVec & x) {
    if constexpr( seasonal ) {
      this->expand_phi(x);
      this->expand_theta(x);
      EigVec res(this->phi.size() + this->theta.size());
      for(size_t i = 0; i < phi.size(); i++) {
        res[i] = this->phi[i];
      }
      for(size_t i = this->phi.size(); i < this->phi.size() + this->theta.size(); i++ ) {
        res[i] = this->theta[i - this->phi.size()];
      }
      return res;
    }
    // otherwise no expansion is necessary
    return x;
  }
private:
  // for this we need a vector that exists for after applying AR to the
  // current series, since we can run updates to the MA part using just that
  void update_ar(const EigVec & x) {
    for (size_t l = this->n_cond; l < this->n; l++) {
      // notice how the p parameters only impact tmp at that point and do not
      // interact with the q parameters - this means the likelihood (or CSS)
      // is separable here - so for updates of q parameters, if you have pre-saved
      // values of this tmp somewhere, you can just fetch them :)
      scalar_t tmp = y_temp(l);
      for (size_t j = 0; j < this->p; j++) {
        // for seasonal models we had to expand
        if constexpr( seasonal ) {
          tmp -= this->phi(j) * y_temp(l - j - 1);
        }
        // for non-seasonal models we can just use unexpanded coefficients
        if constexpr( !seasonal ) {
          tmp -= x(j) * y_temp(l - j - 1);
        }
      }
      this->resid(l) = tmp;
    }
  }
  scalar_t update_ma_loss(const EigVec & x) {
    // this will be responsible for updating the residuals from ma coefficients
    size_t nu = this->n - this->n_cond;
    scalar_t ssq = 0.0;
    for (size_t l = this->n_cond; l < this->n; l++) {
      const size_t ma_offset = min(l - this->n_cond, this->q);
      scalar_t tmp = this->resid(l);
      for (size_t j = 0; j < ma_offset; j++) {
        // see above for why this works
        if constexpr(seasonal) {
          tmp -= this->theta(j) * this->resid(l - j - 1);
        }
        if constexpr(!seasonal) {
          tmp -= x(this->p + j) * this->resid(l - j - 1);
        }
      }
      this->resid(l) = tmp;
      if (!isnan(tmp)) {
        ssq += tmp * tmp;
      }
      else {
        nu--;
      }
    }
    return 0.5 * log(ssq / nu);
  }
  scalar_t arma_loss(const EigVec & x) {
    // if no nans are met, this nu does not need to be changed
    size_t nu = this->n - this->n_cond;
    scalar_t ssq = 0.0;
    for (size_t l = this->n_cond; l < this->n; l++) {
      const size_t ma_offset = min(l - this->n_cond, this->q);
      scalar_t tmp = y_temp(l);
      for (size_t j = 0; j < this->p; j++) {
        // for seasonal models we had to expand
        if constexpr( seasonal ) {
          tmp -= this->phi(j) * y_temp(l - j - 1);
        }
        // for non-seasonal models we can just use unexpanded coefficients
        if constexpr( !seasonal ) {
          tmp -= x(j) * y_temp(l - j - 1);
        }
      }
      for (size_t j = 0; j < ma_offset; j++) {
        // see above for why this works
        if constexpr(seasonal) {
          tmp -= this->theta(j) * this->resid(l - j - 1);
        }
        if constexpr(!seasonal) {
          tmp -= x(this->p + j) * this->resid(l - j - 1);
        }
      }
      this->resid(l) = tmp;
      if (!isnan(tmp)) {
        ssq += tmp * tmp;
      }
      // potential optimization - if this branch is reached,
      // (which should be rare) decrement the pre-allocated nu
      else{
        nu--;
      }
    }
    return 0.5 * log(ssq / nu);
  }
  void apply_xreg(const EigVec & x) {
    if constexpr(has_xreg) {
      // refresh y_temp and load it with original y data
      for (size_t i = 0; i < this->n; i++) this->y_temp(i) = this->y[i];
      this->y_temp = this->y_temp - this->xreg * x.tail(x.size() - this->arma_pars);
      // differencing only needs to be applied after we reran xreg
      // and only if we do not have an intercept
      if( !this->has_intercept ) this->diff_y_temp();
    }
  }
  // these expansions do not operate on coef but use the temporaries
  // stored within the class (i.e. phi and theta)
  void expand_phi(const EigVec & coef) {
    // non-seasonal models are no-ops
    if constexpr(seasonal) {
      const size_t phi_offset = this->mp + this->mq;
      for (size_t i = 0; i < this->mp; i++) this->phi(i) = coef(i);
      for (size_t i = this->mp; i < this->p; i++) this->phi(i) = 0.0;
      for (size_t j = 0; j < this->msp; j++) {
        this->phi((j + 1) * this->ns - 1) += coef(phi_offset + j);
        for (size_t i = 0; i < this->mp; i++) {
          this->phi((j + 1) * this->ns + i) -= coef(i) * coef(phi_offset + j);
        }
      }
    }
  }
  void expand_theta(const EigVec & coef) {
    if constexpr(seasonal) {
      const size_t theta_offset = this->mp + this->mq + this->msp;
      for (size_t i = 0; i < this->mq; i++) this->theta(i) = coef(i + this->mp);
      for (size_t i = this->mq; i < this->q; i++) theta(i) = 0.0;
      for (size_t j = 0; j < this->msq; j++) {
        this->theta((j + 1) * ns - 1) += coef(theta_offset + j);
        for (size_t i = 0; i < this->mq; i++) {
          this->theta((j + 1) * this->ns + i) += coef(i + this->mp) * coef(theta_offset + j);
        }
      }
    }
  }
  void diff_y_temp() {
    // regular differencing
    for (size_t i = 0; i < this->d; i++) {
      for (size_t l = n - 1; l > 0; l--) {
        this->y_temp(l) -= this->y_temp(l - 1);
      }
    }
    // seasonal differencing
    for (size_t i = 0; i < this->D; i++) {
      for (size_t l = n - 1; l >= this->ns; l--) {
        this->y_temp(l) -= this->y_temp(l - ns);
      }
    }
  }
};



#endif
