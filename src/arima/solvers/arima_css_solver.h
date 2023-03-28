#ifndef ARIMA_CSS_SOLVER
#define ARIMA_CSS_SOLVER

#include "utils/utils.h"

#include "arima/structures/structural_model.h"
#include "arima/structures/fitting_method.h"
#include "arima/structures/arima_kind.h"

#include "arima/solvers/arima_css_likelihood.h"
#include "arima/solvers/state_space.h"

#include "arima/utils/transforms.h"
#include "utils/xreg.h"

// include optimizer library
#include "third_party/eigen.h"
#include "third_party/optim.h"

template <const bool has_xreg, const bool seasonal, typename scalar_t=float>
class ARIMA_CSS_PROBLEM : public cppoptlib::function::Function<scalar_t, ARIMA_CSS_PROBLEM<has_xreg, seasonal, scalar_t>> {

  using EigVec = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;
  using EigMat = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;
  using StateXd = cppoptlib::function::State<scalar_t, EigVec, EigMat>;

private:
  const arima_kind kind;
  const bool intercept;
  const size_t n_cond, n;
  const std::vector<scalar_t> y;
  EigMat xreg;
  size_t arma_pars;
  EigVec y_temp, new_x;
  std::vector<scalar_t> residual, temp_phi, temp_theta;
public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  // initialize with a given arima structure
  ARIMA_CSS_PROBLEM(const std::vector<scalar_t> &y,
                    const arima_kind &kind,
                    const bool intercept,
                    const bool drift,
                    std::vector<std::vector<scalar_t>> &xreg)
      : kind(kind), intercept(intercept),
        n_cond(kind.d() + (kind.D() * kind.period()) +
        kind.p() + (kind.P() * kind.period())), n(y.size()), y(y) {
    this->xreg = vec_to_mat(xreg, y.size(), intercept, drift);
    this->arma_pars = kind.p() + kind.P() + kind.q() + kind.Q();
    this->y_temp = EigVec(n);
    for (size_t i = 0; i < n; i++) this->y_temp(i) = y[i];
    if constexpr (!has_xreg) {
      // if you have no intercept, you can do differencing
      // regular differencing
      for (size_t i = 0; i < kind.d(); i++) {
        for (size_t l = n - 1; l > 0; l--) {
          this->y_temp(l) -= this->y_temp(l - 1);
        }
      }
      // seasonal differencing
      size_t ns = kind.period();
      for (size_t i = 0; i < kind.D(); i++) {
        for (size_t l = n - 1; l >= ns; l--) {
          this->y_temp(l) -= this->y_temp(l - ns);
        }
      }
    }
    this->arma_pars = kind.p() + kind.P() + kind.q() + kind.Q();
    // pre-allocate new_x
    this->new_x = EigVec(kind.p() + (kind.P() * kind.period()) +
                         kind.q() + (kind.Q() * kind.period()) +
                         this->xreg.cols());
    // pre-allocate model residuals
    this->residual = std::vector<scalar_t>(n, 0.0);
    // pre-allocate transformation helper vector - this is only necessary
    // for expanding seasonal models
    this->temp_phi = std::vector<scalar_t>(kind.p() + (kind.P() * kind.period()));
    this->temp_theta = std::vector<scalar_t>(kind.q() + (kind.Q() * kind.period()));
  }
  scalar_t operator()(const EigVec &x) {
    if constexpr (has_xreg) {
      // refresh y_temp and load it with original y data
      for (size_t i = 0; i < this->n; i++) this->y_temp(i) = this->y[i];
      this->y_temp -= this->xreg * x.tail(x.size() - this->arma_pars);
      // do differencing here
      if (this->intercept) {
        for (size_t i = 0; i < kind.d(); i++) {
          for (size_t l = n - 1; l > 0; l--) {
            this->y_temp(l) -= this->y_temp(l - 1);
          }
        }
        // seasonal differencing
        size_t ns = kind.period();
        for (size_t i = 0; i < kind.D(); i++) {
          for (size_t l = n - 1; l >= ns; l--) {
            this->y_temp(l) -= this->y_temp(l - ns);
          }
        }
      }
    }
    /* I figured out that I can basically expand this out altogether for non-seasonal
     * models - the compiler should insert an empty function anyways, but just to
     * make sure that this gets compiled away - we can make sure its a dead branch
     */
    // we can expand these in reverse order for simd - which should in principle allow
    //  us to not have to reverse in simd css ssq
    for (size_t i = 0; i < x.size(); i++) this->new_x(i) = x(i);
    if constexpr(seasonal) {
      // the expansion of arima parameters is only necessary for seasonal models
      arima_transform_parameters<EigVec, seasonal, false, true, scalar_t>(
          this->new_x,
          this->kind,
          this->temp_phi,
          this->temp_theta
      );
    }
    scalar_t res = simd_arima_css_ssq<EigVec, scalar_t>(
      this->y_temp, this->new_x, this->kind, this->n_cond, this->residual);
    return 0.5 * log(res);
  }
  // add impl of grad, hessian, eval
  void Gradient(const EigVec &x, EigVec *grad) {
    cppoptlib::utils::ComputeFiniteGradient(*this, x, grad);
  }
  void Hessian(const EigVec &x, EigMat *hessian) {
    cppoptlib::utils::ComputeFiniteHessian(*this, x, hessian);
  }
  StateXd Eval(const EigVec &x, const int order = 1) { //const
    StateXd state(x.rows(), order);
    state.value = this->operator()(x);
    state.x = x;
    if (order >= 1) {
      this->Gradient(x, &state.gradient);
    }
    if ((order >= 2) && (state.hessian)) {
      this->Hessian(x, &*(state.hessian));
    }
    return state;
  }
  void finalize( structural_model<scalar_t> &model,
                 const EigVec & final_pars,
                 scalar_t kappa,
                 SSinit ss_init) {
    // this function creates state space representation and expands it
    // I found out it is easier and cheaper (computationally) to do here
    // do the same for model coefficients
    // finally, make state space model

    // fix reexpansion for reversed coef vector
    for (size_t i = 0; i < final_pars.size(); i++) this->new_x(i) = final_pars(i);
    if constexpr(seasonal) {
      // the expansion of arima parameters is only necessary for seasonal models
      arima_transform_parameters<EigVec, seasonal, false, false, scalar_t>(
          this->new_x,
          this->kind,
          this->temp_phi,
          this->temp_theta
      );
    }

    structural_model<scalar_t> arima_ss = make_arima<EigVec, scalar_t>(
      this->new_x, this->kind, kappa, ss_init);

    model.set(arima_ss);
    // modify y_temp to acount for xreg
    for (size_t i = 0; i < this->n; i++) {
      this->y_temp[i] = this->y[i];
    }
    if constexpr(has_xreg) {
      this->y_temp = this->y_temp -
        this->xreg * final_pars.tail(final_pars.size() - this->arma_pars);
    }
    // get arima steady state values
    arima_steady_state(this->y_temp, model);
  }
};

template <const bool has_xreg, const bool seasonal, typename scalar_t = float>
scalar_t arima_solver_css(std::vector<scalar_t> &y,
                        const arima_kind &kind,
                        structural_model<scalar_t> &model,
                        std::vector<std::vector<scalar_t>> xreg,
                        const bool intercept,
                        const bool drift,
                        std::vector<scalar_t> &coef,
                        const scalar_t kappa,
                        const SSinit ss_init) {
  Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> x(coef.size());
  for(size_t i = 0; i < coef.size(); i++) x(i) = coef[i];
  // initialize solver
  cppoptlib::solver::Bfgs<
    ARIMA_CSS_PROBLEM<has_xreg, seasonal, scalar_t>> solver;
  // and arima problem
  ARIMA_CSS_PROBLEM<has_xreg, seasonal, scalar_t> css_arima_problem(
      y, kind, intercept, drift, xreg);
  // and finally, minimize
  auto [solution, solver_state] = solver.Minimize(css_arima_problem, x);
  // update variance estimate for the arima model - this was passed by reference
  css_arima_problem.finalize( model, solution.x, kappa, ss_init);
  // pass fitted coefficients back to the caller
  for (size_t i = 0; i < coef.size(); i++) coef[i] = solution.x[i];
  return exp(2 * solution.value);
}

#endif

