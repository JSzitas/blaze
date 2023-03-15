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

template <const bool has_xreg, const bool seasonal>
class ARIMA_CSS_PROBLEM : public cppoptlib::function::Function<double, ARIMA_CSS_PROBLEM<has_xreg, seasonal>> {

  using EigVec = Eigen::VectorXd;
  using EigMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
  using StateXd = cppoptlib::function::State<double, Eigen::VectorXd, Eigen::MatrixXd>;

private:
  const arima_kind kind;
  const bool intercept;
  const size_t n_cond, n;
  const std::vector<double> y;
  EigMat xreg;
  size_t arma_pars;
  EigVec y_temp, new_x;
  std::vector<double> residual, temp_phi, temp_theta;
public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  // initialize with a given arima structure
  ARIMA_CSS_PROBLEM(const std::vector<double> &y,
                    const arima_kind &kind,
                    const bool intercept,
                    const bool drift,
                    std::vector<std::vector<double>> &xreg)
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
    this->residual = std::vector<double>(n, 0.0);
    // pre-allocate transformation helper vector - this is only necessary
    // for expanding seasonal models
    this->temp_phi = std::vector<double>(kind.p() + (kind.P() * kind.period()));
    this->temp_theta = std::vector<double>(kind.q() + (kind.Q() * kind.period()));
  }
  double operator()(const EigVec &x) {
    for (size_t i = 0; i < x.size(); i++) this->new_x(i) = x(i);
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
    if constexpr(seasonal) {
      // the expansion of arima parameters is only necessary for seasonal models
      arima_transform_parameters<EigVec, seasonal, false>(
          this->new_x,
          this->kind,
          this->temp_phi,
          this->temp_theta
      );
    }
    double res = arima_css_ssq(this->y_temp, this->new_x, this->kind,
                               this->n_cond, this->residual);
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
  void finalize( structural_model<double> &model,
                 const EigVec & final_pars,
                 double kappa,
                 SSinit ss_init) {
    // this function creates state space representation and expands it
    // I found out it is easier and cheaper (computationally) to do here
    // do the same for model coefficients
    // finally, make state space model
    structural_model<double> arima_ss = make_arima(
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

template <const bool has_xreg, const bool seasonal>
double arima_solver_css(std::vector<double> &y,
                        const arima_kind &kind,
                        structural_model<double> &model,
                        std::vector<std::vector<double>> xreg,
                        const bool intercept,
                        const bool drift,
                        std::vector<double> &coef,
                        const double kappa,
                        const SSinit ss_init) {
  Eigen::VectorXd x(coef.size());
  for(size_t i = 0; i < coef.size(); i++) x(i) = coef[i];
  // initialize solver
  cppoptlib::solver::Bfgs<ARIMA_CSS_PROBLEM<has_xreg, seasonal>> solver;
  // and arima problem
  ARIMA_CSS_PROBLEM<has_xreg, seasonal> css_arima_problem(
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

