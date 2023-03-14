#ifndef ARIMA_ML_SOLVER
#define ARIMA_ML_SOLVER

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

template <const bool has_xreg, const bool seasonal, const bool transform>
class ARIMA_ML_PROBLEM : public cppoptlib::function::Function<double, ARIMA_ML_PROBLEM<has_xreg, seasonal, transform>> {

  using EigVec = Eigen::VectorXd;
  using EigMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
private:
  const arima_kind kind;
  const bool intercept;
  const SSinit ss_init;
  const size_t n, arma_pars;
  const std::vector<double> y;
  EigMat xreg;
  EigVec y_temp, new_x;
  std::vector<double> residual, temp_phi, temp_theta;
  structural_model<double> model;
  double sigma2;
  // add gradient approximation parameters

  // add gradient expansion functionality

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  // initialize with a given arima structure
  ARIMA_ML_PROBLEM(const std::vector<double> &y,
                   const arima_kind &kind,
                   const bool intercept,
                   const bool drift,
                   const std::vector<std::vector<double>> &xreg,
                   const double kappa,
                   const SSinit ss_init = Gardner)
    : kind(kind), intercept(intercept), ss_init(ss_init), n(y.size()),
      arma_pars(kind.p() + kind.P() + kind.q() + kind.Q()), y(y) {
    this->xreg = vec_to_mat<double>(xreg, y.size(), intercept, drift);
    this->y_temp = EigVec(n);
    for (size_t i = 0; i < n; i++) this->y_temp(i) = y[i];
    // pre-allocate new_x
    this->new_x =
      Eigen::VectorXd(kind.p() + (kind.P() * kind.period()) + kind.q() +
      (kind.Q() * kind.period()) + this->xreg.cols());
    // pre-allocate model residuals
    this->residual = std::vector<double>(n);
    // pre-allocate transformation helper vector - this is only necessary
    // for expanding seasonal models
    this->temp_phi = std::vector<double>(kind.p() + (kind.P() * kind.period()));
    this->temp_theta = std::vector<double>(kind.q() + (kind.Q() * kind.period()));
    if constexpr(seasonal || transform) {
      // the expansion of arima parameters is only necessary for seasonal models
      arima_transform_parameters<EigVec, seasonal, transform>(
          this->new_x, this->kind, this->temp_phi, this->temp_theta
      );
    }
    // // initialize state space model
    this->model = make_arima( this->new_x, this->kind, kappa, ss_init);
    this->sigma2 = 0;
  }
  double operator()(const EigVec &x) {
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
    update_arima(this->model, this->new_x, this->kind, this->ss_init);
    // return likelihood - check if this updated model or not (it ideally
    // should not, not here)
    const std::array<double,2> res = arima_likelihood(this->y_temp, this->model);
    this->sigma2 = res[0];
    return res[1];
  }
  // add impl of grad, hessian, eval
  void Gradient(const EigVec &x, EigVec *grad) {
    cppoptlib::utils::ComputeFiniteGradient(*this, x, grad);
  }
  // Computes the Hessian of a function.
  void Hessian(const EigVec &x, EigMat *hessian) {
    cppoptlib::utils::ComputeFiniteHessian(*this, x, hessian);
  }
  // For improved performance, this function will return the state directly.
  // Override this method if you can compute the objective value, gradient and
  // Hessian simultaneously.
  cppoptlib::function::State<double, EigVec, EigMat> Eval(const EigVec &x,
                                            const int order = 1) { //const
    cppoptlib::function::State<double, EigVec, EigMat> state(x.rows(), order);
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
  void finalize() {
    // recompute  and replace V as this has not been updated
    // since initialization in make_arima
    recompute_v(this->model);
  }
  double get_sigma() const { return this->sigma2; }
  structural_model<double> get_structural_model() const { return this->model; }
};

template <const bool has_xreg, const bool seasonal, const bool transform>
double arima_solver_ml(std::vector<double> &y,
                     structural_model<double> &model,
                     const bool intercept,
                     const bool drift,
                     std::vector<std::vector<double>> xreg,
                     const arima_kind &kind,
                     std::vector<double> &coef,
                     const double kappa,
                     const SSinit ss_init) {
  Eigen::VectorXd x(coef.size());
  for (size_t i = 0; i < coef.size(); i++) x(i) = coef[i];
  // initialize solver
  cppoptlib::solver::Bfgs<
    ARIMA_ML_PROBLEM<has_xreg, seasonal, transform>> solver;
  // and arima problem
  ARIMA_ML_PROBLEM<has_xreg, seasonal, transform> ml_arima_problem(
      y, kind, intercept, drift, xreg, kappa, ss_init);
  // and finally, minimize
  auto [solution, solver_state] = solver.Minimize(ml_arima_problem, x);
  // replace model V
  ml_arima_problem.finalize();
  // write back structural model
  model = ml_arima_problem.get_structural_model();
  for (size_t i = 0; i < coef.size(); i++) coef[i] = solution.x[i];
  // return variance estimate for the arima model
  return ml_arima_problem.get_sigma();
}

#endif
