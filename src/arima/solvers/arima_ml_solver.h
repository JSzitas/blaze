#ifndef ARIMA_CSS_SOLVER
#define ARIMA_CSS_SOLVER

#include "utils/utils.h"

#include "arima/structures/structural_model.h"
#include "arima/structures/fitting_method.h"

#include "arima/solvers/arima_css_likelihood.h"
#include "arima/solvers/state_space.h"

#include "arima/utils/transforms.h"
#include "arima/utils/xreg.h"

// include optimizer library
#include "third_party/eigen.h"
#include "third_party/optim.h"

using FunctionXd = cppoptlib::function::Function<double>;

template <const bool has_xreg, const bool seasonal>
class ARIMA_ML_PROBLEM : public FunctionXd {
private:
  arima_kind kind;
  lm_coef<double> xreg_pars;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> xreg;
  size_t arma_pars;
  std::vector<double> y;
  Eigen::VectorXd y_temp;
  size_t n;
  Eigen::VectorXd new_x;
  std::vector<double> residual;
  std::vector<double> transform_temp_phi;
  std::vector<double> transform_temp_theta;
  structural_model<double> &model;
  const SSinit ss_init;
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  // initialize with a given arima structure
  ARIMA_ML_PROBLEM(std::vector<double> &y,
                   const arima_kind &kind,
                   lm_coef<double> &xreg_pars,
                   std::vector<std::vector<double>> &xreg,
                   structural_model<double> &model,
                   std::vector<double> & delta,
                   double kappa,
                   SSinit ss_init = Gardner)
    : kind(kind), xreg_pars(xreg_pars), ss_init(ss_init) {
    // initialize an xreg matrix
    size_t n = y.size();
    std::vector<double> _xreg = flatten_vec(xreg);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> new_mat =
      Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>(
          _xreg.data(), n, xreg.size());
    if (xreg_pars.has_intercept()) {
      new_mat.conservativeResize(Eigen::NoChange, new_mat.cols() + 1);
      new_mat.col(new_mat.cols() - 1) =
        Eigen::Matrix<double, Eigen::Dynamic, 1>::Constant(n, 1, 1);
    }
    this->xreg = new_mat;
    this->arma_pars = kind.p() + kind.P() + kind.q() + kind.Q();
    Eigen::VectorXd y_temp(n);
    for (size_t i = 0; i < n; i++) {
      y_temp[i] = y[i];
    }
    this->y = y;
    this->y_temp = y_temp;
    this->n = n;
    this->xreg = new_mat;
    this->arma_pars = kind.p() + kind.P() + kind.q() + kind.Q();
    // pre-allocate new_x
    this->new_x =
      Eigen::VectorXd(kind.p() + (kind.P() * kind.period()) + kind.q() +
      (kind.Q() * kind.period()) + this->xreg.cols());
    // pre-allocate model residuals
    this->residual = std::vector<double>(n);
    // pre-allocate transformation helper vector - this is only necessary
    // for expanding seasonal models
    this->transform_temp_phi =
      std::vector<double>(kind.p() + (kind.P() * kind.period()));
    this->transform_temp_theta =
      std::vector<double>(kind.q() + (kind.Q() * kind.period()));
    this->model = model;
    if constexpr(seasonal) {
      // the expansion of arima parameters is only necessary for seasonal models
      arima_transform_parameters<seasonal, false>(this->new_x, this->kind,
                                                  this->transform_temp_phi,
                                                  this->transform_temp_theta);
    }
    // initialize state space model
    structural_model<double> arima_ss = make_arima( this->new_x,
                                                    delta, this->kind,
                                                    kappa, ss_init);
    this->model.set(arima_ss);
  }
  double operator()(const Eigen::VectorXd &x) {
    for (int i = 0; i < x.size(); i++) {
      this->new_x[i] = x[i];
    }
    if constexpr (has_xreg) {
      // refresh y_temp and load it with original y data
      for (size_t i = 0; i < this->n; i++) {
        this->y_temp[i] = this->y[i];
      }
      this->y_temp =
        this->y_temp - this->xreg * x.tail(x.size() - this->arma_pars);
    }
    /* I figured out that I can basically expand this out altogether for non-seasonal
     * models - the compiler should insert an empty function anyways, but just to
     * make sure that this gets compiled away - we can make sure its a dead branch
     */
    if constexpr(seasonal) {
      // the expansion of arima parameters is only necessary for seasonal models
      arima_transform_parameters<seasonal, false>(this->new_x, this->kind,
                                                  this->transform_temp_phi,
                                                  this->transform_temp_theta);
    }
    // update arima
    update_arima(this->model, this->new_x, this->ss_init);
    // return likelihood - check if this updated model or not (it ideally
    // should not, not here)
    return arima_likelihood(this->y_temp, this->model);
  }
};

template <const bool has_xreg, const bool seasonal>
void arima_solver_ml(std::vector<double> &y,
                     structural_model<double> &model,
                     lm_coef<double> xreg_coef,
                     std::vector<std::vector<double>> xreg,
                     const arima_kind &kind,
                     std::vector<double> &coef,
                     std::vector<double> &delta, const double kappa,
                     const SSinit ss_init, double &sigma2) {

  auto vec_size = kind.p() + kind.q() + kind.P() + kind.Q() + xreg_coef.size();
  auto arma_size = kind.p() + kind.q() + kind.P() + kind.Q();
  Eigen::VectorXd x(vec_size);
  for (auto &val : x) {
    val = 0;
  }
  // initialize to all zeroes except for xreg
  for (size_t i = arma_size; i < vec_size; i++) {
    x[i] = xreg_coef.coef[i - arma_size];
  }
  // initialize solver
  using Solver = cppoptlib::solver::Bfgs<ARIMA_ML_PROBLEM<has_xreg, seasonal>>;
  Solver solver;
  // and arima problem
  ARIMA_ML_PROBLEM<has_xreg, seasonal> ml_arima_problem(y, kind, xreg_coef,
                                                          xreg, model,
                                                          delta, kappa,
                                                          ss_init);
  // and finally, minimize
  auto [solution, solver_state] = solver.Minimize(ml_arima_problem, x);
  // update variance estimate for the arima model - this was passed by reference
  sigma2 = exp(2 * solution.value);
  // pass fitted coefficients back to the caller
  for (size_t i = 0; i < vec_size; i++) {
    coef[i] = solution.x[i];
  }
}

#endif
