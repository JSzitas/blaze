#ifndef ARIMA_LASSO_SOLVER
#define ARIMA_LASSO_SOLVER

#include "utils/utils.h"
#include "utils/xreg.h"

#include "arima/structures/structures.h"
#include "arima/solvers/state_space.h"

#include "arima/utils/transforms.h"

// include optimizer library
#include "third_party/eigen.h"
#include "third_party/optim.h"

template <typename C, typename scalar_t> scalar_t sum_abs(const C &x) {
  scalar_t res = 0;
  for( auto& val:x ) res += std::abs(val);
  return res;
}

template <typename C, typename scalar_t> scalar_t sum_abs(
    const C &x, const scalar_t from, const scalar_t to) {
  scalar_t res = 0;
  for(size_t i = from; i < to; i++) res += std::abs(x[i]);
  return res;
}

template <typename C, typename scalar_t,
          const bool arima_only=false> scalar_t arima_css_ssq_lasso(
    const C & y, const C & pars, const arima_kind kind,
    const scalar_t penalty,
    const size_t n_cond, std::vector<scalar_t> resid) {
  const size_t n = y.size(), p = kind.p() + kind.period() * kind.P(),
    q = kind.q() + kind.period() * kind.Q();
  // prepare the residuals - possibly move this out and never allocate here?
  int ma_offset, nu = n-n_cond;
  scalar_t ssq = 0.0, tmp = 0.0;
  for (size_t l = n_cond; l < n; l++) {
    ma_offset = min(l - n_cond, q);
    tmp = y[l];
    for (size_t j = 0; j < p; j++) {
      tmp -= pars[j] * y[l - j - 1];
    }
    // to offset that this is all in one vector, we need to
    // start at p and go to p + q
    for (size_t j = 0; j < ma_offset; j++) {
      tmp -= pars[p + j] * resid[l - j - 1];
    }
    resid[l] = tmp;
    if (!isnan(tmp)) {
      ssq += tmp * tmp;
    }
    else {
      nu--;
    }
  }
  if constexpr(arima_only) {
    return (ssq/nu) + (penalty * sum_abs<C, scalar_t>(pars, 0, p+q));
  }
  if constexpr(!arima_only) {
    return (ssq/nu) + (penalty * sum_abs<C, scalar_t>(pars));
  }
}

template <typename scalar_t> scalar_t lambda_max(
    const std::vector<scalar_t> &y,
    const size_t max_lag = 3 ) {
  const size_t n = y.size();
  scalar_t result = 0;
  for(size_t lag=1; lag < max_lag; lag++) {
    // get absolute scale magnitude for given lag
    scalar_t sum = 0;
    for( size_t j = lag; j < n; j++ ) {
      sum += y[j] * y[j-lag];
    }
    // take absolute value of sum
    sum = std::abs(sum);
    // if sum  is largest value so far, assign to result, otherwise keep result
    result = result < sum ? sum : result;
  }
  // divide by n to normalize
  return result/n;
}

// build lambda path for a given time series
template <typename scalar_t,
          const bool reverse=false> std::vector<scalar_t> lambda_path(
    const std::vector<scalar_t> &y,
    const size_t max_lag = 3,
    const size_t K = 50,
    const scalar_t epsilon = 0.005) {
  scalar_t max_lambda = lambda_max<scalar_t>(y, max_lag);
  // build sequence on log scale, then exponentiate
  std::vector<scalar_t> lambda_path = regular_sequence<scalar_t, reverse>(
    std::log(max_lambda), std::log(max_lambda*epsilon), K);
  // exponentiate back from log sequence
  for( auto & val:lambda_path) {
    val = std::exp(val);
  }
  return lambda_path;
}

// return true if this lambda
template <typename C, typename scalar_t> bool lambda_calibration_check(
    const C &low_coefs,
    const C &high_coefs,
    const scalar_t higher_lambda,
    const scalar_t lower_lambda,
    const scalar_t constant) {

  const size_t n = low_coefs.size();
  scalar_t largest_coef = 0;
  for( size_t i = 0; i < n; i++ ) {
    scalar_t current_coef = std::abs(low_coefs[i] - high_coefs[i]);
    // set to largest coef
    largest_coef = largest_coef > current_coef ? largest_coef : current_coef;
  }
  return ((largest_coef/(higher_lambda + lower_lambda)) - constant) <= 0;
}

template <const bool has_xreg, const bool seasonal, typename scalar_t = double,
          const bool arima_only=false>
class ARIMA_LASSO_CSS_PROBLEM :
  public cppoptlib::function::Function<
    scalar_t,
    ARIMA_LASSO_CSS_PROBLEM<has_xreg, seasonal, scalar_t, arima_only>> {

  using EigVec = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;
  using EigMat = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;
  using StateXd = cppoptlib::function::State<scalar_t, EigVec, EigMat>;

private:
  const arima_kind kind;
  const bool intercept;
  const size_t n_cond, n;
  const std::vector<scalar_t> y;
  scalar_t lambda;
  // initialized in constructor
  EigMat xreg;
  size_t arma_pars;
  EigVec y_temp, new_x;
  std::vector<scalar_t> residual, temp_phi, temp_theta;
public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  // initialize with a given arima structure
  ARIMA_LASSO_CSS_PROBLEM(const std::vector<scalar_t> &y,
                    const arima_kind &kind,
                    const bool intercept,
                    const bool drift,
                    std::vector<std::vector<scalar_t>> &xreg,
                    const scalar_t lambda = -1)
    : kind(kind), intercept(intercept),
      n_cond(kind.d() + (kind.D() * kind.period()) +
             kind.p() + (kind.P() * kind.period())),
      n(y.size()), y(y), lambda(lambda) {
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
    this->temp_phi = std::vector<scalar_t>(
      kind.p() + (kind.P() * kind.period()));
    this->temp_theta = std::vector<scalar_t>(
      kind.q() + (kind.Q() * kind.period()));
  }
  scalar_t operator()(const EigVec &x) {
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
      arima_transform_parameters<EigVec, seasonal, false, false, scalar_t>(
          this->new_x,
          this->kind,
          this->temp_phi,
          this->temp_theta
      );
    }
    scalar_t res = arima_css_ssq_lasso<EigVec, scalar_t, arima_only>(
      this->y_temp, this->new_x, this->kind,
      this->lambda, this->n_cond, this->residual);
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
  void reset_lambda(const scalar_t lambda) {
    this->lambda = lambda;
  }
  void finalize( structural_model<scalar_t> &model,
                 const EigVec & final_pars,
                 scalar_t kappa,
                 SSinit ss_init,
                 std::vector<scalar_t> &residuals) {
    // this function creates state space representation and expands it
    // I found out it is easier and cheaper (computationally) to do here
    // do the same for model coefficients
    // finally, make state space model
    structural_model<scalar_t> arima_ss = make_arima(
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
    for(size_t i=0; i < residuals.size(); i++) {
      residuals[i] = this->residual[i];
    }
  }
};

template <const bool has_xreg, const bool seasonal, typename scalar_t>
scalar_t lasso_solver(
    std::vector<scalar_t> &y,
    const arima_kind &kind,
    structural_model<scalar_t> &model,
    std::vector<std::vector<scalar_t>> xreg,
    const bool intercept,
    const bool drift,
    std::vector<scalar_t> &coef,
    const scalar_t kappa,
    const SSinit ss_init,
    std::vector<scalar_t> &residuals) {
  using EigVec = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;
  EigVec x(coef.size());
  // specify arima lasso path - in reverse (high to low lambda)
  auto lam_path = lambda_path<scalar_t, false>(
    y, std::max(kind.p(), kind.q()), 100);
  // std::cout << "lambda path: " << std::endl;
  // print_vector(lam_path);
  std::vector<EigVec> coefs_along_path(100, EigVec(coef.size()));
  for(size_t i = 0; i < coef.size(); i++) x(i) = coef[i];
  // initialize solver
  cppoptlib::solver::Bfgs<
    ARIMA_LASSO_CSS_PROBLEM<has_xreg, seasonal, scalar_t, true>> solver;
  // and arima problem
  ARIMA_LASSO_CSS_PROBLEM<has_xreg, seasonal, scalar_t, true> css_arima_problem(
      y, kind, intercept, drift, xreg);
  css_arima_problem.reset_lambda(lam_path[0]);
  // and finally, minimize
  auto [solution, solver_state] = solver.Minimize(css_arima_problem, x);
  // set up current value vector
  std::vector<scalar_t> values(100);
  // copy current solution parameters
  coefs_along_path[0] = solution.x;
  values[0] = solution.value;
  size_t j = 1;
  for(; j < 100; j++) {
    css_arima_problem.reset_lambda(lam_path[j]);
    // and finally, minimize
    auto [solution, solver_state] = solver.Minimize(css_arima_problem, x);
    // copy current solution parameters and solution value
    coefs_along_path[j] = solution.x;
    values[j] = solution.value;
    // std::cout << "Current lambda: " << lam_path[j] << std::endl;
    // copy current values
    if(!lambda_calibration_check(coefs_along_path[j-1], coefs_along_path[j],
                      // 0.75 is hardcoded, but this is value recommended by
                      // the original paper
                      lam_path[j-1], lam_path[j], 0.75 )) {
      // std::cout << "Chosen lambda: " << lam_path[j-1] << std::endl;
      break;
    }
  }
  // finalize
  css_arima_problem.finalize(
    model, coefs_along_path[j-1], kappa, ss_init, residuals);
  // pass fitted coefficients back to the caller
  for (size_t i = 0; i < coef.size(); i++) coef[i] = coefs_along_path[j-1][i];
  return exp(2 * values[j-1]);
}

#endif
