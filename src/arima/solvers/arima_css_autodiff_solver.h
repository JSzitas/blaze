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

// include autodiff library
#include "third_party/autodiff.h"
using namespace autodiff;


using FunctionXd = cppoptlib::function::Function<double>;
using StateXd = cppoptlib::function::State<double, Eigen::VectorXd, Eigen::MatrixXd>;

template <const bool has_xreg, const bool seasonal>
class ARIMA_CSS_PROBLEM : public FunctionXd {
private:
  arima_kind kind;
  size_t n_cond;
  lm_coef<double> xreg_pars;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> xreg;
  size_t arma_pars;
  size_t n;
  autodiff::VectorXvar new_x, transform_temp_phi, transform_temp_theta, xreg_temp;
  Eigen::VectorXd y, y_temp, residual;
public:
  // debug only:
  size_t f_evals;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    // initialize with a given arima structure
    ARIMA_CSS_PROBLEM(std::vector<double> &y, const arima_kind &kind,
                      lm_coef<double> &xreg_pars,
                      std::vector<std::vector<double>> &xreg, size_t n_cond)
      : kind(kind), n_cond(n_cond), xreg_pars(xreg_pars) {
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
      this->y_temp = Eigen::VectorXd(n);
      for (size_t i = 0; i < n; i++) {
        this->y_temp(i) = y[i];
      }
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
      this->y = Eigen::VectorXd(n);
      for( size_t i = 0; i < n; i++ ) {
        this->y(i) = this->y_temp(i);
      }
      this->n = n;
      this->arma_pars = kind.p() + kind.P() + kind.q() + kind.Q();
      // pre-allocate new_x
      this->new_x =
        autodiff::VectorXvar(kind.p() + (kind.P() * kind.period()) + kind.q() +
        (kind.Q() * kind.period()) + this->xreg.cols());
      this->xreg_temp = autodiff::VectorXvar(this->xreg.cols());
        //Eigen::VectorXd(this->xreg.cols());
      // pre-allocate model residuals
      this->residual = Eigen::VectorXd(n);
      // pre-allocate transformation helper vector - this is only necessary
      // for expanding seasonal models
      this->transform_temp_phi =
        autodiff::VectorXvar(kind.p() + (kind.P() * kind.period()));
      this->transform_temp_theta =
        autodiff::VectorXvar(kind.q() + (kind.Q() * kind.period()));
      // debug only:
      this->f_evals = 0;
    }

    double operator()(const Eigen::VectorXd &x) {
      return 3.0;
    }
    autodiff::var compute_loss( autodiff::VectorXvar &x ) {

        for (size_t i = 0; i < x.size(); i++) {
          this->new_x(i) = x(i);
        }
        // thse variablees track all operations
        autodiff::var ssq_xreg = 0.0, tmp_xreg = 0.0;

        if constexpr (has_xreg) {
          // refresh y_temp and load it with original y data
          for (size_t i = 0; i < this->n; i++) {
            this->y_temp(i) = this->y(i);
          }
          // load only xreg parameters here
          for( size_t i = 0; i < this->xreg_temp.size(); i++ ) {
            this->xreg_temp(i) = x(this->arma_pars + i);
          }
          // Eigen::VectorXd
          autodiff::MatrixXvar temp_xreg_mat = this->xreg.eval();
          autodiff::VectorXvar xreg_res = temp_xreg_mat * this->xreg_temp;
          for(size_t i = 0; i < this->n; i++) {
             tmp_xreg = y_temp(i) - xreg_res(i);
             y_temp(i) = y_temp(i) - autodiff::val(xreg_res(i));
             ssq_xreg = ssq_xreg + (tmp_xreg * tmp_xreg);
          }
          // this->y_temp = this->y_temp - xreg_res;
          // do differencing here
          if (!this->xreg_pars.has_intercept()) {
            for (size_t i = 0; i < this->kind.d(); i++) {
              for (size_t l = this->n - 1; l > 0; l--) {
                this->y_temp(l) -= this->y_temp(l - 1);

                // this->y_temp(l) = this->y_temp(l) - this->y_temp(l - 1);
              }
            }
            // seasonal differencing
            size_t ns = this->kind.period();
            for (size_t i = 0; i < this->kind.D(); i++) {
              for (size_t l = this->n - 1; l >= ns; l--) {
                this->y_temp(l) -= this->y_temp(l - ns);
                // this->y_temp(l) = this->y_temp(l) - this->y_temp(l - ns);
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
          arima_transform_parameters<autodiff::VectorXvar,
                                     true, false>(this->new_x, this->kind,
                                                      this->transform_temp_phi,
                                                      this->transform_temp_theta);
        }

        int ma_offset, nu = 0;
        const size_t p = kind.p() + kind.period() * kind.P(),
          q = kind.q() + kind.period() * kind.Q();

        // reseet tmp to zero since we used it to compute intermediate products previously
        autodiff::var ssq = 0.0, tmp = 0.0;
        for (size_t l = this->n_cond; l < this->n; l++) {
          ma_offset = min(l - this->n_cond, q);
          tmp = y_temp(l);
          for (size_t j = 0; j < p; j++) {
            tmp = tmp - x(j) * y_temp(l - j - 1);
          }
          // to offset that this is all in one vector, we need to
          // start at p and go to p + q
          for (size_t j = 0; j < ma_offset; j++) {
            tmp = tmp - x(p + j) * this->residual(l - j - 1);
          }
          this->residual(l) = autodiff::val(tmp);
          if (!isnan(autodiff::val(tmp))) {
            nu++;
            ssq = ssq + (tmp * tmp);
          }
        }
        autodiff::var result = 0.5 * log((ssq + ssq_xreg)/nu);
        return result;
    }
    StateXd Eval(const Eigen::VectorXd &x,
                 const int order = 2) {

      StateXd state(x.rows(), order);

      autodiff::VectorXvar new_vec(x.size());
      for( size_t i = 0; i < x.size(); i++) {
        new_vec(i) = x(i);
      }

      autodiff::var res = this->compute_loss(new_vec);
      Eigen::VectorXd grad;
      Eigen::MatrixXd hess = autodiff::hessian( res, new_vec, grad );

      std::cout << "g: " << grad << std::endl;
      std::cout << "H: " << hess << std::endl;

      state.x = x;
      state.value = (double)res;
      state.gradient = grad;
      state.hessian = hess;
      this->f_evals++;
      return state;
    }
    void finalize( structural_model<double> &model,
                   const Eigen::VectorXd & final_pars,
                   std::vector<double> & delta,
                   double kappa,
                   SSinit ss_init) {
      // this function creates state space representation and expands it
      // I found out it is easier and cheaper (computationally) to do here
      // do the same for model coefficients
      // finally, make state space model
      structural_model<double> arima_ss = make_arima( this->new_x,
                                                      delta, this->kind,
                                                      kappa, ss_init);
      model.set(arima_ss);
      // modify y_temp to acount for xreg
      for (size_t i = 0; i < this->n; i++) {
        this->y_temp(i) = this->y(i);
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
void arima_solver_css(std::vector<double> &y, structural_model<double> &model,
                      lm_coef<double> xreg_coef,
                      std::vector<std::vector<double>> xreg,
                      const arima_kind &kind, std::vector<double> &coef,
                      std::vector<double> &delta, const int n_cond,
                      const int n_available, const double kappa,
                      const SSinit ss_init, double &sigma2) {

  size_t vec_size = kind.p() + kind.q() + kind.P() + kind.Q() + xreg_coef.size();
  size_t arma_size = kind.p() + kind.q() + kind.P() + kind.Q();
  Eigen::VectorXd x(vec_size);
  for (auto &val : x) {
    val = 0.3;
  }
  // initialize to all zeroes except for xreg
  for (size_t i = arma_size; i < vec_size; i++) {
    x[i] = xreg_coef.coef[i - arma_size];
  }
  // initialize solver
  using Solver = cppoptlib::solver::Bfgs<ARIMA_CSS_PROBLEM<has_xreg, seasonal>>;
  Solver solver;
  // and arima problem
  ARIMA_CSS_PROBLEM<has_xreg, seasonal> css_arima_problem(y, kind, xreg_coef,
                                                          xreg, n_cond);
  // and finally, minimize
  auto [solution, solver_state] = solver.Minimize(css_arima_problem, x);
  // update variance estimate for the arima model - this was passed by reference
  sigma2 = exp(2 * solution.value);
  std::cout << " Function evaluations taken: " << css_arima_problem.f_evals << std::endl;
  std::cout << " Solver iterations: " << solver_state.num_iterations << std::endl;
  for (size_t i = 0; i < vec_size; i++) {
    std::cout << solution.x[i] << ",";
  }


  // solver_state.Status == IterationLimit ||      GradientNormViolation,  // Minimum norm in gradient vector has been reached.
  //     HessianConditionViolation  // Maximum condition number of hessian_t has been reached.
  // };

  // enum class Status {
  //   NotStarted = -1,
  //     Continue = 0,     // Optimization should continue.
  //     IterationLimit,   // Maximum of allowed iterations has been reached.
  //     XDeltaViolation,  // Minimum change in parameter vector has been reached.
  //     FDeltaViolation,  // Minimum chnage in cost function has been reached.
  //     GradientNormViolation,  // Minimum norm in gradient vector has been reached.
  //     HessianConditionViolation  // Maximum condition number of hessian_t has been
  //   // reached.
  // };
  // if we are to return standard errors as well
  // if constexpr(return_hessian) {
  //   // first get the computed numerical hessian
  //   auto state = css_arima_problem.Eval(solution.x);
  //   //std::cout << *(state.hessian) << std::endl;
  //   // available under state.hessian
  //   // next multiply by -1 and the number of available observations
  //   auto est_hessian = state.hessian * -1 * (double)n_available;
  //   // next, invert this
  //   // the result is the variance covariance matrix of coefficient estimates, except
  //   // for the fact that xreg coefficients (including intercept) have significantly
  //   // understated effects
  //   // since for typical use (e.g. coefficient standard errors) you only want the
  //   // diagonal elements anyways, you can simply multiply the diagonal entries
  //   // on those elements by the sample variance used in originally scaling the data
  // }
  // css_arima_problem.finalize( model, solution.x, delta, kappa, ss_init);
  // pass fitted coefficients back to the caller
  for (size_t i = 0; i < vec_size; i++) {
    coef[i] = solution.x[i];
  }
}

#endif
