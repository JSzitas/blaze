#include <Rcpp.h>
using namespace Rcpp;

// include optimizer library
#include "third_party/eigen.h"
#include "third_party/optim.h"

#include "iostream"
#include "optional"

// [[Rcpp::plugins(cpp17)]]


using FunctionXd = cppoptlib::function::Function<double>;

class Rosenbrock : public FunctionXd {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  double operator()(const Eigen::VectorXd &x) const {
    const double t1 = (1 - x[0]);
    const double t2 = (x[1] - x[0] * x[0]);
    return   t1 * t1 + 100 * t2 * t2;
  }
};

// [[Rcpp::export]]
void test_solver(double x0, double x1) {

  using Solver = cppoptlib::solver::Bfgs<Rosenbrock>;

  Rosenbrock f;
  Eigen::VectorXd x(2);
  x << x0, x1;

  // Evaluate
  auto state = f.Eval(x);
  // std::cout << f(x) << " = " << state.value << std::endl;
  // std::cout << state.x << std::endl;
  // std::cout << state.gradient << std::endl;
  Solver solver;

  auto [solution, solver_state] = solver.Minimize(f, x);
  if (state.hessian) {
    std::cout << *(state.hessian) << std::endl;
  }
  // std::cout << "argmin " << solution.x.transpose() << std::endl;
  // std::cout << "f in argmin " << solution.value << std::endl;
  // std::cout << "iterations " << solver_state.num_iterations << std::endl;
  // std::cout << "solver status " << solver_state.status << std::endl;

  return;
}
