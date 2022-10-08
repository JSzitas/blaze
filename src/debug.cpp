#include <Rcpp.h>
using namespace Rcpp;

#include "boxcox.h"
#include "solvers.h"
#include "utils.h"

// [[Rcpp::export]]
double some_fun( std::vector<double> &x,
                              double lambda,
                              int period = 2) {
  return lambda_coef(x, lambda);
}


// [[Rcpp::export]]
double particle_solver_bc( std::vector<double> x,
                           double lower = -3,
                           double upper = 3,
                           int n_particles = 20 ) {
  // we pass this rather than the raw function so we can get 'impure'
  // objective function - ie it holds its own data - and thus just use
  // the map function without having to worry about passing additional
  // arguments. additionally, maybe we just shift to passing this objective
  // wrapper as an argument to this function
  box_cox_obj_wrap<double> f_wrapper( x, 0 );
  return particle_drainer( f_wrapper, x, lower, upper, n_particles);
}

// void debug_all_same( std::vector<double> x) {
//   auto val = all_const(x, 0.0001);
//
//   // std::cout << "Abs: " << abs(x) << std::endl;
//
//   // auto val = is_same( x,y );
//
//   std::cout << "Are all same? " << val << std::endl;
// }


