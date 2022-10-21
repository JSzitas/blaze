#include <Rcpp.h>
using namespace Rcpp;

#include "boxcox.h"
#include "solvers.h"
#include "utils.h"

#include "poly.h"

// #include "lm.h"

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

// void print_mat( small_lm::Matrix<double> x ) {
//   int p=0;
//   for( int i=0; i < x.nrow(); i++ ) {
//     for( int j=0; j< x.ncol(); j++ ) {
//       std::cout << x.elem(p) <<", ";
//       p++;
//     }
//     std::cout << std::endl;
//   }
// }
//
// void print_vector( std::vector<double> x ) {
//   for( auto &item:x ) {
//     std::cout << item << ", ";
//   }
//   std::cout << std::endl;
// }
//
// // [[Rcpp::export]]
// void test_inn_prod( std::vector<double> x,
//                     std::vector<double> y ) {
//   auto res = small_lm::inner_product(x,y);
//   std::cout << res <<std::endl;
// }
//
//
// // [[Rcpp::export]]
// void test_mat( std::vector<double> x,
//                std::vector<double> y,
//                int nrow = 5,
//                int ncol = 3,
//                int ncol2 = 7) {
//   small_lm::Matrix<double> mymat(x, nrow, ncol);
//   small_lm::Matrix<double> mymat2(y, ncol, ncol2);
//
//   // print_vector( mymat.col(0) );
//   // print_vector( mymat2.row(0));
//
//
//   auto res = mymat * mymat2;
//   print_mat(res);
// }

// [[Rcpp::export]]
std::vector<double> test_polyroot_abs( std::vector<double> &x ) {
  auto res = polyroot( x );
  std::vector<double> result(res.size());

  for( int i = 0; i < res.size(); i++ ) {
    result[i] = abs(res[i]);
  }

  // auto result = abs(res);
  return result;
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


