#include <Rcpp.h>
using namespace Rcpp;

#include "utils/seasonality_utils.h"
#include "utils/utils.h"

#include <random>

// [[Rcpp::export]]
std::vector<size_t> find_seasons(std::vector<double>& x) {
  return find_seasonalities(x);
}

// [[Rcpp::export]]
size_t find_period(std::vector<double>& x) {
  return period(x);
}

// [[Rcpp::export]]
double draw(std::vector<double>& x) {
  
  auto twister = std::mersenne_twister_engine<size_t, 64, 312, 156, 31,
                                              0xb5026f5aa96619e9, 29,
                                              0x5555555555555555, 17,
                                              0x71d67fffeda60000, 37,
                                              0xfff7eee000000000, 43, 6364136223846793005>{std::random_device{}()};
  
  // auto twister = std::mt19937{std::random_device{}()};
  return draw_from_vec(x, twister);
}

