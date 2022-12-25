#ifndef BOXCOX
#define BOXCOX

#include <cmath>

template <typename T> T small_max(T a, T b) {
  if (a > b) {
    return a;
  }
  return b;
}

template <typename NumKind = double>
NumKind lambda_coef(std::vector<NumKind> &x, NumKind lambda, int period = 2) {
  // pre-compute the number of windows to compute over
  // this is integer division so the if the last period does not exactly overlap
  // it should only contain however many items are available
  int num_buckets = ceil((float)x.size() / period);

  std::vector<NumKind> mu_ests(num_buckets);
  std::vector<NumKind> sigma_ests(num_buckets);
  // split data into buckets defined by period and compute averages
  for (int i = 0; i < num_buckets; i++) {
    for (int j = 0; j < period; j++) {
      // i indexes the seasonal period - the bucket
      // j index step within seasonal period
      // sum within the the estimate
      mu_ests[i] = mu_ests[i] + x[(i * period) + j];
    }
    // now divide this by the number of seasonal periods to get the mean for mus
    mu_ests[i] = mu_ests[i] / period;
  }
  // now compute the standard errors
  for (int i = 0; i < num_buckets; i++) {
    for (int j = 0; j < period; j++) {
      // i indexes the seasonal period - the bucket
      // j index step within seasonal period
      // sum the squares within the estimate
      sigma_ests[i] = sigma_ests[i] + pow(x[(i * period) + j] - mu_ests[i], 2);
    }
    // now divide this by the number of seasonal periods - 1 ( the N-1 formula
    // for variance ) (optionally 1 in case we only have a period of 1)
    sigma_ests[i] = sigma_ests[i] / small_max(period - 1, 1);
    // finally we need to take the square root of this
    sigma_ests[i] = sqrt(sigma_ests[i]);
  }
  // now compute the ratios
  for (int i = 0; i < num_buckets; i++) {
    // we can very happy reuse the mu_ests without having to allocate more
    // memory
    mu_ests[i] = sigma_ests[i] / pow(mu_ests[i], 1.0 - lambda);
  }
  // compute standard deviation divided by the mean
  NumKind final_mu = 0.0, final_sigma = 0.0;
  for (int i = 0; i < num_buckets; i++) {
    final_mu = final_mu + mu_ests[i];
    final_sigma = final_sigma + pow(mu_ests[i], 2);
  }
  final_mu = final_mu / (NumKind)num_buckets;
  final_sigma = final_sigma / (NumKind)num_buckets;
  // // subtract mean
  final_sigma = final_sigma - pow(final_mu, 2);
  final_sigma = sqrt(final_sigma);
  return final_sigma / final_mu;
}

template <typename T> struct box_cox_obj_wrap {
  box_cox_obj_wrap(std::vector<T> &data, T par) : x(data), par(par) {}
  T operator()(T par_val) {
    // when you call this with par, it provides an evaluation of par
    // using the objective and the data. I am sure this can be turned to
    // a variable argument function with template metaprogramming...
    // but that misses the point
    return lambda_coef(x, par_val);
  }
  std::vector<T> &x;
  T par;
};

#endif
