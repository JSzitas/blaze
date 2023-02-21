#ifndef BOXCOX
#define BOXCOX

#include <cmath>



template <typename scalar_t = double>
scalar_t lambda_coef(std::vector<scalar_t> &x, scalar_t lambda, const size_t period = 2) {
  // pre-compute the number of windows to compute over
  // this is integer division so the if the last period does not exactly overlap
  // it should only contain however many items are available
  const size_t num_buckets = ceil((float)x.size() / period);

  std::vector<scalar_t> mu_ests(num_buckets);
  std::vector<scalar_t> sigma_ests(num_buckets);
  // split data into buckets defined by period and compute averages
  for (size_t i = 0; i < num_buckets; i++) {
    for (size_t j = 0; j < period; j++) {
      // i indexes the seasonal period - the bucket
      // j index step within seasonal period
      // sum within the the estimate
      mu_ests[i] = mu_ests[i] + x[(i * period) + j];
    }
    // now divide this by the number of seasonal periods to get the mean for mus
    mu_ests[i] = mu_ests[i] / period;
  }
  // now compute the standard errors
  for (size_t i = 0; i < num_buckets; i++) {
    for (size_t j = 0; j < period; j++) {
      // i indexes the seasonal period - the bucket
      // j index step within seasonal period
      // sum the squares within the estimate
      sigma_ests[i] = sigma_ests[i] + pow(x[(i * period) + j] - mu_ests[i], 2);
    }
    // now divide this by the number of seasonal periods - 1 ( the N-1 formula
    // for variance ) (optionally 1 in case we only have a period of 1)
    sigma_ests[i] = sigma_ests[i] / std::max(period - 1, (size_t)1);
    // finally we need to take the square root of this
    sigma_ests[i] = sqrt(sigma_ests[i]);
  }
  // now compute the ratios
  for (size_t i = 0; i < num_buckets; i++) {
    // we can very happy reuse the mu_ests without having to allocate more
    // memory
    mu_ests[i] = sigma_ests[i] / pow(mu_ests[i], 1.0 - lambda);
  }
  // compute standard deviation divided by the mean
  scalar_t final_mu = 0.0, final_sigma = 0.0;
  for (size_t i = 0; i < num_buckets; i++) {
    final_mu = final_mu + mu_ests[i];
    final_sigma = final_sigma + pow(mu_ests[i], 2);
  }
  final_mu = final_mu / (scalar_t)num_buckets;
  final_sigma = final_sigma / (scalar_t)num_buckets;
  // subtract mean
  final_sigma = final_sigma - pow(final_mu, 2);
  final_sigma = sqrt(final_sigma);
  return final_sigma / final_mu;
}

template <typename scalar_t> struct box_cox_obj_wrap {
  box_cox_obj_wrap(std::vector<scalar_t> &data, scalar_t par) : x(data), par(par) {}
  scalar_t operator()(scalar_t par_val) {
    return lambda_coef(x, par_val);
  }
  std::vector<scalar_t> &x;
  scalar_t par;
};

#endif
