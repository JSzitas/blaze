#ifndef BOXCOX
#define BOXCOX

#include <cmath>

template <typename scalar_t=float> class BoxCoxTransformer{
  using Vec = std::vector<scalar_t>;
  // estimated lambda
  scalar_t lambda;
public:
  BoxCoxTransformer<scalar_t>(const Vec &x,
                              const size_t period = 2,
                              const scalar_t lambda = -99999,
                              const scalar_t lower = -2,
                              const scalar_t upper = 3,
                              const scalar_t tol = 1e-6,
                              const scalar_t eps = 1e-8,
                              const size_t max_iter = 200) {
    if( lambda == -99999) {
      this->lambda = solve_brent(x, period, lower, upper, tol, eps, max_iter);
    }
    else{
      this->lambda = lambda;
    }
  }
  const Vec transform(const Vec &y) {
    Vec result(y.size());
    if( lambda < 0)
      std::transform(y.cbegin(), y.cend(),
                     std::back_inserter(result),
                     [this](scalar_t val){
                       if( val < 0 ) val = std::nan("");
                       return val;
                     });
    else if(lambda == 0)
      std::transform(y.cbegin(), y.cend(),
                     std::back_inserter(result),
                     [this](scalar_t val){
                       return std::log(val);
                     });
    else
      std::transform(y.cbegin(), y.cend(), std::back_inserter(result),
        [this](scalar_t val){
          return (sgn(val) * std::pow(std::abs(val), this->lambda - 1))/this->lambda;
          });
    return result;
  }
  const Vec inverse_transform(const Vec &y) {
    std::vector<scalar_t> result;
    if(lambda < 0)
      std::transform(y.cbegin(), y.cend(),
                     std::back_inserter(result),
                     [this](scalar_t val){
        if(val > (-1/this->lambda) ) val = std::nan("");
        return val;
      });
    else if(lambda == 0)
      std::transform(y.cbegin(), y.cend(),
                     std::back_inserter(result),
                     [this](scalar_t val){
                      return std::exp(val);
                     });
    else
      std::transform(y.cbegin(), y.cend(),
                     std::back_inserter(result),
                     [this](scalar_t val){
        val *= (this->lambda + 1);
        return sgn(val) * std::pow(std::abs(val),(1/this->lambda));
      });
    return result;
  }
  const scalar_t inverse_transform(const scalar_t y) {
    if(lambda < 0) {
      if(y > (-1/this->lambda) )
        return std::nan("");
    }
    else if(lambda == 0)
      return std::exp(y);
    else {
      y *= (this->lambda + 1);
      return sgn(y) * std::pow(std::abs(y),(1/this->lambda));
    }
  }
private:
  const scalar_t sgn( const scalar_t x ) {
    return (scalar_t(0) < x) - (x < scalar_t(0));
  }
  const scalar_t ab_sign(const scalar_t a, const scalar_t b) {
    return b >= 0.0 ? std::abs(a) : -std::abs(a);
  }
  // Brent optimizer specialized for solving box cox
  // Adapted from Numerical Recipes in C
  scalar_t solve_brent(const Vec &x, const size_t period,
                       const scalar_t lower, const scalar_t upper,
                       const scalar_t tol, const scalar_t eps,
                       const size_t max_iter) {
    // values
    scalar_t a = lower, b = upper, c = upper, d, e, min1, min2;
    // function evaluations
    scalar_t fa = lambda_loss(x, period, a), fb = lambda_loss(x, period, b),
      fc, p, q, r, s, tol1, xm;
    // this is the case where the initial values are dumb
    if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) return -999999999;
    fc = fb;
    for (size_t iter=0; iter < max_iter; iter++) {
      if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
        c = a; fc = fa; e = d = b - a;
      }
      if (std::abs(fc) < std::abs(fb)) {
        a = b; b = c; c = a; fa = fb; fb = fc; fc = fa;
      }
      //Convergence check.
      tol1 = 2.0 * eps * std::abs(b) + 0.5 * tol;
      xm = 0.5 * (c - b);
      if (std::abs(xm) <= tol1 || fb == 0.0) return b;
      if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb)) {
        //Attempt inverse quadratic interpolation.
        s = fb / fa;
        if (a == c) {
          p = 2.0 * xm * s; q = 1.0 - s;
        } else {
          q = fa / fc; r = fb / fc;
          p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
          q = (q - 1.0) * (r - 1.0) * (s - 1.0);
        }
        //Check whether in bounds.
        if (p > 0.0) q = -q;
        p = std::abs(p);
        min1 = 3.0 * xm * q - std::abs(tol1 * q);
        min2 = std::abs(e * q);
        if (2.0 * p < (min1 < min2 ? min1 : min2)) {
          //Accept interpolation.
          e = d; d = p / q;
        } else {
          //Interpolation failed, use bisection.
          d = xm; e = d;
        }
      } else { //Bounds decreasing too slowly, use bisection.
        d = xm; e = d;
      }
      a = b; //Move last best guess to a.
      fa = fb;
      //Evaluate new trial root.
      if (std::abs(d) > tol1) b += d;
      else b += ((xm >= 0.0) ? std::abs(tol1) : -std::abs(tol1));
      fb = lambda_loss(x, period, b);
    }
    return -999999999; //Never get here.
  }
  scalar_t lambda_loss(const Vec &x, const size_t period, const scalar_t lambda) {
    // pre-compute the number of windows to compute over
    // this is integer division so the if the last period does not exactly overlap
    // it should only contain however many items are available
    const size_t num_buckets = ceil((scalar_t)x.size() / period);
    Vec mu_ests(num_buckets), sigma_ests(num_buckets);
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
};

#endif
