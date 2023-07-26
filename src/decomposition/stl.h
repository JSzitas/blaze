#ifndef BLAZE_STL
#define BLAZE_STL

#include "utils/median.h"
#include "utils/seasonality_utils.h"

#include "decomposition/decompose_result.h"

template <typename scalar_t> class STL {
  std::vector<scalar_t> y;
  // these have to be ints (probably)
  int n, period, seasonal, trend_size,
  low_pass_size, seasonal_degree, trend_degree, low_pass_degree,
  seasonal_jump, trend_jump, low_pass_jump;
  // these must be size_t
  size_t inner_iter, outer_iter;
  bool robust;
  std::vector<scalar_t> season, trend, residual, robust_weights;
  // temporaries
  std::vector<scalar_t> detrended, cycle_sub, cycle, seasonal_s, rw_temp;
public:
  STL<scalar_t>(){}
  STL<scalar_t>(const std::vector<scalar_t> &y, 
                const int period = 0, 
                const int seasonal = 7, 
                const int trend_size = 0, 
                const int low_pass_size = 0, 
                const int seasonal_degree = 1, 
                const int trend_degree = 1, 
                const int low_pass_degree = 1,
                const int seasonal_jump = 1, 
                const int trend_jump = 1, 
                const int low_pass_jump = 1, 
                const bool robust = false) : 
  y(y), n(y.size()), inner_iter(5 - robust*3), outer_iter(1 + robust*14),
  period(period), seasonal(seasonal), trend_size(trend_size),
  low_pass_size(low_pass_size), seasonal_degree(seasonal_degree),
  trend_degree(trend_degree), low_pass_degree(low_pass_degree),
  seasonal_jump(seasonal_jump), trend_jump(trend_jump),
  low_pass_jump(low_pass_jump), robust(robust),
  season(std::vector<scalar_t>(y.size(),0.0)),
  trend(std::vector<scalar_t>(y.size(),0.0)),
  residual(std::vector<scalar_t>(y.size(),0.0)),
  robust_weights(std::vector<scalar_t>(y.size(), 1.0)) {
    if(period == 0) this->period = find_seasonalities(y, 1)[0];
    else this->period = period;
    if(seasonal < 7) {
      std::cout << "Length of seasonal period provided smaller than 7; setting to 7." << std::endl;
      this->seasonal = 7;
    } else {
      if(seasonal % 2 == 0) this->seasonal = 7;
    }
    if(trend_size < 3 || trend_size < this->period) {
      std::cout << "Trend smoothing window size must be greater than 3," <<
        " greater equal to seasonal period, and odd; setting to " <<
          "seasonal period (or next odd number) instead." << std::endl;
      this->trend_size = this->period + ((this->period % 2) == 0);
    }
    if(trend_size == 0) {
      this->trend_size = static_cast<int>( 
        std::ceil(1.5 * this->period/(1 - 1.5/this->seasonal))
      );
      // ensure that the trend smoothing window size is odd
      this->trend_size += (this->trend_size % 2) == 0;
    }
    if(low_pass_size < 3 || low_pass_degree < this->period) {
      this->low_pass_size = this->period + 1 + (this->period % 2);
    }
    this->detrended = std::vector<scalar_t>(this->n + (2 * this->period), 0.0);
    this->cycle_sub = std::vector<scalar_t>(this->n + (2 * this->period), 0.0);
    this->cycle = std::vector<scalar_t>(this->n + (2 * this->period), 0.0);
    this->seasonal_s = std::vector<scalar_t>(this->n + (2 * this->period), 0.0);
    this->rw_temp = std::vector<scalar_t>(this->n + (2 * this->period), 0.0);
  }
  void fit() {
    if(n == 1) {
      std::cout << "Only 1 valid observation found - setting trend equal to it" 
                << std::endl;
      this->residual[0] = 0; 
      this->season[0] = 0; 
      this->trend[0] = this->y[0];
      return;
    }
    // convert user supplied runtime argument to precompiled code 
    if(this->robust) inner_step<true>();
    else inner_step<false>();
    // in case of additional runs, it is basically certain we are running 
    // with robustness weights 
    for(size_t i = 1; i < this->outer_iter; i++) {
      inner_step<true>();
      robustness_weights_update();
    }
    for(size_t i = 0; i < this->n; i++) {
      this->residual[i] = this->y[i] - this->season[i] - this->trend[i];
    }
    return;
  }
  const std::vector<scalar_t> get_trend() const { return this->trend; }
  const std::vector<scalar_t> get_season() const { return this->season; }
  const std::vector<scalar_t> get_remainder() const { return this->residual; }
  const DecompositionResult<scalar_t> get_decomposition_result() const {
    return DecompositionResult<scalar_t>(this->trend, this->season, this->residual);
  }
  void print_summary() {
    std::cout << "Period: " << period << " | seasonal: " << seasonal << 
      " | trend_size: " << trend_size << " | low_pass_size: " <<
        low_pass_size << std::endl;
  }
private:
  template <const bool robust_update=false> void inner_step() {
    for(size_t j = 0; j < inner_iter; j++) {
      // step 1 - detrending
      for(int i = 0; i < this->n; i++) this->detrended[i] = y[i] - trend[i];
      // step 2 - cycle-subseries smoothing
      subseries_smoothing<robust_update>();
      // step 3 - low pass filtering of smoothed cycle subseries
      low_pass_filtering_moving_averages();
      // and the corresponding loess smoothing
      loess_smoothing<false>(this->cycle, this->seasonal_s, this->detrended, this->rw_temp, 
                             n, low_pass_size, low_pass_degree, low_pass_jump);
      // step 4/step 5 - detrend smoothed cycle subseries and deseazonalize
      for(int i = 0; i < this->n; i++) {
        season[i] = this->cycle_sub[period + i] - this->detrended[i];
        this->detrended[i] = y[i] - season[i];
      }
      // step 6 - trend smoothing
      loess_smoothing<robust_update>(this->detrended, robust_weights, trend, this->cycle,
                                     n, trend_size, trend_degree, trend_jump);
    }
  }
  void moving_average( 
      std::vector<scalar_t> & average,
      std::vector<scalar_t> &x, 
      const int n,
      const int seasonal_length) {
    // this is advantageous since you can sum the first seasonal_length values 
    // within the same register
    float v = 0.0;
    for(int i = 0; i < seasonal_length; i++) v += x[i];
    average[0] = v;
    // and here we avoid 2 dereferences on average[j] and average[j+1]
    for(int j = 0; j < n - seasonal_length; j++) {
      v += x[seasonal_length + j] - x[j];
      average[j+1] = v;
    }
    // this should be pretty easy to auto-vectorize; godbolt indicates that 
    // e.g. clang pretty much always manages for floats and doubles
    for(auto &val:average) val /= seasonal_length;
  }
  void low_pass_filtering_moving_averages() {
    // three distinct moving averages 
    moving_average(this->cycle, this->cycle_sub, this->n + 2 * this->period, period);
    moving_average(this->detrended, this->cycle, this->n + this->period + 1, period);
    moving_average(this->cycle, this->detrended, this->n + 2, 3);
  }
  void robustness_weights_update() {
    for(int i = 0; i < this->y.size(); i++) {
      this->robust_weights[i] = std::abs(this->y[i] - this->trend[i] - this->season[i]);
    }
    scalar_t cmad = 6 * median(robust_weights);
    // if effectively zero 
    if(cmad < 0.0001) {
      for(auto &weight:robust_weights) weight = 1.0;
      return;
    }
    scalar_t c9 = 0.999 * cmad, c1 = 0.001 * cmad;
    for(auto &weight:robust_weights) {
      weight = weight <= c1 ? 1.0 : 
               weight <= c9 ? std::pow(1.0 - std::pow(weight/cmad, 2), 2) :
               0.0;
    }
  }
  template <const bool robust_update = false, const size_t ys_offset = 0>
  void loess_smoothing(const std::vector<scalar_t> &y, 
                       const std::vector<scalar_t> &robustness_weights,
                       std::vector<scalar_t> & ys,
                       std::vector<scalar_t> &res,
                       const int n,
                       const int seasonal_length, 
                       const int polynomial_degree,
                       const int n_jump) {
    int n_left = 0, n_right = std::min(seasonal_length, n);
    const int newnj = std::min(n_jump, n - 1);
    if(seasonal_length >= n) {
      for(int i = 0; i < n; i += newnj) {
        ys[i+ ys_offset] = estimate_location<robust_update>(
          y, res, robustness_weights, n, seasonal_length, polynomial_degree,
          i+1, n_left, n_right);
      }
    }
    else if(newnj == 1) {
      const int nsh = (seasonal_length + 2)/2;
      for(int i = 0; i < n; i++) {
        if((i + 1) > nsh && n_right != n) {
          n_left++; n_right++;
        }
        ys[i + ys_offset] = estimate_location<robust_update>(
          y, res, robustness_weights, n, seasonal_length, polynomial_degree,
          i + 1, n_left, n_right);
      }
    } else {
      const int nsh = (seasonal_length + 1)/2;
      for(int i = 0; i < n; i+= newnj) {
        // these if else statements are just indexing updates 
        if((i + 1) < nsh) {}
        else if((i+1) < (n - nsh + 1)) {
          n_left = i + 1 - nsh;
          n_right = seasonal_length + i + 1 - nsh;
        }
        else {
          n_left = n - seasonal_length;
          n_right = n;
        }
        ys[i + ys_offset] = estimate_location<robust_update>(
          y, res, robustness_weights, n, seasonal_length, polynomial_degree,
          i+1, n_left, n_right);
      }
    }
    if(std::isnan(ys[n + ys_offset-1])) ys[n + ys_offset-1] = y[n-1];
    // if we did not jumps we are done 
    if(newnj == 1) return;
    // otherwise we have a fun interpolating bit 
    for(int i = 0; i < (n - newnj); i += newnj) {
      scalar_t delta = (ys[i + newnj] - ys[i]) / newnj;
      for(int j = i; j < i + newnj; j++) {
        ys[j] = ys[i] + delta * ((j + 1) - (i + 1));
      }
    }
    const int k = (static_cast<int>(n - 1) / static_cast<int>(newnj)) * newnj + 1;
    if(k != n) {
      ys[n - 1] = estimate_location<robust_update>(
        y, res, robustness_weights, n, seasonal_length, polynomial_degree, n, n_left, n_right);
      if(std::isnan(ys[n - 1])) ys[n - 1] = y[n - 1];
      if(k != (n - 1)) {
        scalar_t delta = (ys[n - 1] - ys[k - 1]) / (n - k);
        for(int j = k; j < n; j++) ys[j] = ys[k - 1] + delta * ((j + 1) - k);
      }
    }
  }
  // estimate the value for a given location 
  template <const bool robust_update = false>
  scalar_t estimate_location(const std::vector<scalar_t> &y,
                             std::vector<scalar_t> &weights,
                             const std::vector<scalar_t> & robustness_weights,
                             const int n,
                             const int seasonal_length, const int polynomial_degree,
                             const int xs, const int n_left,
                             const int n_right){
    // n_left is actually nleft-1 if comparing to original implementation
    scalar_t h = std::max(xs - (n_left + 1), n_right - xs);
    if(seasonal_length > n) h += static_cast<scalar_t>((seasonal_length - n)/2);
    const scalar_t h9 = .999 * h, h1 = .001 * h;
    scalar_t a = 0.0;
    for(int j = n_left; j < n_right; j++) {
      weights[j] = 0.0;
      const scalar_t r = std::abs(j + 1 - xs);
      if(r <= h9) {
        if(r <= h1) weights[j] = 1.0;
        else weights[j] = std::pow(1.0 - std::pow(r/h, 3), 3);
        if constexpr(robust_update) weights[j] *= robustness_weights[j];
        a += weights[j];
      }
    }
    if(a <= 0) return std::nan("0.0");
    for(int j = n_left; j < n_right; j++) weights[j] /= a;
    if(h > 0 && polynomial_degree > 0) {
      scalar_t a = 0.0, c = 0.0; 
      for(int j = n_left; j < n_right; j++) a += weights[j] * (j + 1);
      for(int j = n_left; j < n_right; j++) {
        c += weights[j] * std::pow((j + 1 - a), 2); 
      }
      if(sqrt(c) > (0.001 * (n - 1.0))) {
        scalar_t b = (xs - a)/c;
        for(int j = n_left; j < n_right; j++) {
          weights[j] *= ((b * (j + 1 - a)) + 1.0);
        } 
      }
    }
    scalar_t res = 0.0;
    for(int j = n_left; j < n_right; j++) res += weights[j] * y[j];
    return res;
  }
  template <const bool robust_update = false> void subseries_smoothing() {
    // for each seasonal period, carry out subseries based smoothing
    for(int j = 0; j < period; j++) {
      const int k = ((n - (j + 1))/period) + 1;
      for(int i = 0; i < k; i++)
        this->cycle[i] = this->detrended[i * period + j];
      if constexpr(robust_update) {
        for(int i = 0; i < k; i++) 
          this->rw_temp[i] = this->robust_weights[i * period + j];
      }
      // note that the ys_offset is set to 1
      loess_smoothing<robust_update, 1>(
          this->cycle, this->rw_temp, this->seasonal_s, season, k, seasonal,
          seasonal_degree, seasonal_jump);
      // ensuring the first observation is correct
      this->seasonal_s[0] = estimate_location<robust_update>(
        this->cycle, season, this->rw_temp, k, seasonal, seasonal_degree, 0,
        0, min(seasonal, k));
      if(std::isnan(this->seasonal_s[0]))
        this->seasonal_s[0] = this->seasonal_s[1];
      // same for the last one 
      this->seasonal_s[k + 1] = estimate_location<robust_update>(
        this->cycle, season, this->rw_temp, 
        k, seasonal, seasonal_degree, k+1, std::max(0, k - seasonal),  k);
      if(std::isnan(this->seasonal_s[k + 1]))
        this->seasonal_s[k + 1] = this->seasonal_s[k];
      // write back to the bigger temporary
      for(int m = 0; m < k+2; m++)
        this->cycle_sub[m * period + j] = this->seasonal_s[m];
    }
  }
};

#endif                                                                                                                                                            
