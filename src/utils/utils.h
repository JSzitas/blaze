#ifndef GENERAL_UTILS
#define GENERAL_UTILS

// included mainly for isnan()
#include <math.h>
#include <unordered_set>
#include <limits>
#include <random>

template <typename T> T abs(T x) {
  if (x < 0.0) {
    x = x - 2 * x;
  }
  return x;
}

template <typename T> bool is_same(T &a, T &b, T tol = 0.000001) {
  return abs(a - b) < tol;
}

template <class T, typename U> bool all_const(T &a, U tol) {
  for (size_t i = 1; i < a.size(); i++) {
    if (!is_same(a[i], a[0], tol)) {
      return false;
    }
  }
  return true;
}

template <typename T, class Generator>
T runif(T &lower, T &upper, Generator &gen) {
  T sample = (T)gen.yield();
  return ((upper - lower) * sample) + lower;
}

template <class T, class F> T map(T &x, F &fun) {
  T result;
  result.reserve(x.size());
  for (auto &item : x) {
    result.push_back(fun(item));
  }
  return result;
}

template <class T> void print_vector(T &a) {
  if( a.size() < 1 ) {
    return;
  }
  for (size_t i = 0; i < a.size(); i++) {
    std::cout << a[i] << ", ";
  }
  std::cout << " " << std::endl;
}

template <typename T> T max(T &a, T &b) { return a < b ? b : a; }
template <typename T> T max(T a, T b) { return a < b ? b : a; }
template <typename T> T min(T &a, T &b) { return a > b ? b : a; }
template <typename T> T min(T a, T b) { return a > b ? b : a; }

template <typename T, typename U> T max(T a, U b) { return (T)(a < b ? b : a); }
template <typename T, typename U> T max(T &a, U &b) { return (T)(a < b ? b : a); }

template <typename T, typename U> T min(T a, U b) { return (T)(a > b ? b : a); }
template <typename T, typename U> T min(T &a, U &b) { return (T)(a > b ? b : a); }

template <typename T> T max(const std::vector<T> &x) {
  T max = -1 * std::numeric_limits<T>::infinity();
  for(size_t j=0; j< x.size(); j++) {
    if(x[j] > max) {
      max = x[j];
    }
  }
  return max;
}
template <typename T> T min(const std::vector<T> &x) {
  T min = std::numeric_limits<T>::infinity();
  for(size_t j=0; j< x.size(); j++) {
    if(x[j] < min) {
      min = x[j];
    }
  }
  return min;
}
template <typename T> size_t max_at(const std::vector<T> &x) {
  size_t max_at = 0;
  for(size_t j=0; j< x.size(); j++) {
    if(x[j] > x[max_at]) {
      max_at = j;
    }
  }
  return max_at;
}
template <typename T> size_t min_at(const std::vector<T> &x) {
  size_t min_at = 0;
  for(size_t j=0; j< x.size(); j++) {
    if(x[j] < x[min_at]) {
      min_at = j;
    }
  }
  return min_at;
}
template <typename T> size_t max_at(
    const std::vector<T> &x, const size_t start, const size_t end) {
  size_t max_at = start;
  for(size_t j=start+1; j < end; j++) {
    if(x[j] > x[max_at]) {
      max_at = j;
    }
  }
  return max_at;
}
template <typename T> size_t min_at(
    const std::vector<T> &x, const size_t start, const size_t end) {
  size_t min_at = start;
  for(size_t j=start+1; j < end; j++) {
    if(x[j] < x[min_at]) {
      min_at = j;
    }
  }
  return min_at;
}

template <typename T> void pop_front(std::vector<T> &x) { x.erase(x.begin()); }

template <class T> T diff(T &a, size_t lag = 1, size_t d = 1) {
  size_t p = a.size();
  T result = a;
  for (size_t j = 0; j < d; j++) {
    for (size_t i = lag; i < p - j; i++) {
      result[i - lag] = result[i] - result[i - lag];
    }
  }
  result.resize(p - (d * lag));
  return result;
}

template <class U>
std::vector<std::vector<U>> diff(std::vector<std::vector<U>> &a,
                                 size_t lag = 1,
                                 size_t d = 1) {
  std::vector<std::vector<U>> result(a.size());
  for (size_t i = 0; i < a.size(); i++) {
    result[i] = diff(a[i], lag, d);
  }
  return result;
}

// defining this as void and modifying the first element is probably
// pretty sensible, seeing how this operator works
template <typename U = double>
void operator-=(std::vector<U> &a, std::vector<U> &b) {
  for (size_t i = 0; i < a.size(); i++) {
    a[i] -= b[i];
  }
}
// this is for usage with mixes of st::vector and Eigen::vector
template <class T, class U> void operator-=(T &a, U &b) {
  for (size_t i = 0; i < a.size(); i++) {
    a[i] -= b[i];
  }
}

template <typename U=double> size_t count_na(std::vector<U> &a) {
  size_t count = 0;
  for (auto &item : a) {
    if (std::isnan(item)) {
      count++;
    }
  }
  return count;
}

template <typename scalar_t> scalar_t mean(
    const std::vector<scalar_t> & x) {
  scalar_t result = 0;
  size_t n = x.size();
  for( auto &val:x ) {
    if(!std::isnan(val)) {
      result += val;
    }
    else {
      n--;
    }
  }
  return result/(scalar_t)n;
}


template <typename U = double> std::vector<size_t> find_na(std::vector<U> &a) {
  std::vector<size_t> result;
  result.reserve(a.size());
  for (size_t i = 0; i < a.size(); i++) {
    if (std::isnan(a[i])) {
      result.push_back(i);
    }
  }
  return result;
}

template <typename T>
std::vector<T> intersect(std::vector<T> &a, std::vector<T> &b) {
  std::unordered_set<T> comp_set(a.begin(), a.end());
  std::vector<T> result;
  result.reserve(min(a.size(), b.size()));
  for (auto &val : b) {
    if (comp_set.count(val)) {
      result.push_back(val);
    }
  }
  return result;
}

template <typename T>
std::vector<T> intersect(std::vector<T> &a, std::vector<T> b) {
  std::unordered_set<T> comp_set(a.begin(), a.end());
  std::vector<T> result;
  result.reserve(min(a.size(), b.size()));
  for (auto &val : b) {
    if (comp_set.count(val)) {
      result.push_back(val);
    }
  }
  return result;
}

template <typename U = double>
std::vector<U> flatten_vec(std::vector<std::vector<U>> x) {
 size_t size = 0;
  for (size_t i = 0; i < x.size(); i++) {
    size += x[i].size();
  }
  size_t p = 0;
  std::vector<U> result(size);
  for (size_t i = 0; i < x.size(); i++) {
    for (size_t j = 0; j < x[i].size(); j++) {
      result[p] = x[i][j];
      p++;
    }
  }
  return result;
}

template <typename U=double> struct StandardScaler{
  StandardScaler<U>(){
    this->mean = 0;
    this->sd = 0;
  }
  StandardScaler<U>(std::vector<U> &x) {
    U mean_val = 0.0;
    size_t i = 0;
    for(; i < x.size(); i++) {
      mean_val += x[i];
    }
    mean_val /= (U)i;
    i = 0;
    U sd_val = 0.0;
    for(;i <x.size(); i++) {
      sd_val += std::pow(x[i] - mean_val,2);
    }
    sd_val /= (U)(i-1);
    sd_val = sqrt(sd_val);
    this->mean = mean_val;
    this->sd = sd_val;
  }
  void scale(std::vector<U> &x) const {
    for( size_t i=0; i < x.size(); i++ ) {
      x[i] = (x[i] - this->mean)/this->sd;
    }
  }
  void rescale(std::vector<U> &x) const {
    for(size_t i=0; i < x.size(); i++) {
      x[i] = this->mean + (x[i]*this->sd);
    }
  }
  void rescale_w_mean(std::vector<U> &x) const {
    for(size_t i=0; i < x.size(); i++) {
      x[i] += this->mean;
    }
  }
  void rescale_w_sd(std::vector<U> &x) const {
    for(size_t i=0; i < x.size(); i++) {
      x[i] *= this->sd;
    }
  }
  U scale_val( const U x) const {
    return (x - this->mean)/this->sd;;
  }
  U rescale_val( const U x ) const {
    return this->mean + (x*this->sd);
  }
  U rescale_val_w_mean( const U x ) const {
    return this->mean + x;
  }
  const U get_mean() const {
    return this->mean;
  }
  const U get_sd() const {
    return this->sd;
  }
  const std::pair<U,U> summary() const {
    return std::pair(this->mean, this->sd);
  }
private:
  U mean, sd;
};

template <typename U=double> struct MinMaxScaler{
  MinMaxScaler<U>(std::vector<U> &x, U a=0, U b=1) {
    U min_val = x[0];
    U max_val = x[0];
    for( size_t i =1; i < x.size(); i++) {
      min_val = min_val < x[i] ? min_val : x[i];
      max_val = max_val > x[i] ? max_val : x[i];
    }
    this->min_v = min_val;
    this->max_v = max_val;
    this->a = a;
    this->b = b;
    this->spread = this->max_v - this->min_v;
    this->range = this->b - this->a;
  }
  void scale(std::vector<U> &x) const {
    for( size_t i=0; i < x.size(); i++ ) {
      x[i] = (x[i] - this->min_v)/this->spread;
      x[i] = (x[i] * this->range) + this->a;
    }
  }
  void rescale(std::vector<U> &x) const {
    for(size_t i=0; i < x.size(); i++) {
      x[i] = this->a + (x[i]/this->range);
      x[i] = (x[i] * this->spread) + this->min_v;
    }
  }
  U rescale_val( const U x ) const {
    auto y = this->a + (x/this->range);
    return (y * this->spread) + this->min_v;
  }
private:
  U spread, range, a,b, min_v, max_v;
};

template <typename scalar_t=double> size_t continguous_len(
  const std::vector<scalar_t> & y) {
  // find continguous non-missing size of y
  size_t first=0, last=0;
  for( size_t i = 0; i < y.size(); i++ ) {
    if(!isnan(y[i])) {
      last++;
    } else {
      first = i;
      last = i;
    }
  }
  return last - first;
}

template<typename scalar_t> bool all_is_nan(const std::vector<scalar_t> &x) {
  for( auto &val:x ) {
    if(!std::isnan(val))
      return false;
  }
  return true;
}

template<typename scalar_t> size_t size_omit_nan(const std::vector<scalar_t> &x) {
  size_t n = 0;
  for( auto &val:x ) {
    if(!std::isnan(val))
      n++;
  }
  return n;
}

template <typename scalar_t=double> bool is_constant(
  std::vector<scalar_t> &x,
  const scalar_t tol = 1e-7) {
  scalar_t current_val;
  for( auto & val:x ) {
    // find non-nan value to start with
    if(!isnan(val)) {
      current_val = val;
      // break out of loop
      break;
    }
  }
  for( auto & val:x ) {
    // if we find out that any two values are different, the vector is not
    // constant
    if( std::abs(val - current_val) > tol ) {
      return false;
    }
  }
  return true;
}

template <typename scalar_t> scalar_t sign(const scalar_t x) {
  return x > 0.0 ? 1.0 : -1.0;
}

template <typename scalar_t> scalar_t trunc(const scalar_t x) {
  return sign(x) * std::floor(std::abs(x));
}

template <typename scalar_t,
          const bool truncate = true>
std::vector<std::vector<scalar_t>> get_lags(
  const std::vector<scalar_t> &y,
  const size_t lags = 2) {

  const size_t n = y.size();
  std::vector<std::vector<scalar_t>> result(
      lags+1, std::vector<scalar_t>(n-lags));
  // for each lag 
  for( size_t i = 0; i < lags+1; i++ ){
    // every result starts at n_lag (since we skip rows where any column is missing)
    // thus what changes is which value of y is at that position
    for(size_t j = lags; j < n; j++) {
      result[i][j-lags] = y[j-i];
    }
  }
  return result;
}

template <typename scalar_t, const bool reverse=false>
std::vector<scalar_t> regular_sequence(
    scalar_t seq_min,
    scalar_t seq_max,
    const size_t size = 100) {
  scalar_t increment = (seq_min + seq_max)/size;
  if constexpr(reverse) {
    // swap seq_min and seq_max
    std::swap(seq_min, seq_max);
    increment *= -1.0;
  }
  std::vector<scalar_t> result(size+1, seq_min);
  for( size_t j = 0; j < size+1; j++ ) {
    result[j] += j*increment;
  }
  return result;
}

template <typename T> std::vector<T> cummulative_product(
    const std::vector<T> &y) {
  std::vector<T> result(y.size());
  for(size_t i=1; i < y.size(); i++) {
    result[i] *= result[i-1];
  }
  return result;
}

template <typename T> std::vector<T> cummulative_sum(
    const std::vector<T> &y) {
  std::vector<T> result(y.size());
  for(size_t i=1; i < y.size(); i++) {
    result[i] += result[i-1];
  }
  return result;
}

template <typename ForwardIt> void cummulative_sum(
    ForwardIt begin, ForwardIt end) {
  for(;begin != (end-1); begin++) {
    *(begin+1) = *(begin+1) + *(begin);
  }
}

template <typename ForwardIt> void reverse_scaling(
    ForwardIt begin, ForwardIt end) {
  auto c_size = end-begin;
  auto size = 1;
  for(;begin != end; begin++) {
    *(begin) = *(begin)/(c_size/size);
    size++;
  }
}

template <typename ForwardIt, typename T> void add_to_vec(
    ForwardIt begin, ForwardIt end, T val) {
  for(;begin != end; begin++) {
    *(begin) = *(begin) + val;
  }
}

template <typename scalar_t> scalar_t crossprod(
    const std::vector<scalar_t> &x) {
  scalar_t result = 0;
  for(size_t j = 0; j < x.size(); j++) {
    result += std::pow(x[j], 2);
  }
  return result;
}

template <typename T> T draw_from_vec(
    const std::vector<T> &x, std::mt19937& twister) {
  
  auto res = twister();
  constexpr std::uint64_t range = 2147483647;
  
  return x[x.size() * res/range];
}

// mostly just to have a clear implementation
template <typename scalar_t> std::array<scalar_t, 2> 
welfords_algorithm(const std::vector<scalar_t> &y) {
  std::array<scalar_t, 2> result;
  const size_t n = y.size();
  // mean, standard error
  result[0] = y[0];
  result[1] = 0;
  scalar_t prev_value = y[0];
  for(size_t j = 1; j < y.size(); j++) {
    prev_value = result[0];
    result[0] += (y[j] - prev_value)/static_cast<scalar_t>(j);
    result[1] += ((y[j] - prev_value) * (y[j] - result[0]));
  }
  result[1] = std::sqrt(result[1]/(n-1));
  return result;
}

template <typename scalar_t> void welfords_algorithm(
    std::array<scalar_t, 2> &current_state,
    const scalar_t current_val, const size_t index) {
  if(std::isnan(current_val)) return;
  scalar_t prev_value = current_state[0];
  current_state[0] += (current_val - prev_value)/static_cast<scalar_t>(index);
  current_state[1] += ((current_val - prev_value) * (current_val - current_state[0]));
}

template <typename scalar_t> void welfords_algorithm(
    scalar_t &current_mean, scalar_t &current_std_err,
    const scalar_t current_val, const size_t index) {
  if(std::isnan(current_val)) return;
  scalar_t prev_value = current_mean;
  current_mean += (current_val - prev_value)/static_cast<scalar_t>(index);
  current_std_err += ((current_val - prev_value) * (current_val - current_mean));
}


#endif
