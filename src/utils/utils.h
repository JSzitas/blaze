#ifndef UTILS
#define UTILS

// included mainly for isnan()
#include <math.h>

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

template <typename T> T mean(std::vector<T> &x) {
  T result = 0;
  for (auto &item : x) {
    result += item;
  }
  return result / (T)x.size();
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

// template <class T> void print_vec_as_sqr_mat(T &a) {
//
//   int nrow = sqrt(a.size());
//   int ncol = nrow;
//
//   for (int i = 0; i < nrow; i++) {
//     for (int j = 0; j < ncol; j++) {
//       std::cout << a[i + (ncol * j)] << ", ";
//     }
//     std::cout << "|" << std::endl;
//   }
// }

template <typename T> T max(T &a, T &b) { return a < b ? b : a; }
template <typename T> T max(T a, T b) { return a < b ? b : a; }
template <typename T> T min(T &a, T &b) { return a > b ? b : a; }
template <typename T> T min(T a, T b) { return a > b ? b : a; }

template <typename T, typename U> T max(T a, U b) { return (T)(a < b ? b : a); }
template <typename T, typename U> T max(T &a, U &b) { return (T)(a < b ? b : a); }

template <typename T, typename U> T min(T a, U b) { return (T)(a > b ? b : a); }
template <typename T, typename U> T min(T &a, U &b) { return (T)(a > b ? b : a); }


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
  for (size_t i = 0; i < size; i++) {
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
    U mean_val = 0;
    size_t i = 0;
    for(; i < x.size(); i++) {
      mean_val += x[i];
    }
    mean_val /= (U)i;
    i = 0;
    U sd_val = 0;
    for(;i <x.size();i++) {
      sd_val += pow(x[i] - mean_val,2);
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

#endif
