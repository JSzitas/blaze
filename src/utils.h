#ifndef UTILS
#define UTILS

template <typename T> T abs( T x) {
  if( x < 0.0 ) {
    x = x - 2*x;
  }
  return x;
}

template <typename T> bool is_same( T &a, T &b, T tol = 0.000001 ) {
  return abs(a - b) < tol;
}

template <class T, typename U> bool all_const( T &a, U tol ) {
  for( long long unsigned int i=1; i < a.size(); i++) {
    if( !is_same(a[i], a[0], tol) ) {
      return false;
    }
  }
  return true;
}

template <typename T, class Generator> T runif( T& lower, T& upper, Generator& gen ) {
  T sample = (T)gen.yield();
  return ((upper - lower) * sample) + lower;
}

template <class T, class F> T map( T& x, F& fun ) {
  T result;
  result.reserve(x.size());
  for( auto &item:x ) {
    result.push_back(fun(item));
  }
  return result;
}

template <typename T> T mean( std::vector<T> & x ) {
  T result = 0;
  for( auto &item:x ) {
    result += item;
  }
  return result/(T)x.size();
}

template <class T> void print_vector( T& a ) {
  for(int i=0; i < (a.size()-1); i++) {
    std::cout << a[i] << ", ";
  }
  std::cout << a[(a.size()-1)] << std::endl;
}

template <typename T> T max(T&a, T&b) {
  return a < b ? b : a;
}

template <typename T> T max(T a, T b) {
  return a < b ? b : a;
}

template <typename T> T min(T&a, T&b) {
  return a > b ? b : a;
}

template <typename T> T min(T a, T b) {
  return a > b ? b : a;
}

#endif
