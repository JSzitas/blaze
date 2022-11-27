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

template <class T> void print_vec_as_sqr_mat( T& a ) {

  int nrow = sqrt(a.size());
  int ncol = nrow;

  for(int i=0; i < nrow; i++) {
    for(int j=0; j < ncol; j++) {
      std::cout << a[i+(ncol * j)] << ", ";
    }
    std::cout << "|" << std::endl;
  }
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

template <typename T> void pop_front(std::vector<T>& x)
{
  x.erase(x.begin());
}

template <class T> T diff( T&a, int lag = 1, int d = 1 ) {
  int p = a.size();
  T result = a;
  for( int j=0; j < d; j++ ) {
    for( int i = lag; i < p-j; i++ ) {
      result[i-lag] = result[i] - result[i-lag];
    }
  }
  result.resize(p - (d*lag));
  return result;
}

template <class U> std::vector<std::vector<U>> diff(
    std::vector<std::vector<U>> &a,
    int lag = 1,
    int d = 1
  ) {
  std::vector<std::vector<U>> result( a.size() );
  for(int i=0; i < a.size(); i++) {
    result[i] = diff( a[i], lag, d );
  }
  return result;
}

// defining this as void and modifying the first element is probably
// pretty sensible, seeing how this operator works
template <typename U=double> void operator-=(std::vector<U> &a,
                                             std::vector<U> &b) {
  for( int i = 0; i < a.size(); i++ ) {
    a[i] -= b[i];
  }
}

#endif
