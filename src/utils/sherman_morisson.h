#ifndef SHERMAN_MORISSON
#define SHERMAN_MORISSON

// A O(n^2) implementation of the Sherman Morisson update formula 
template <typename scalar_t, const size_t max_size = 20> void
  sherman_morrison_update(
    std::vector<scalar_t> &A,
    const std::vector<scalar_t> &g) {
    const size_t p = g.size();
    std::array<scalar_t, max_size> temp_left, temp_right;
    scalar_t denominator = 1.0;
    // this is actually just a multiplication of 2 distinct vectors - Ag' and gA
    for(size_t j = 0; j < p; j++) {
      scalar_t tmp = 0.0, denom_tmp = 0.0;
      for(size_t i = 0; i < p; i++) {
        // dereference once and put on the stack - hopefully faster than 
        // using two dereferences via operator []
        scalar_t g_i = g[i];
        tmp += A[(i*p) + j] * g_i;
        denom_tmp += A[(j*p) + i] * g_i;
      }
      denominator += denom_tmp * g[j];
      // this is the first vector
      temp_left[j] = tmp;
      temp_right[j] = denom_tmp;
    }
    // this loop is only not superfluous since we do not know the denominator
    for(size_t j = 0; j < p; j++) {
      // likewise avoid extra dereferences via operator []
      const scalar_t tmp = temp_right[j];
      for(size_t i = 0; i < p; i++) {
        A[(p*j) + i] -= (temp_left[i] * tmp)/denominator;
      }
    }
}
// overload
template <typename scalar_t, const size_t max_size = 20> void
  sherman_morrison_update(
    std::vector<scalar_t> &A,
    const std::array<scalar_t, max_size> &g) {
    const size_t p = std::sqrt(A.size());
    std::array<scalar_t, max_size> temp_left, temp_right;
    scalar_t denominator = 1.0;
    // this is actually just a multiplication of 2 distinct vectors - Ag' and gA
    for(size_t j = 0; j < p; j++) {
      scalar_t tmp = 0.0, denom_tmp = 0.0;
      for(size_t i = 0; i < p; i++) {
        // dereference once and put on the stack - hopefully faster than 
        // using two dereferences via operator []
        scalar_t g_i = g[i];
        tmp += A[(i*p) + j] * g_i;
        denom_tmp += A[(j*p) + i] * g_i;
      }
      denominator += denom_tmp * g[j];
      // this is the first vector
      temp_left[j] = tmp;
      temp_right[j] = denom_tmp;
    }
    // this loop is only not superfluous since we do not know the denominator
    for(size_t j = 0; j < p; j++) {
      // likewise avoid extra dereferences via operator []
      const scalar_t tmp = temp_right[j];
      for(size_t i = 0; i < p; i++) {
        A[(p*j) + i] -= (temp_left[i] * tmp)/denominator;
      }
    }
  }


#endif
