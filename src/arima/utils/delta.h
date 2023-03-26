#ifndef DELTA
#define DELTA

template <typename scalar_t = float> std::vector<scalar_t> make_delta(
  const size_t n_diff,
  const size_t seas_period = 1,
  const size_t n_seas_diff = 0) {

  const size_t diff_size = n_diff + 1;
  std::vector<scalar_t> a(diff_size + (n_seas_diff * seas_period));
  a[0] = 1;
  std::vector<scalar_t> temp(diff_size + (n_seas_diff * seas_period));
  for (size_t k = 0; k < n_diff; k++) {
    for (size_t i = 0; i <= k; i++) {
      // the array extend is always 2, hence we can always just do these two
      // operations first this is temp[i+0] += a[i] * 1;
      temp[i] += a[i]; // * 1
      // and this is the other operation, e.g. a[i] * -1 == -= a[i];
      temp[i + 1] -= a[i];
    }
    // move all of the elements of temp to a - but temp has constant size,
    // so we can just use k+2
    for (size_t i = 0; i < k + 2; i++) {
      a[i] = std::move(temp[i]);
    }
    std::fill(temp.begin(), temp.end(), 0);
  }
  // seasonal differences:
  for (size_t k = 0; k < n_seas_diff; k++) {
    for (size_t i = 0; i < diff_size + (k * seas_period); i++) {
      /* we know that this only operates on the first and last element of
       * the vector - it adds a[i] * 1 to the first element and adds
       * a[i] * -1 to the last - which is effectively adding and subtracting
       * a[i] at various points, i.e.
       * temp[i+0] += a[i] * 1; */
      temp[i] += a[i];
      // and temp[i+seas_period] += a[i] * -1;
      temp[i + seas_period] -= a[i];
    }
    for (size_t i = 0; i < temp.size(); i++) {
      a[i] = std::move(temp[i]);
    }
    std::fill(temp.begin(), temp.end(), 0);
  }
  // remove leading coefficient and flip signs
  pop_front(a);
  for (size_t i = 0; i < a.size(); i++) {
    a[i] = -a[i];
  }
  return a;
}

#endif
