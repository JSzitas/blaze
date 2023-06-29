#ifndef BLAZE_QUANTILE
#define BLAZE_QUANTILE

#include "utils/utils.h"
#include "cmath"

template <typename scalar_t, const size_t n_pieces = 5> 
std::array<scalar_t, n_pieces> interpolate_piecewise(
  const scalar_t midpoint = 0.5, 
  const scalar_t start = 0, 
  const scalar_t end = 1) {
  static_assert(
    n_pieces % 2 != 0,
    "Piecewise interpolation not supported for even lengths of interpolated sequences"
    );
  // generate a piecewise interpolation between start and midpoint, and 
  // midpoint and end - particularly useful for generating skewed sequences
  constexpr size_t mid_size = n_pieces/2;
  // first rate of interpolation (basically the slope)
  scalar_t increment = (midpoint - start)/mid_size;
  std::array<scalar_t, n_pieces> result;
  size_t i = 0;
  for(; i < mid_size; i++) {
    result[i] = start + (i * increment);
  }
  increment = (end - midpoint)/mid_size;
  for(size_t j = 0; j < mid_size+1; j++) {
    result[i] = midpoint + (increment * j);
    i++;
  }
  return result;
}

// this in general requires O(n_landmarks * 4) storage, so it is much more 
// memory friendly than O(n) solution typically required. 
// to estimate the quantile from n samples should require something on the order 
// of O(n_landmarks * n), as each update can in worst case iterate over all 
// landmarks
template <typename scalar_t, const size_t n_landmarks = 27> class OnlineQuantile {
  scalar_t p;
  bool fitted;
  // arrays for storage of temporaries 
  std::array<scalar_t, n_landmarks> landmarks, ratio, n_at_position;
  std::array<int, n_landmarks> position;
  size_t samples_seen = 0;
  
  static_assert(
    n_landmarks % 2 != 0, 
    "Number of landmarks (template parameter n_landmarks) must be odd, not even."
    );
public:  
  OnlineQuantile<scalar_t, n_landmarks>(const scalar_t p = 0.5) : p(p), fitted(false) {
    static constexpr size_t n_minus1 = n_landmarks-1;
    static constexpr size_t n_steps = (n_minus1)/2;
    for(size_t i = 0; i < n_landmarks; i++) {
     this->landmarks[i] = 0.0; 
     this->position[i] = static_cast<int>(i+1);
    }
    this->ratio = interpolate_piecewise<scalar_t, n_landmarks>(p, 0, 1);
    this->n_at_position = interpolate_piecewise<scalar_t, n_landmarks>(p, 0, 1);
    // correct offsetting
    for(size_t i = 0; i < n_landmarks; i++) {
      this->n_at_position[i] *= (n_landmarks-1);
      this->n_at_position[i] += 1;
    }
  }
  void update(const scalar_t x) {
    if(samples_seen < n_landmarks) {
      this->landmarks[samples_seen++] = x;
      return;
    }
    samples_seen++;
    if(!this->fitted) {
      // sort landmarks seen so far
      std::sort(std::begin(this->landmarks),
                std::begin(this->landmarks) + this->samples_seen-1);
      this->fitted = true;
    }
    int k = 0;
    // smaller than smallest
    if(x < this->landmarks[0]) this->landmarks[0] = x;
    else {
      // either largest landmark or gets overwriten within the loop
      k = n_landmarks-1;
      for(size_t i = 1; i < n_landmarks-1; i++) {
        if(this->landmarks[i - 1] <= x && x < this->landmarks[i]) {
          k = i;
          break;
        }
      }
      if(this->landmarks[n_landmarks-1] < x) {
        // overwrite largest value
        this->landmarks[n_landmarks-1] = x;
      }
    }
    // increment observations at position greater than k
    for(size_t i = static_cast<size_t>(k+1); i < n_landmarks; i++) {
      this->position[i]++;
    }
    // update n available at position
    for(size_t i=0; i < n_landmarks; i++) n_at_position[i] += ratio[i];
    for(size_t i = 1; i < n_landmarks-1; i++) {
      const int n = this->position[i];
      const scalar_t q = this->landmarks[i];
      int d = this->n_at_position[i];
      d -= n;
      int pos_p1_n = this->position[i + 1] - n, pos_m1_n = this->position[i - 1] - n;
      if((d >= 1 && pos_p1_n > 1) || (d <= -1 && pos_m1_n < -1) ) {
        d = copysign(d);
        // next and previous landmarks
        scalar_t qp1 = this->landmarks[i+1], qm1 = this->landmarks[i-1];
        // and counts at those positions
        int np1 = this->position[i+1], nm1 = this->position[i-1];
        // interpolating step
        const scalar_t qn = compute_p2(qp1, q, qm1, d, np1, n, nm1);
        // std::cout << "Previous: " << this->landmarks[i] << " Index: " << i;
        if(qm1 < qn && qn < qp1) {
          // std::cout << " | At midpoint |";
          this->landmarks[i] = qn;
        } else {
          // std::cout << " | Not at midpoint - linear interpolation | " << std::endl;
          this->landmarks[i] = q + d * (this->landmarks[i + d] - q) /
            (this->position[i + d] - n);
        }
        this->position[i] = n + d;
        // std::cout << " Updated: " << this->landmarks[i] << std::endl;
      }
    }
  }
  const scalar_t quantile() {
    // return midpoint landmark
    if(this->fitted) return this->landmarks[n_landmarks/2];
    return std::nan("0.0"); 
  }
  void print_summary() {
    std::cout << "Landmarks: ";
    print_vector(this->landmarks);
    std::cout << "N at position: ";
    print_vector(this->n_at_position);
    std::cout << "Position: ";
    print_vector(this->position);
    std::cout << "Ratio: ";
    print_vector(this->ratio);
    std::cout << "p: " << this->p << std::endl;
  }
private:
  scalar_t compute_p2(const scalar_t qp1, const scalar_t q, const scalar_t qm1,
                      const scalar_t d, const int np1, const scalar_t n, const int nm1) {
    const scalar_t range = d/(np1 - nm1),
      left = (n - nm1 + d) * (qp1 - q )/(np1 - n),
      right = (np1 - n - d) * (q - qm1 )/(n - nm1);
    return q + range * (left + right);
  }
};
// functional wrapper for memory savings
template <typename scalar_t> scalar_t quantile(
  const std::vector<scalar_t> &x,
  const scalar_t p = 0.5) {
  std::vector<scalar_t> y = x;
  std::sort(y.begin(), y.end());
  const size_t index = p * (x.size()-1);
  return y[index];
}

#endif                                
                                