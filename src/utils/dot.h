#ifndef SIMD_AR_DOT_HEADER
#define SIMD_AR_DOT_HEADER

template<typename T>
inline T dot(const T* x, const T* y, int f) {
  T s = 0;
  for (int z = 0; z < f; z++) {
    s += (*x) * (*y);
    x++;
    y++;
  }
  return s;
}

// inspired by annoylib, see https://github.com/spotify/annoy/blob/main/src/annoylib.h
// #if !defined(NO_MANUAL_VECTORIZATION) && defined(__GNUC__) && (__GNUC__ >6) && defined(__AVX512F__) 
// #define DOT_USE_AVX512
#if !defined(NO_MANUAL_VECTORIZATION) && defined(__AVX__) && defined (__SSE__) && defined(__SSE2__) && defined(__SSE3__)
#define DOT_USE_AVX
#else
#endif

#if defined(DOT_USE_AVX) || defined(DOT_USE_AVX512)
#if defined(_MSC_VER)
#include <intrin.h>
#elif defined(__GNUC__)
#include <x86intrin.h>
#endif
#endif

#ifdef DOT_USE_AVX
// Horizontal single sum of 256bit vector.
inline float hsum256_ps_avx(__m256 v) {
  const __m128 x128 = _mm_add_ps(_mm256_extractf128_ps(v, 1), _mm256_castps256_ps128(v));
  const __m128 x64 = _mm_add_ps(x128, _mm_movehl_ps(x128, x128));
  const __m128 x32 = _mm_add_ss(x64, _mm_shuffle_ps(x64, x64, 0x55));
  return _mm_cvtss_f32(x32);
}

inline double hsum256_pd_avx(__m256d v) {
  __m128d vlow  = _mm256_castpd256_pd128(v);
  __m128d vhigh = _mm256_extractf128_pd(v, 1); // high 128
  vlow  = _mm_add_pd(vlow, vhigh);     // reduce down to 128
  __m128d high64 = _mm_unpackhi_pd(vlow, vlow);
  return  _mm_cvtsd_f64(_mm_add_sd(vlow, high64));  // reduce to scalar
}

// overload
template<>
inline float dot<float>(const float* x, const float *y, int f) {
  float result = 0;
  if (f > 7) {
    __m256 d = _mm256_setzero_ps();
    for (; f > 7; f -= 8) {
      d = _mm256_add_ps(d, _mm256_mul_ps(_mm256_loadu_ps(x), _mm256_loadu_ps(y)));
      x += 8;
      y += 8;
    }
    // Sum all floats in dot register.
    result += hsum256_ps_avx(d);
  }
  // Don't forget the remaining values.
  for (; f > 0; f--) {
    result += *x * *y;
    x++;
    y++;
  }
  return result;
}

// second overload
template<>
inline double dot<double>(const double* x, const double *y, int f) {
  double result = 0;
  if (f > 3) {
    __m256d d = _mm256_setzero_pd();
    for (; f > 3; f -= 4) {
      d = _mm256_add_pd(d, _mm256_mul_pd(_mm256_loadu_pd(x), _mm256_loadu_pd(y)));
      x += 4;
      y += 4;
    }
    // Sum all floats in dot register.
    result += hsum256_pd_avx(d);
  }
  // Don't forget the remaining values.
  for (; f > 0; f--) {
    result += *x * *y;
    x++;
    y++;
  }
  return result;
}
#endif

#endif
