#pragma once
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef __x86_64__
#include <immintrin.h>
#endif

#ifdef __aarch64__
#include "sse2neon.h"
#include <arm_neon.h>
#endif

#define FILTER_NAN(x) ((x) == (x) ? (x) : 0)
#define BOXCAR_MIN_ITER 3
#define BOXCAR_MAX_ITER 6
#define TOLERANCE 0.00001

void optimal_filter_size_dbl(const double sigma, size_t *filter_radius,
                             size_t *n_iter) {
  *n_iter = 0;
  *filter_radius = 0;
  double tmp = -1.0;
  size_t i;

  for (i = BOXCAR_MIN_ITER; i <= BOXCAR_MAX_ITER; ++i) {
    const double radius = sqrt((3.0 * sigma * sigma / i) + 0.25) - 0.5;
    const double diff = fabs(radius - floor(radius + 0.5));

    if (tmp < 0.0 || diff < tmp) {
      tmp = diff;
      *n_iter = i;
      *filter_radius = (size_t)(radius + 0.5);
    }
  }
  return;
}

void filter_boxcar_1d_flt(float *data, float *data_copy, const size_t size,
                          const size_t filter_radius) {
  // Define filter size
  const size_t filter_size = 2 * filter_radius + 1;
  const float inv_filter_size = 1.0 / filter_size;
  size_t i;

  // Make copy of data, taking care of NaN
  for (i = size; i--;)
    data_copy[filter_radius + i] = FILTER_NAN(data[i]);

  // Fill overlap regions with 0
  for (i = filter_radius; i--;)
    data_copy[i] = data_copy[size + filter_radius + i] = 0.0;

  // Apply boxcar filter to last data point
  data[size - 1] = 0.0;
  for (i = filter_size; i--;)
    data[size - 1] += data_copy[size + i - 1];
  data[size - 1] *= inv_filter_size;

  // Recursively apply boxcar filter to all previous data points
  for (i = size - 1; i--;)
    data[i] = data[i + 1] +
              (data_copy[i] - data_copy[filter_size + i]) * inv_filter_size;

  return;
}

#if !defined(NSSE) && defined(__SSE__)
static inline __m128 filter_nan_sse(__m128 data_v) {
  // Filter nan
  __m128 nan_mask = _mm_cmpord_ps(data_v, data_v);
  __m128 zero_v = _mm_setzero_ps();
  __m128 data_filtered_v = _mm_blendv_ps(zero_v, data_v, nan_mask);
  return data_filtered_v;
}

void print_vector(__m128 vec) {
  float res[4];
  _mm_storeu_ps(res, vec);
  for (int i = 0; i < 4; i++) {
    printf("%.1f ", res[i]);
  }
  printf("\n");
}

void filter_simd_sse(float *data, const size_t size, const size_t stride,
                     const size_t filter_radius) {
  const size_t filter_size = 2 * filter_radius + 1;
  size_t i;

  const float inv_filter_size = 1.0 / filter_size;
  __m128 inv_filter_size_v = _mm_set1_ps(inv_filter_size);
  __m128 zero_v = _mm_setzero_ps();

  // __m128 *data_copy = (__m128 *)malloc(sizeof(__m128) * (size + 2 *
  // filter_radius));

  __m128 *data_copy = (__m128 *)_mm_malloc(
      sizeof(__m128) * (size + 2 * filter_radius), sizeof(__m128));

  for (i = size; i--;)
    data_copy[filter_radius + i] =
        filter_nan_sse(_mm_loadu_ps(data + (stride * i)));

  for (i = filter_radius; i--;)
    data_copy[i] = data_copy[size + filter_radius + i] = zero_v;

  // Calculate last point
  __m128 last_pt = zero_v;
  for (i = filter_size; i--;) {
    last_pt = _mm_add_ps(last_pt, data_copy[size + i - 1]);
  }
  last_pt = _mm_mul_ps(last_pt, inv_filter_size_v);
  _mm_storeu_ps(data + (stride * (size - 1)), last_pt);

  __m128 next_pt = last_pt;
  for (int col = size - 1; col--;) {
    __m128 current_pt =
        _mm_sub_ps(data_copy[col], data_copy[filter_size + col]);
    current_pt = _mm_mul_ps(current_pt, inv_filter_size_v);
    current_pt = _mm_add_ps(next_pt, current_pt);

    _mm_storeu_ps(data + (stride * col), current_pt);

    next_pt = current_pt;
  }

  free(data_copy);

  return;
}
#endif

#if !defined(NARM_NEON) && defined(__ARM_NEON__)
static inline float32x4_t filter_nan_neon(float32x4_t data_v) {
  // Filter nan
  uint32x4_t nan_mask = vceqq_f32(data_v, data_v);
  float32x4_t zero_v = vdupq_n_f32(0);
  float32x4_t data_filtered_v = vbslq_f32(nan_mask, data_v, zero_v);
  return data_filtered_v;
}

void filter_simd_neon(float *data, const size_t size, const size_t stride,
                      const size_t filter_radius) {
  const size_t filter_size = 2 * filter_radius + 1;
  size_t i;

  const float inv_filter_size = 1.0 / filter_size;
  float32x4_t inv_filter_size_v = vdupq_n_f32(inv_filter_size);

  float32x4_t zero_v = vdupq_n_f32(0);

  float32x4_t *data_copy =
      (float32x4_t *)malloc(sizeof(float32x4_t) * (size + 2 * filter_radius));

  for (i = size; i--;)
    data_copy[filter_radius + i] =
        filter_nan_neon(vld1q_f32(data + (stride * i)));

  for (i = filter_radius; i--;)
    data_copy[i] = data_copy[size + filter_radius + i] = zero_v;

  // Calculate last point
  float32x4_t last_pt = zero_v;
  for (i = filter_size; i--;) {
    last_pt = vaddq_f32(last_pt, data_copy[size + i - 1]);
  }
  last_pt = vmulq_f32(last_pt, inv_filter_size_v);
  vst1q_f32(data + (stride * (size - 1)), last_pt);

  float32x4_t next_pt = last_pt;
  for (int col = size - 1; col--;) {
    float32x4_t current_pt =
        vsubq_f32(data_copy[col], data_copy[filter_size + col]);
    current_pt = vmulq_f32(current_pt, inv_filter_size_v);
    current_pt = vaddq_f32(next_pt, current_pt);

    vst1q_f32(data + (stride * col), current_pt);

    next_pt = current_pt;
  }

  free(data_copy);

  return;
}
#endif

#if !defined(NAVX2) && defined(__AVX2__)
static inline __m256 filter_nan_avx(__m256 data_v) {
  // Filter nan
  __m256 nan_mask = _mm256_cmp_ps(data_v, data_v, _CMP_ORD_Q);
  __m256 zero_v = _mm256_setzero_ps();
  __m256 data_filtered_v = _mm256_blendv_ps(zero_v, data_v, nan_mask);
  return data_filtered_v;
}

void filter_simd_avx(float *data, const size_t size, const size_t stride,
                     const size_t filter_radius) {
  const size_t filter_size = 2 * filter_radius + 1;
  size_t i;

  const float inv_filter_size = 1.0 / filter_size;
  __m256 inv_filter_size_v = _mm256_set1_ps(inv_filter_size);
  __m256 zero_v = _mm256_setzero_ps();

  // __m256 *data_copy = (__m256 *)malloc(sizeof(__m256) * (size + 2 *
  // filter_radius));

  __m256 *data_copy = (__m256 *)_mm_malloc(
      sizeof(__m256) * (size + 2 * filter_radius), sizeof(__m256));

  for (i = size; i--;)
    data_copy[filter_radius + i] =
        filter_nan_avx(_mm256_loadu_ps(data + (stride * i)));

  for (i = filter_radius; i--;)
    data_copy[i] = data_copy[size + filter_radius + i] = zero_v;

  // Calculate last point
  __m256 last_pt = zero_v;
  for (i = filter_size; i--;) {
    last_pt = _mm256_add_ps(last_pt, data_copy[size + i - 1]);
  }
  last_pt = _mm256_mul_ps(last_pt, inv_filter_size_v);
  _mm256_storeu_ps(data + (stride * (size - 1)), last_pt);

  __m256 next_pt = last_pt;
  for (int col = size - 1; col--;) {
    __m256 current_pt =
        _mm256_sub_ps(data_copy[col], data_copy[filter_size + col]);
    current_pt = _mm256_mul_ps(current_pt, inv_filter_size_v);
    current_pt = _mm256_add_ps(next_pt, current_pt);

    _mm256_storeu_ps(data + (stride * col), current_pt);

    next_pt = current_pt;
  }

  free(data_copy);

  return;
}
#endif


#if !defined(NAVX2) && defined(__AVX512F__) 
static inline __m512 filter_nan_avx_512(__m512 data_v) {
  // Filter nan
  __mmask16 nan_mask = _mm512_cmp_ps_mask(data_v, data_v, _CMP_ORD_Q);
  __m512 zero_v = _mm512_setzero_ps();
  __m512 data_filtered_v = _mm512_mask_blend_ps(nan_mask, zero_v, data_v);
  return data_filtered_v;
}

void filter_simd_avx_512(float *data, const size_t size, const size_t stride,
                     const size_t filter_radius) {
  const size_t filter_size = 2 * filter_radius + 1;
  size_t i;

  const float inv_filter_size = 1.0 / filter_size;
  __m512 inv_filter_size_v = _mm512_set1_ps(inv_filter_size);
  __m512 zero_v = _mm512_setzero_ps();

  // __m256 *data_copy = (__m256 *)malloc(sizeof(__m256) * (size + 2 *
  // filter_radius));

  __m512 *data_copy = (__m512 *)_mm_malloc(
      sizeof(__m512) * (size + 2 * filter_radius), sizeof(__m512));

  for (i = size; i--;)
    data_copy[filter_radius + i] =
        filter_nan_avx_512(_mm512_loadu_ps(data + (stride * i)));

  for (i = filter_radius; i--;)
    data_copy[i] = data_copy[size + filter_radius + i] = zero_v;

  // Calculate last point
  __m512 last_pt = zero_v;
  for (i = filter_size; i--;) {
    last_pt = _mm512_add_ps(last_pt, data_copy[size + i - 1]);
  }
  last_pt = _mm512_mul_ps(last_pt, inv_filter_size_v);
  _mm512_storeu_ps(data + (stride * (size - 1)), last_pt);

  __m512 next_pt = last_pt;
  for (int col = size - 1; col--;) {
    __m512 current_pt =
        _mm512_sub_ps(data_copy[col], data_copy[filter_size + col]);
    current_pt = _mm512_mul_ps(current_pt, inv_filter_size_v);
    current_pt = _mm512_add_ps(next_pt, current_pt);

    _mm512_storeu_ps(data + (stride * col), current_pt);

    next_pt = current_pt;
  }

  free(data_copy);

  return;
}
#endif

void assert_array(float *expected, float *actual, size_t size_x,
                  size_t size_y) {
  for (size_t y = 0; y < size_y; y++) {
    for (size_t x = 0; x < size_x; x++) {
      float expected_value = expected[y * size_x + x];
      float actual_value = actual[y * size_x + x];

      float diff = fabs(expected_value - actual_value);
      if (diff > TOLERANCE) {
        printf("Error asserting point %ld %ld, expected: %f actual %f\n", x, y,
               expected_value, actual_value);
        exit(-1);
        // assert(deviation < TOLERANCE);
      }
    }
  }
}
