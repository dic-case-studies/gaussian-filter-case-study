#pragma once
#include "common.h"
#include <stdlib.h>

/* Code taken from SoFiA https://github.com/SoFiA-Admin/SoFiA-2.git */

#define FILTER_NAN(x) ((x) == (x) ? (x) : 0)
#define BOXCAR_MIN_ITER 3
#define BOXCAR_MAX_ITER 6

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

void filter_gauss_2d_flt(float *data, float *data_copy, float *data_row,
                         float *data_col, const size_t size_x,
                         const size_t size_y, const size_t n_iter,
                         const size_t filter_radius) {
  // Set up a few variables
  const size_t size_xy = size_x * size_y;
  float *ptr = data + size_xy;
  float *ptr2;

  // Run row filter (along x-axis)
  // This is straightforward, as the data are contiguous in x.
  while (ptr > data) {
    ptr -= size_x;
    for (size_t i = n_iter; i--;)
      filter_boxcar_1d_flt(ptr, data_row, size_x, filter_radius);
  }

  // Run column filter (along y-axis)
  // This is more complicated, as the data are non-contiguous in y.
  for (size_t x = size_x; x--;) {
    // Copy data into column array
    ptr = data + size_xy - size_x + x;
    ptr2 = data_copy + size_y;
    while (ptr2-- > data_copy) {
      *ptr2 = *ptr;
      ptr -= size_x;
    }

    // Apply filter
    for (size_t i = n_iter; i--;)
      filter_boxcar_1d_flt(data_copy, data_col, size_y, filter_radius);

    // Copy column array back into data array
    ptr = data + size_xy - size_x + x;
    ptr2 = data_copy + size_y;
    while (ptr2-- > data_copy) {
      *ptr = *ptr2;
      ptr -= size_x;
    }
  }

  return;
}
