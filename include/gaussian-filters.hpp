#ifndef GAUSSIAN_FILTERS_HPP
#define GAUSSIAN_FILTERS_HPP

#include "common.h"
#include "fits-helper.h"
#include "reference.h"
#include <time.h>

#ifdef __x86_64__
#include <immintrin.h>
#endif

#ifdef __aarch64__
#include "sse2neon.h"
#include <arm_neon.h>
#endif

float *filter_golden(const char *infilename, const char *outfilename,
                     size_t &size_x, size_t &size_y) {
  int status = 0;

  size_t n_iter;
  size_t filter_radius;
  double sigma = 3.5;
  optimal_filter_size_dbl(sigma, &filter_radius, &n_iter);
  printf("sigma: %f filter_radius: %ld niter: %ld\n", sigma, filter_radius,
         n_iter);

  fitsfile *fptr;
  fits_data_t fits;
  if (!fits_open_image(&fptr, infilename, READONLY, &status)) {
    status = extract_data_from_fits(fptr, &fits);
  } else {
    fits_report_error(stderr, status);
    abort();
  }

  float *data = fits.data;
  size_x = fits.naxes[0];
  size_y = fits.naxes[1];

  float *column = (float *)malloc(sizeof(float) * size_y);
  float *data_col =
      (float *)malloc(sizeof(float) * (size_y + 2 * filter_radius));
  float *data_row =
      (float *)malloc(sizeof(float) * (size_x + 2 * filter_radius));

  clock_t begin = clock();

  filter_gauss_2d_flt(data, column, data_row, data_col, size_x, size_y, 1,
                      filter_radius);

  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("Time elapsed GOLDEN: %f s\n", time_spent);

  printf("Storing result in FITS\n");
  status = store_data_copying_fits_header(fptr, outfilename, &fits);

  // destroy_fits(&fits);
  fits_close_file(fptr, &status);

  if (status) {
    fits_report_error(stderr, status);
  }

  free(column);
  free(data_col);
  free(data_row);

  return fits.data;
}

#if !defined(NARM_NEON) && defined(__ARM_NEON__)
void filter_gauss_2d_neon(float *data, const size_t size_x, const size_t size_y,
                          const size_t n_iter, const size_t filter_radius) {
  // Set up a few variables
  const size_t size_xy = size_x * size_y;
  float *ptr = data + size_xy;

  float *data_row =
      (float *)malloc(sizeof(float) * (size_x + 2 * filter_radius));
  // Run row filter (along x-axis)
  // This is straightforward, as the data are contiguous in x.
  while (ptr > data) {
    ptr -= size_x;
    // for (size_t i = n_iter; i--;)
    filter_boxcar_1d_flt(ptr, data_row, size_x, filter_radius);
  }

  free(data_row);

  size_t quot = size_x / 4;
  size_t limit = quot * 4;
  for (size_t x = 0; x < limit; x += 4) {
    // Apply filter
    for (size_t i = n_iter; i--;)
      filter_simd_neon(data + x, size_y, size_x, filter_radius);
  }

  float *ptr2;
  float *data_copy = (float *)malloc(sizeof(float) * size_y);
  float *data_col =
      (float *)malloc(sizeof(float) * (size_y + 2 * filter_radius));

  for (size_t x = size_x; x > limit; x--) {
    // Copy data into column array
    ptr = data + size_xy - size_x + x - 1;
    ptr2 = data_copy + size_y;
    while (ptr2-- > data_copy) {
      *ptr2 = *ptr;
      ptr -= size_x;
    }

    // Apply filter
    // for (size_t i = n_iter; i--;)
    filter_boxcar_1d_flt(data_copy, data_col, size_y, filter_radius);

    // Copy column array back into data array
    ptr = data + size_xy - size_x + x - 1;
    ptr2 = data_copy + size_y;
    while (ptr2-- > data_copy) {
      *ptr = *ptr2;
      ptr -= size_x;
    }
  }

  free(data_copy);
  free(data_col);

  return;
}

float *filter_neon(const char *infilename, const char *outfilename) {
  int status = 0;

  size_t n_iter;
  size_t filter_radius;
  double sigma = 3.5;
  optimal_filter_size_dbl(sigma, &filter_radius, &n_iter);
  printf("sigma: %f filter_radius: %ld niter: %ld\n", sigma, filter_radius,
         n_iter);

  fitsfile *fptr;
  fits_data_t fits;
  if (!fits_open_image(&fptr, infilename, READONLY, &status)) {
    status = extract_data_from_fits(fptr, &fits);
  } else {
    fits_report_error(stderr, status);
    abort();
  }

  float *data = fits.data;
  size_t size_x = fits.naxes[0];
  size_t size_y = fits.naxes[1];

  clock_t begin = clock();

  filter_gauss_2d_neon(data, size_x, size_y, 1, filter_radius);

  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("Time elapsed NEON: %f s\n", time_spent);

  printf("Storing result in FITS\n");
  status = store_data_copying_fits_header(fptr, outfilename, &fits);

  // destroy_fits(&fits);
  fits_close_file(fptr, &status);

  if (status) {
    fits_report_error(stderr, status);
  }
  return fits.data;
}
#endif

#if !defined(NSSE) && defined(__SSE__)
void filter_gauss_2d_sse(float *data, const size_t size_x, const size_t size_y,
                         const size_t n_iter, const size_t filter_radius) {
  // Set up a few variables
  const size_t size_xy = size_x * size_y;
  float *ptr = data + size_xy;

  float *data_row = (float *)_mm_malloc(
      (size_x + 2 * filter_radius) * sizeof(float), sizeof(float));

  // Run row filter (along x-axis)
  // This is straightforward, as the data are contiguous in x.
  while (ptr > data) {
    ptr -= size_x;
    // for (size_t i = n_iter; i--;)
    filter_boxcar_1d_flt(ptr, data_row, size_x, filter_radius);
  }

  free(data_row);

  size_t quot = size_x / 4;
  size_t limit = quot * 4;

  for (size_t x = 0; x < limit; x += 4) {
    // Apply filter
    for (size_t i = n_iter; i--;)
      filter_simd_sse(data + x, size_y, size_x, filter_radius);
  }

  float *ptr2;
  // data_copy = (float *)malloc(sizeof(float) * size_y);
  float *data_copy = (float *)_mm_malloc(size_y * sizeof(float), sizeof(float));

  // data_col = (float *)malloc(sizeof(float) * (size_y + 2 * filter_radius));
  float *data_col = (float *)_mm_malloc(
      (size_y + 2 * filter_radius) * sizeof(float), sizeof(float));

  for (size_t x = size_x; x > limit; x--) {
    // Copy data into column array
    ptr = data + size_xy - size_x + x - 1;
    ptr2 = data_copy + size_y;
    while (ptr2-- > data_copy) {
      *ptr2 = *ptr;
      ptr -= size_x;
    }

    // Apply filter
    // for (size_t i = n_iter; i--;)
    filter_boxcar_1d_flt(data_copy, data_col, size_y, filter_radius);

    // Copy column array back into data array
    ptr = data + size_xy - size_x + x - 1;
    ptr2 = data_copy + size_y;
    while (ptr2-- > data_copy) {
      *ptr = *ptr2;
      ptr -= size_x;
    }
  }

  free(data_copy);
  free(data_col);

  return;
}

float *filter_sse(const char *infilename, const char *outfilename) {
  int status = 0;

  size_t n_iter;
  size_t filter_radius;
  double sigma = 3.5;
  optimal_filter_size_dbl(sigma, &filter_radius, &n_iter);
  printf("sigma: %f filter_radius: %ld niter: %ld\n", sigma, filter_radius,
         n_iter);

  fitsfile *fptr;
  fits_data_t fits;
  if (!fits_open_image(&fptr, infilename, READONLY, &status)) {
    status = extract_data_from_fits(fptr, &fits);
  } else {
    fits_report_error(stderr, status);
    abort();
  }

  float *data = fits.data;
  size_t size_x = fits.naxes[0];
  size_t size_y = fits.naxes[1];

  clock_t begin = clock();

  filter_gauss_2d_sse(data, size_x, size_y, 1, filter_radius);

  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("Time elapsed SSE: %f s\n", time_spent);

  printf("Storing result in FITS\n");
  status = store_data_copying_fits_header(fptr, outfilename, &fits);

  // destroy_fits(&fits);
  fits_close_file(fptr, &status);

  if (status) {
    fits_report_error(stderr, status);
  }
  return fits.data;
}
#endif

#if !defined(NAVX2) && defined(__AVX2__)
void filter_gauss_2d_avx(float *data, const size_t size_x, const size_t size_y,
                         const size_t n_iter, const size_t filter_radius) {
  // Set up a few variables
  const size_t size_xy = size_x * size_y;
  float *ptr = data + size_xy;

  float *data_row = (float *)_mm_malloc(
      (size_x + 2 * filter_radius) * sizeof(float), sizeof(float));

  // Run row filter (along x-axis)
  // This is straightforward, as the data are contiguous in x.
  while (ptr > data) {
    ptr -= size_x;
    // for (size_t i = n_iter; i--;)
    filter_boxcar_1d_flt(ptr, data_row, size_x, filter_radius);
  }

  free(data_row);

  size_t quot = size_x / 8;
  size_t limit = quot * 8;

  for (size_t x = 0; x < limit; x += 8) {
    // Apply filter
    for (size_t i = n_iter; i--;)
      filter_simd_avx(data + x, size_y, size_x, filter_radius);
  }

  float *ptr2;
  // data_copy = (float *)malloc(sizeof(float) * size_y);
  float *data_copy = (float *)_mm_malloc(size_y * sizeof(float), sizeof(float));
  // data_col = (float *)malloc(sizeof(float) * (size_y + 2 * filter_radius));
  float *data_col = (float *)_mm_malloc(
      (size_y + 2 * filter_radius) * sizeof(float), sizeof(float));

  for (size_t x = size_x; x > limit; x--) {
    // Copy data into column array
    ptr = data + size_xy - size_x + x - 1;
    ptr2 = data_copy + size_y;
    while (ptr2-- > data_copy) {
      *ptr2 = *ptr;
      ptr -= size_x;
    }

    // Apply filter
    // for (size_t i = n_iter; i--;)
    filter_boxcar_1d_flt(data_copy, data_col, size_y, filter_radius);

    // Copy column array back into data array
    ptr = data + size_xy - size_x + x - 1;
    ptr2 = data_copy + size_y;
    while (ptr2-- > data_copy) {
      *ptr = *ptr2;
      ptr -= size_x;
    }
  }

  free(data_copy);
  free(data_col);

  return;
}

float *filter_avx(const char *infilename, const char *outfilename) {
  int status = 0;

  size_t n_iter;
  size_t filter_radius;
  double sigma = 3.5;
  optimal_filter_size_dbl(sigma, &filter_radius, &n_iter);
  printf("sigma: %f filter_radius: %ld niter: %ld\n", sigma, filter_radius,
         n_iter);

  fitsfile *fptr;
  fits_data_t fits;
  if (!fits_open_image(&fptr, infilename, READONLY, &status)) {
    status = extract_data_from_fits(fptr, &fits);
  } else {
    fits_report_error(stderr, status);
    abort();
  }

  float *data = fits.data;
  size_t size_x = fits.naxes[0];
  size_t size_y = fits.naxes[1];

  clock_t begin = clock();

  filter_gauss_2d_avx(data, size_x, size_y, 1, filter_radius);

  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("Time elapsed AVX: %f s\n", time_spent);

  printf("Storing result in FITS\n");
  status = store_data_copying_fits_header(fptr, outfilename, &fits);

  // destroy_fits(&fits);
  fits_close_file(fptr, &status);

  if (status) {
    fits_report_error(stderr, status);
  }
  return fits.data;
}
#endif

#if !defined(NAVX2) && defined(__AVX512F__)
void filter_gauss_2d_avx_512(float *data, const size_t size_x,
                             const size_t size_y, const size_t n_iter,
                             const size_t filter_radius) {
  // Set up a few variables
  const size_t size_xy = size_x * size_y;
  float *ptr = data + size_xy;

  float *data_row = (float *)_mm_malloc(
      (size_x + 2 * filter_radius) * sizeof(float), sizeof(float));

  // Run row filter (along x-axis)
  // This is straightforward, as the data are contiguous in x.
  while (ptr > data) {
    ptr -= size_x;
    // for (size_t i = n_iter; i--;)
    filter_boxcar_1d_flt(ptr, data_row, size_x, filter_radius);
  }

  free(data_row);

  size_t quot = size_x / 16;
  size_t limit = quot * 16;

  for (size_t x = 0; x < limit; x += 16) {
    // Apply filter
    // for (size_t i = n_iter; i--;) {
    filter_simd_avx_512(data + x, size_y, size_x, filter_radius);
    // }
  }

  float *ptr2;
  // data_copy = (float *)malloc(sizeof(float) * size_y);
  float *data_copy = (float *)_mm_malloc(size_y * sizeof(float), sizeof(float));
  // data_col = (float *)malloc(sizeof(float) * (size_y + 2 * filter_radius));
  float *data_col = (float *)_mm_malloc((size_y + 2 * filter_radius) * sizeof(float),
                                 sizeof(float));

  for (size_t x = size_x; x > limit; x--) {
    // Copy data into column array
    ptr = data + size_xy - size_x + x - 1;
    ptr2 = data_copy + size_y;
    while (ptr2-- > data_copy) {
      *ptr2 = *ptr;
      ptr -= size_x;
    }

    // Apply filter
    // for (size_t i = n_iter; i--;)
    filter_boxcar_1d_flt(data_copy, data_col, size_y, filter_radius);

    // Copy column array back into data array
    ptr = data + size_xy - size_x + x - 1;
    ptr2 = data_copy + size_y;
    while (ptr2-- > data_copy) {
      *ptr = *ptr2;
      ptr -= size_x;
    }
  }

  free(data_copy);
  free(data_col);

  return;
}

float *filter_avx_512(const char *infilename, const char *outfilename) {
  int status = 0;

  size_t n_iter;
  size_t filter_radius;
  double sigma = 3.5;
  optimal_filter_size_dbl(sigma, &filter_radius, &n_iter);
  printf("sigma: %f filter_radius: %ld niter: %ld\n", sigma, filter_radius,
         n_iter);

  fitsfile *fptr;
  fits_data_t fits;
  if (!fits_open_image(&fptr, infilename, READONLY, &status)) {
    status = extract_data_from_fits(fptr, &fits);
  } else {
    fits_report_error(stderr, status);
    abort();
  }

  float *data = fits.data;
  size_t size_x = fits.naxes[0];
  size_t size_y = fits.naxes[1];

  clock_t begin = clock();

  filter_gauss_2d_avx_512(data, size_x, size_y, n_iter, filter_radius);

  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("Time elapsed AVX512: %f s\n", time_spent);

  printf("Storing result in FITS\n");
  status = store_data_copying_fits_header(fptr, outfilename, &fits);

  // destroy_fits(&fits);
  fits_close_file(fptr, &status);

  if (status) {
    fits_report_error(stderr, status);
  }
  return fits.data;
}
#endif

#endif // GAUSSIAN_FILTERS_HPP