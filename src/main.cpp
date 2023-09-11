#include "gaussian-filters.hpp"
#include <cstddef>
#include <string>

int main(int argc, char **argv) {
  if (argc != 2) {
    printf("Usage: %s <inimage>\n", argv[0]);
    abort();
  }

  std::string infilename = argv[1];
  std::string outfilename_expected = "out/golden_output_file.fits";
  std::string outfilename_neon = "out/neon_output_file.fits";
  std::string outfilename_sse = "out/sse_output_file.fits";
  std::string outfilename_avx = "out/avx_output_file.fits";
  std::string outfilename_avx_512 = "out/avx_512_output_file.fits";

  size_t size_x;
  size_t size_y;

#ifndef NGOLDEN
  float *expected;
  expected =
      filter_golden(infilename.c_str(), outfilename_expected.c_str(), size_x, size_y);
#endif

#if !defined(NARM_NEON) && defined(__ARM_NEON__)
  float *actual_neon;
  actual_neon = filter_neon(infilename.c_str(), outfilename_neon.c_str());

#ifdef ASSERT
  assert_array(expected, actual_neon, size_x, size_y);
  printf("Assertions passed NEON\n");
#endif
  if (actual_neon)
    free(actual_neon);
#endif

#if !defined(NSSE) && defined(__SSE__)
  float *actual_sse;
  actual_sse = filter_sse(infilename.c_str(), outfilename_sse.c_str());

#ifdef ASSERT
  assert_array(expected, actual_sse, size_x, size_y);
  printf("Assertions passed SSE\n");
#endif
  if (actual_sse)
    free(actual_sse);
#endif

#if !defined(NAVX2) && defined(__AVX2__)
  float *actual_avx;
  actual_avx = filter_avx(infilename.c_str(), outfilename_avx.c_str());

#ifdef ASSERT
  assert_array(expected, actual_avx, size_x, size_y);
  printf("Assertions passed AVX\n");
#endif
  if (actual_avx)
    free(actual_avx);
#endif

// find a more generic AVX512 macro for this
#if !defined(NAVX512) && defined(__AVX512F__)
  float *actual_avx_512;
  actual_avx_512 =
      filter_avx_512(infilename.c_str(), outfilename_avx_512.c_str());

#ifdef ASSERT
  assert_array(expected, actual_avx_512, size_x, size_y);
  printf("Assertions passed AVX512\n");
#endif
  if (actual_avx_512)
    free(actual_avx_512);
#endif

  if (expected)
    free(expected);
}
