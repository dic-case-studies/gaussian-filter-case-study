#include "gaussian-filters.hpp"
#include <cstddef>

int main(int argc, char **argv) {
  if (argc != 2) {
    printf("Usage: imgauss <inimage>\n");
    abort();
  }

  char *infilename = argv[1];
  char *outfilename_expected = "out/golden_output_file.fits";
  char *outfilename_neon = "out/neon_output_file.fits";
  char *outfilename_sse = "out/sse_output_file.fits";
  char *outfilename_avx = "out/avx_output_file.fits";

  size_t size_x = 15606;
  size_t size_y = 13680;
  float *expected, *actual_neon, *actual_sse, *actual_avx;

#ifndef NGOLDEN
  expected = golden(infilename, outfilename_expected);
#endif

#if !defined(NARM_NEON) && defined(__ARM_NEON__)
  actual_neon = filter_neon(infilename, outfilename_neon);

#ifdef ASSERT
  assert_array(expected, actual_neon, size_x, size_y);
  printf("Assertions passed NEON\n");
#endif
  if (actual_neon)
    free(actual_neon);
#endif

#if !defined(NSSE) && defined(__SSE__)
  actual_sse = filter_sse(infilename, outfilename_sse);

#ifdef ASSERT
  assert_array(expected, actual_sse, size_x, size_y);
  printf("Assertions passed SSE\n");
#endif
  if (actual_sse)
    free(actual_sse);
#endif

#if !defined(NAVX2) && defined(__AVX2__)
  actual_avx = filter_avx(infilename, outfilename_avx);

#ifdef ASSERT
  assert_array(expected, actual_avx, size_x, size_y);
  printf("Assertions passed AVX\n");
#endif
  if (actual_avx)
    free(actual_avx);
#endif

  if (expected)
    free(expected);
}
