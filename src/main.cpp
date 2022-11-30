#include <cstddef>
#include "gaussian-filters.hpp"

int main(int argc, char **argv)
{
  if (argc != 2)
  {
    printf("Usage: imgauss <inimage>\n");
    abort();
  }

  char *infilename = argv[1];
  char *outfilename_expected = "golden_output_file.fits";
  char *outfilename_neon = "neon_output_file.fits";
  char *outfilename_sse = "sse_output_file.fits";
  char *outfilename_avx = "avx_output_file.fits";

  size_t size_x = 1024;
  size_t size_y = 1024;
  float *expected, *actual_neon, *actual_sse, *actual_avx;

  expected = golden(infilename, outfilename_expected);

#ifdef __ARM_NEON__
  actual_neon = filter_neon(infilename, outfilename_neon);

  assert_array(expected, actual_neon, size_x, size_y);
  printf("Assertions passed NEON\n");
#endif

#ifdef __SSE__
  actual_sse = filter_sse(infilename, outfilename_sse);

  assert_array(expected, actual_sse, size_x, size_y);
  printf("Assertions passed SSE\n");
#endif

#ifdef __AVX2__
  actual_avx = filter_avx(infilename, outfilename_avx);

  assert_array(expected, actual_avx, size_x, size_y);
  printf("Assertions passed AVX\n");
#endif

  if (expected)
    free(expected);
  if (actual_neon)
    free(actual_neon);
  if (actual_sse)
    free(actual_sse);
  if (actual_avx)
    free(actual_avx);
}
