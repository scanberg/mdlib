INCLUDE(CheckCSourceRuns)

SET(HAVE_AVX_EXTENSIONS)
SET(HAVE_AVX2_EXTENSIONS)
SET(HAVE_FMA_EXTENSIONS)

# Check AVX 2
SET(CMAKE_REQUIRED_FLAGS)
IF(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  SET(CMAKE_REQUIRED_FLAGS "-mavx2 -mfma")
ELSEIF(MSVC AND NOT CMAKE_CL_64)  # reserve for WINDOWS
  SET(CMAKE_REQUIRED_FLAGS "/arch:AVX2")
ENDIF()

CHECK_C_SOURCE_RUNS("
#include <immintrin.h>
int main()
{
    __m256i a = _mm256_set_epi32 (-1, 2, -3, 4, -1, 2, -3, 4);
    __m256i result = _mm256_abs_epi32 (a);
    return 0;
}" HAVE_AVX2_EXTENSIONS)

CHECK_C_SOURCE_RUNS("
#include <immintrin.h>
int main()
{
    __m256 a = _mm256_set1_ps (1.0f);
    __m256 b = _mm256_set1_ps (2.0f);
    __m256 c = _mm256_set1_ps (3.0f);
    __m256 result = _mm256_fmadd_ps(a, b, c);
    return 0;
}" HAVE_FMA_EXTENSIONS)

# Check AVX
SET(CMAKE_REQUIRED_FLAGS)
IF(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    SET(CMAKE_REQUIRED_FLAGS "-mavx")
ELSEIF(MSVC AND NOT CMAKE_CL_64)
    SET(CMAKE_REQUIRED_FLAGS "/arch:AVX")
endif()

CHECK_C_SOURCE_RUNS("
#include <immintrin.h>
int main()
{
    __m256 a = _mm256_set_ps (-1.0f, 2.0f, -3.0f, 4.0f, -1.0f, 2.0f, -3.0f, 4.0f);
    __m256 b = _mm256_set_ps (1.0f, 2.0f, 3.0f, 4.0f, 1.0f, 2.0f, 3.0f, 4.0f);
    __m256 result = _mm256_add_ps (a, b);
    return 0;
}" HAVE_AVX_EXTENSIONS)