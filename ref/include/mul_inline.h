#include "emmintrin.h"
#include "immintrin.h"
#define SWAP_32 0xb1

// Approximately divide each coefficient by t. Carry the results.
inline void reduce_step_narrow_i(
  residue_narrow_t *result, const residue_narrow_t *x) {

  __attribute__((__aligned__(64)))
  static const int32_t carry_permute[16] = {
    0xa, 0x0, 0x1, 0x2, 0x3, 0x4, 0x5, 0x6,
    0x7, 0x8, 0x9, 0xb, 0xc, 0xd, 0xe, 0xf,
  };

  __m512i accum;

  __m512i carry_permute_512 = _mm512_load_si512((__m512i *) carry_permute);
  __m512i mask = _mm512_set1_epi32(0x3ffffff);

  accum = _mm512_load_si512((__m512i*) &x->limbs[0]);

  __m512i error = _mm512_srai_epi32(accum, 26);
  __m512i carry = _mm512_permutexvar_epi32(carry_permute_512, error);
  accum = _mm512_and_si512(accum, mask);
  accum = _mm512_sub_epi32(accum, error);
  __m512i shift_error = _mm512_slli_epi32(error, 4);
  accum = _mm512_add_epi32(accum, shift_error);
  accum = _mm512_add_epi32(accum, carry);
  _mm512_store_si512((__m512i*) &result->limbs[0], accum);
}

// Approximately divide each coefficient by t. Carry the results.
inline void reduce_step_wide_i(
  residue_wide_t *result, const residue_wide_t *x) {

  __attribute__((__aligned__(64)))
  static const int64_t carry_permute_0[8] = {
    0xa, 0x0, 0x1, 0x2, 0x3, 0x4, 0x5, 0x6,
  };

  __attribute__((__aligned__(64)))
  static const int64_t carry_permute_8[8] = {
    0x7, 0x8, 0x9, 0xb, 0xc, 0xd, 0xe, 0xf,
  };

  __m512i accum0;
  __m256i accum8;

  __m512i carry_permute_0_512 = _mm512_load_si512(carry_permute_0);
  __m512i carry_permute_8_512 = _mm512_load_si512(carry_permute_8);
  __m512i mask = _mm512_set1_epi64(0x3ffffff);

  accum0 = _mm512_load_si512((__m512i*) &x->limbs[0]);
//   print8x64(accum0, "accum0");
  __m512i error0 = _mm512_srai_epi64(accum0, 26);
//   print8x64(error0, "error0");
  accum0 = _mm512_and_si512(accum0, mask);
//   print8x64(accum0, "accum0");
  __m512i shift_error0 = _mm512_slli_epi64(error0, 4);
//   print8x64(shift_error0, "shift_error0");

  accum8 = _mm256_load_si256((__m256i*) &x->limbs[8]);
//   print4x64(accum8, "accum8");
  __m256i error8 = _mm256_srai_epi64(accum8, 26);
//   print4x64(error8, "error8");

  // Compute the carries right now that we have the necessary inputs.
  __m512i carry0 = _mm512_permutex2var_epi64(
      error0, carry_permute_0_512, _mm512_castsi256_si512(error8));
//   print8x64(carry0, "carry0");
  __m256i carry8 = _mm512_castsi512_si256(
      _mm512_permutex2var_epi64(
          error0, carry_permute_8_512, _mm512_castsi256_si512(error8)));
//   print4x64(carry8, "carry8");

  accum8 = _mm256_and_si256(accum8, _mm512_castsi512_si256(mask));
//   print4x64(accum8, "accum8");
  __m256i shift_error8 = _mm256_slli_epi64(error8, 4);
//   print4x64(shift_error8, "shift_error8");

  accum0 = _mm512_sub_epi64(accum0, error0);
//   print8x64(accum0, "accum0");
  accum0 = _mm512_add_epi64(accum0, shift_error0);
//   print8x64(accum0, "accum0");
  accum8 = _mm256_sub_epi64(accum8, error8);
//   print4x64(accum8, "accum8");
  accum8 = _mm256_add_epi64(accum8, shift_error8);
//   print4x64(accum8, "accum8");

  accum0 = _mm512_add_epi64(accum0, carry0);
//   print8x64(accum0, "accum0");
  _mm512_store_si512((__m512i*) &result->limbs[0], accum0);
  accum8 = _mm256_add_epi64(accum8, carry8);
//   print4x64(accum8, "accum8");
  _mm256_store_si256((__m256i*) &result->limbs[8], accum8);
}
// Multiply two narrow residues and produce a wide result. The result is reduced
// to 32 bits, but not narrowed for performance reasons.
inline void mul_narrow_i(
  residue_narrow_t *result, const residue_narrow_t *x,
  const residue_narrow_t *y) {

  residue_wide_t temp;

  __m512i lhs_source = _mm512_load_si512((__m512i *) &x->limbs[0]);
  __m512i rhs_source = _mm512_load_si512((__m512i *) &y->limbs[0]);

  // General overview:
  // Take advantage of the symmetry in the cyclical convolution structure in the
  // original Granger Moss Paper
  // Here's the table:
  // Each pair of numbers n,m is shorthand for (x_n - x_m)(y_m - y_n)
  // 10 1    9  2    8  3    7  4    6  5
  // 5  7    4  8    3  9    2 10    1  0
  // 0  2    10 3    9  4    8  5    7  6
  // 6  8    5  9    4 10    3  0    2  1
  // 1  3    0  4    10 5    9  6    8  7
  // 7  9    6 10    5  0    4  1    3  2
  // 2  4    1  5    0  6    10 7    9  8
  // 8 10    7  0    6  1    5  2    4  3
  // --- Invisible break here ---
  // 3  5    2  6    1  7    0  8    10 9
  // 9  0    8  1    7  2    6  3    5  4
  // 4  6    3  7    2  8    1  9    0 10

  // Note that the sequence is the same in the columns. Wherever 10 appears, 4
  // is above it and 5 is below.
  // The invisible break is at the first 8 elements: 512 bits for 64 bit results
  // The column starting with 0 doesn't appear at the top, and the column
  // starting with 4 doesn't appear after the break. It happens that if we
  // alternate top/bottom/top/bottom, starting with the pairing using column 4,
  // we go until we end up with the bottom pairing using column 0:
  // 4  7    1 10    9  2    6  5    3  8
  // 10 2    7  5    4  8    1  0    9  3
  // 5  8    2  0    10 3    7  6    4  9
  // 0  3    8  6    3  9    2  1    10 4      -- Middle 4 rows elided
  //     7  1    10 9    2  6    5  3    8  0
  //     2  7    5  4    8  1    0  9    3  6
  //     8  2    0 10    3  7    6  4    9  1

  // Low slots: 4, 10, 5, 0, 6, 1, 7, 2
  // High slots: 1, 7, 2, 8, 3, 9, 4, 10
  __attribute__((__aligned__(64)))
  static const int32_t permute_4_then_1[16] = {
    4, 1, 10, 7, 5, 2, 0, 8, 6, 3, 1, 9, 7, 4, 2, 10,
  };

  // Low slots: 7, 2, 8, 3, 9, 4, 10, 5
  // High slots: 10, 5, 0, 6, 1, 7, 2, 8
  __attribute__((__aligned__(64)))
  static const int32_t permute_7_then_10[16] = {
    7, 10, 2, 5, 8, 0, 3, 6, 9, 1, 4, 7, 10, 2, 5, 8,
  };

  // Low slots: 9, 4, 10, 5, 0, 6, 1, 7
  // High slots: 6, 1, 7, 2, 8, 3, 9, 4
  __attribute__((__aligned__(64)))
  static const int32_t permute_9_then_6[16] = {
    9, 6, 4, 1, 10, 7, 5, 2, 0, 8, 6, 3, 1, 9, 7, 4,
  };

  // Low slots: 2, 8, 3, 9, 4, 10, 5, 0
  // High slots: 5, 0, 6, 1, 7, 2, 8, 3
  __attribute__((__aligned__(64)))
  static const int32_t permute_2_then_5[16] = {
    2, 5, 8, 0, 3, 6, 9, 1, 4, 7, 10, 2, 5, 8, 0, 3
  };

  // Low slots: 3, 9, 4, 10, 5, 0, 6, 1
  // High slots: 0, 6, 1, 7, 2, 8, 3, 9
  __attribute__((__aligned__(64)))
  static const int32_t permute_3_then_0[16] = {
    3, 0, 9, 6, 4, 1, 10, 7, 5, 2, 0, 8, 6, 3, 1, 9
  };

  // Low slots: 8, 3, 9, 4, 10, 5, 0, 6
  // High slots: (Don't care) 1-15 odds. Passthrough.
  __attribute__((__aligned__(64)))
  static const int32_t permute_8[16] = {
    8, 1, 3, 3, 9, 5, 4, 7, 10, 9, 5, 11, 0, 13, 6, 15,
  };

  // Permutation to collapse the reduction result.
  __attribute__((__aligned__(64)))
  static const int32_t permute_final_result[16] = {
    0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30,
  };

  __m512i _permute_4_then_1 = _mm512_load_si512(permute_4_then_1);
  __m512i _permute_7_then_10 = _mm512_load_si512(permute_7_then_10);
  __m512i _permute_9_then_6 = _mm512_load_si512(permute_9_then_6);
  __m512i _permute_2_then_5 = _mm512_load_si512(permute_2_then_5);
  __m512i _permute_3_then_0 = _mm512_load_si512(permute_3_then_0);
  __m512i _permute_8 = _mm512_load_si512(permute_8);
  __m512i _permute_final_result = _mm512_load_si512(permute_final_result);

  __m512i lhs_4 = _mm512_permutexvar_epi32(_permute_4_then_1, lhs_source);
  __m512i lhs_7 = _mm512_permutexvar_epi32(_permute_7_then_10, lhs_source);
  __m512i rhs_4 = _mm512_permutexvar_epi32(_permute_4_then_1, rhs_source);
  __m512i rhs_7 = _mm512_permutexvar_epi32(_permute_7_then_10, rhs_source);

  __m512i lhs_9 = _mm512_permutexvar_epi32(_permute_9_then_6, lhs_source);
  __m512i rhs_9 = _mm512_permutexvar_epi32(_permute_9_then_6, rhs_source);

  __m512i sublhs_512 = _mm512_sub_epi32(lhs_4, lhs_7);
  __m512i subrhs_512 = _mm512_sub_epi32(rhs_7, rhs_4);

  __m512i accum0 = _mm512_mul_epi32(sublhs_512, subrhs_512);
  // We use a permute and a shift to accomplish the same thing so that they can
  // execute in parallel.
  __m512i lhs_1 = _mm512_srli_epi64(lhs_4, 32);
  __m512i rhs_1 = _mm512_srli_epi64(rhs_4, 32);

  __m256i sublhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(lhs_1),
                                        _mm512_castsi512_si256(lhs_7));
  __m256i subrhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(rhs_7),
                                        _mm512_castsi512_si256(rhs_1));
  __m256i accum8 = _mm256_mul_epi32(sublhs_256, subrhs_256);

  __m512i lhs_10 = _mm512_srli_epi64(lhs_7, 32);
  __m512i rhs_10 = _mm512_srli_epi64(rhs_7, 32);

  // Interleaved here for latency hiding.
  __m512i lhs_2 = _mm512_permutexvar_epi32(_permute_2_then_5, lhs_source);
  sublhs_512 = _mm512_sub_epi32(lhs_10, lhs_1);
  subrhs_512 = _mm512_sub_epi32(rhs_1, rhs_10);

  __m512i rhs_2 = _mm512_permutexvar_epi32(_permute_2_then_5, rhs_source);
  __m512i accum0_temp = _mm512_mul_epi32(sublhs_512, subrhs_512);

  sublhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(lhs_10),
                                _mm512_castsi512_si256(lhs_9));
  subrhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(rhs_9),
                                _mm512_castsi512_si256(rhs_10));
  __m256i accum8_temp = _mm256_mul_epi32(sublhs_256, subrhs_256);

  accum0 = _mm512_add_epi64(accum0, accum0_temp);
  accum8 = _mm256_add_epi64(accum8, accum8_temp);

  sublhs_512 = _mm512_sub_epi32(lhs_9, lhs_2);
  subrhs_512 = _mm512_sub_epi32(rhs_2, rhs_9);
  accum0_temp = _mm512_mul_epi32(sublhs_512, subrhs_512);

  __m512i lhs_6 = _mm512_srli_epi64(lhs_9, 32);
  __m512i rhs_6 = _mm512_srli_epi64(rhs_9, 32);

  __m512i lhs_3 = _mm512_permutexvar_epi32(_permute_3_then_0, lhs_source);
  sublhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(lhs_2),
                                _mm512_castsi512_si256(lhs_6));
  subrhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(rhs_6),
                                _mm512_castsi512_si256(rhs_2));

  __m512i rhs_3 = _mm512_permutexvar_epi32(_permute_3_then_0, rhs_source);
  accum8_temp = _mm256_mul_epi32(sublhs_256, subrhs_256);
  accum0 = _mm512_add_epi64(accum0, accum0_temp);
  accum8 = _mm256_add_epi64(accum8, accum8_temp);

  __m512i lhs_5 = _mm512_srli_epi64(lhs_2, 32);
  __m512i rhs_5 = _mm512_srli_epi64(rhs_2, 32);

  __m512i lhs_8 = _mm512_permutexvar_epi32(_permute_8, lhs_source);
  sublhs_512 = _mm512_sub_epi32(lhs_6, lhs_5);
  subrhs_512 = _mm512_sub_epi32(rhs_5, rhs_6);
  __m512i rhs_8 = _mm512_permutexvar_epi32(_permute_8, rhs_source);
  accum0_temp = _mm512_mul_epi32(sublhs_512, subrhs_512);

  sublhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(lhs_3),
                                _mm512_castsi512_si256(lhs_5));
  subrhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(rhs_5),
                                _mm512_castsi512_si256(rhs_3));
  accum8_temp = _mm256_mul_epi32(sublhs_256, subrhs_256);
  accum0 = _mm512_add_epi64(accum0, accum0_temp);
  accum8 = _mm256_add_epi64(accum8, accum8_temp);

  sublhs_512 = _mm512_sub_epi32(lhs_8, lhs_3);
  subrhs_512 = _mm512_sub_epi32(rhs_3, rhs_8);
  accum0_temp = _mm512_mul_epi32(sublhs_512, subrhs_512);

  __m512i lhs_0 = _mm512_srli_epi64(lhs_3, 32);
  __m512i rhs_0 = _mm512_srli_epi64(rhs_3, 32);

  sublhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(lhs_0),
                                _mm512_castsi512_si256(lhs_8));
  subrhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(rhs_8),
                                _mm512_castsi512_si256(rhs_0));
  accum8_temp = _mm256_mul_epi32(sublhs_256, subrhs_256);
  accum0 = _mm512_add_epi64(accum0, accum0_temp);
  accum8 = _mm256_add_epi64(accum8, accum8_temp);

  _mm512_store_si512((__m512i*) &temp.limbs[0], accum0);
  _mm256_store_si256((__m256i*) &temp.limbs[8], accum8);
  reduce_step_wide_i(&temp, &temp);
  reduce_step_wide_i(&temp, &temp);
  accum0 = _mm512_load_si512((__m512i*) &temp.limbs[0]);
  accum8 = _mm256_load_si256((__m256i*) &temp.limbs[8]);
  __m512i final_result = _mm512_permutex2var_epi32(
    accum0, _permute_final_result, _mm512_castsi256_si512(accum8));
  _mm512_store_si512((__m512i*) &result->limbs[0], final_result);
}

// Square a narrow residue and produce a narrow result. The result is reduced to
// 32 bits.
inline void square_narrow_i(
  residue_narrow_t *result, const residue_narrow_t *x) {

  residue_wide_t temp;

  __m512i lhs_source = _mm512_load_si512((__m512i *) &x->limbs[0]);

  // General overview:
  // Take advantage of the symmetry in the cyclical convolution structure in the
  // original Granger Moss Paper
  // Here's the table:
  // Each pair of numbers n,m is shorthand for (x_n - x_m)(y_m - y_n)
  // 10 1    9  2    8  3    7  4    6  5
  // 5  7    4  8    3  9    2 10    1  0
  // 0  2    10 3    9  4    8  5    7  6
  // 6  8    5  9    4 10    3  0    2  1
  // 1  3    0  4    10 5    9  6    8  7
  // 7  9    6 10    5  0    4  1    3  2
  // 2  4    1  5    0  6    10 7    9  8
  // 8 10    7  0    6  1    5  2    4  3
  // --- Invisible break here ---
  // 3  5    2  6    1  7    0  8    10 9
  // 9  0    8  1    7  2    6  3    5  4
  // 4  6    3  7    2  8    1  9    0 10

  // Note that the sequence is the same in the columns. Wherever 10 appears, 4
  // is above it and 5 is below.
  // The invisible break is at the first 8 elements: 512 bits for 64 bit results
  // The column starting with 0 doesn't appear at the top, and the column
  // starting with 4 doesn't appear after the break. It happens that if we
  // alternate top/bottom/top/bottom, starting with the pairing using column 4,
  // we go until we end up with the bottom pairing using column 0:
  // 4  7    1 10    9  2    6  5    3  8
  // 10 2    7  5    4  8    1  0    9  3
  // 5  8    2  0    10 3    7  6    4  9
  // 0  3    8  6    3  9    2  1    10 4      -- Middle 4 rows elided
  //     7  1    10 9    2  6    5  3    8  0
  //     2  7    5  4    8  1    0  9    3  6
  //     8  2    0 10    3  7    6  4    9  1

  // Low slots: 4, 10, 5, 0, 6, 1, 7, 2
  // High slots: 1, 7, 2, 8, 3, 9, 4, 10
  __attribute__((__aligned__(64)))
  static const int32_t permute_4_then_1[16] = {
    4, 1, 10, 7, 5, 2, 0, 8, 6, 3, 1, 9, 7, 4, 2, 10,
  };

  // Low slots: 7, 2, 8, 3, 9, 4, 10, 5
  // High slots: 10, 5, 0, 6, 1, 7, 2, 8
  __attribute__((__aligned__(64)))
  static const int32_t permute_7_then_10[16] = {
    7, 10, 2, 5, 8, 0, 3, 6, 9, 1, 4, 7, 10, 2, 5, 8,
  };

  // Low slots: 9, 4, 10, 5, 0, 6, 1, 7
  // High slots: 6, 1, 7, 2, 8, 3, 9, 4
  __attribute__((__aligned__(64)))
  static const int32_t permute_9_then_6[16] = {
    9, 6, 4, 1, 10, 7, 5, 2, 0, 8, 6, 3, 1, 9, 7, 4,
  };

  // Low slots: 2, 8, 3, 9, 4, 10, 5, 0
  // High slots: 5, 0, 6, 1, 7, 2, 8, 3
  __attribute__((__aligned__(64)))
  static const int32_t permute_2_then_5[16] = {
    2, 5, 8, 0, 3, 6, 9, 1, 4, 7, 10, 2, 5, 8, 0, 3
  };

  // Low slots: 3, 9, 4, 10, 5, 0, 6, 1
  // High slots: 0, 6, 1, 7, 2, 8, 3, 9
  __attribute__((__aligned__(64)))
  static const int32_t permute_3_then_0[16] = {
    3, 0, 9, 6, 4, 1, 10, 7, 5, 2, 0, 8, 6, 3, 1, 9
  };

  // Low slots: 8, 3, 9, 4, 10, 5, 0, 6
  // High slots: (Don't care) 1-15 odds. Passthrough.
  __attribute__((__aligned__(64)))
  static const int32_t permute_8[16] = {
    8, 1, 3, 3, 9, 5, 4, 7, 10, 9, 5, 11, 0, 13, 6, 15,
  };

  // Permutation to collapse the reduction result.
  __attribute__((__aligned__(64)))
  static const int32_t permute_final_result[16] = {
    0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30,
  };

  __m512i _permute_4_then_1 = _mm512_load_si512(permute_4_then_1);
  __m512i _permute_7_then_10 = _mm512_load_si512(permute_7_then_10);
  __m512i _permute_9_then_6 = _mm512_load_si512(permute_9_then_6);
  __m512i _permute_2_then_5 = _mm512_load_si512(permute_2_then_5);
  __m512i _permute_3_then_0 = _mm512_load_si512(permute_3_then_0);
  __m512i _permute_8 = _mm512_load_si512(permute_8);
  __m512i _permute_final_result = _mm512_load_si512(permute_final_result);

  __m512i lhs_4 = _mm512_permutexvar_epi32(_permute_4_then_1, lhs_source);
  __m512i lhs_7 = _mm512_permutexvar_epi32(_permute_7_then_10, lhs_source);

  __m512i lhs_9 = _mm512_permutexvar_epi32(_permute_9_then_6, lhs_source);
  __m512i lhs_2 = _mm512_permutexvar_epi32(_permute_2_then_5, lhs_source);

  __m512i sublhs_512 = _mm512_sub_epi32(lhs_4, lhs_7);

  __m512i accum0 = _mm512_setzero();
  __m512i accum0_temp = _mm512_mul_epi32(sublhs_512, sublhs_512);
  accum0 = _mm512_sub_epi64(accum0, accum0_temp);
  __m512i lhs_1 = _mm512_srli_epi64(lhs_4, 32);

  __m256i sublhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(lhs_1),
                                        _mm512_castsi512_si256(lhs_7));
  __m256i accum8 = _mm256_setzero_si256();
  __m256i accum8_temp = _mm256_mul_epi32(sublhs_256, sublhs_256);

  accum8 = _mm256_sub_epi64(accum8, accum8_temp);
  __m512i lhs_10 = _mm512_srli_epi64(lhs_7, 32);

  sublhs_512 = _mm512_sub_epi32(lhs_10, lhs_1);

  accum0_temp = _mm512_mul_epi32(sublhs_512, sublhs_512);

  sublhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(lhs_10),
                                _mm512_castsi512_si256(lhs_9));
  accum8_temp = _mm256_mul_epi32(sublhs_256, sublhs_256);

  accum0 = _mm512_sub_epi64(accum0, accum0_temp);
  accum8 = _mm256_sub_epi64(accum8, accum8_temp);

  sublhs_512 = _mm512_sub_epi32(lhs_9, lhs_2);
  accum0_temp = _mm512_mul_epi32(sublhs_512, sublhs_512);

  __m512i lhs_6 = _mm512_srli_epi64(lhs_9, 32);
  __m512i lhs_3 = _mm512_permutexvar_epi32(_permute_3_then_0, lhs_source);
  sublhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(lhs_2),
                                _mm512_castsi512_si256(lhs_6));

  accum8_temp = _mm256_mul_epi32(sublhs_256, sublhs_256);
  accum0 = _mm512_sub_epi64(accum0, accum0_temp);
  accum8 = _mm256_sub_epi64(accum8, accum8_temp);

  __m512i lhs_5 = _mm512_srli_epi64(lhs_2, 32);

  __m512i lhs_8 = _mm512_permutexvar_epi32(_permute_8, lhs_source);
  sublhs_512 = _mm512_sub_epi32(lhs_6, lhs_5);
  accum0_temp = _mm512_mul_epi32(sublhs_512, sublhs_512);

  sublhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(lhs_3),
                                _mm512_castsi512_si256(lhs_5));
  accum8_temp = _mm256_mul_epi32(sublhs_256, sublhs_256);
  accum0 = _mm512_sub_epi64(accum0, accum0_temp);
  accum8 = _mm256_sub_epi64(accum8, accum8_temp);

  sublhs_512 = _mm512_sub_epi32(lhs_8, lhs_3);
  accum0_temp = _mm512_mul_epi32(sublhs_512, sublhs_512);

  __m512i lhs_0 = _mm512_srli_epi64(lhs_3, 32);

  sublhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(lhs_0),
                                _mm512_castsi512_si256(lhs_8));
  accum8_temp = _mm256_mul_epi32(sublhs_256, sublhs_256);
  accum0 = _mm512_sub_epi64(accum0, accum0_temp);
  accum8 = _mm256_sub_epi64(accum8, accum8_temp);

  _mm512_store_si512((__m512i*) &temp.limbs[0], accum0);
  _mm256_store_si256((__m256i*) &temp.limbs[8], accum8);
  reduce_step_wide_i(&temp, &temp);
  reduce_step_wide_i(&temp, &temp);
  accum0 = _mm512_load_si512((__m512i*) &temp.limbs[0]);
  accum8 = _mm256_load_si256((__m256i*) &temp.limbs[8]);
  __m512i final_result = _mm512_permutex2var_epi32(
    accum0, _permute_final_result, _mm512_castsi256_si512(accum8));
  _mm512_store_si512((__m512i*) &result->limbs[0], final_result);
}

// Subtract 2 12x32-bit residues.
inline void sub_narrow_i(
  residue_narrow_t *result, const residue_narrow_t * __restrict x,
  const residue_narrow_t * __restrict y) {

  __m512i lhs = _mm512_load_si512((__m512i*) &x->limbs[0]);
  __m512i rhs = _mm512_load_si512((__m512i*) &y->limbs[0]);
  __m512i sub = _mm512_sub_epi32(lhs, rhs);
  _mm512_store_si512((__m512i*) &result->limbs[0], sub);
}

// negate a 12x32-bit residue.
inline void negate_narrow_i(
  residue_narrow_t *result, const residue_narrow_t *x) {

  __m512i lhs = _mm512_load_si512((__m512i*) &x->limbs[0]);
  __m512i zero = _mm512_setzero();
  __m512i neg = _mm512_sub_epi32(zero, lhs);
  _mm512_store_si512((__m512i*) &result->limbs[0], neg);
}

// Add 2 12x32-bit residues.
inline void add_narrow_i(
  residue_narrow_t *result, const residue_narrow_t * __restrict x,
  const residue_narrow_t * __restrict y) {

  __m512i lhs = _mm512_load_si512((__m512i*) &x->limbs[0]);
  __m512i rhs = _mm512_load_si512((__m512i*) &y->limbs[0]);
  __m512i add = _mm512_add_epi32(lhs, rhs);
  _mm512_store_si512((__m512i*) &result->limbs[0], add);
}

// Scale a narrow residue by 2.
inline void double_narrow_i(
  residue_narrow_t *result, const residue_narrow_t *x) {

  __m512i lhs = _mm512_load_si512((__m512i*) &x->limbs[0]);
  __m512i dub = _mm512_slli_epi32(lhs, 1);
  _mm512_store_si512((__m512i*) &result->limbs[0], dub);
}

