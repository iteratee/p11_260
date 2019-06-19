#include <stdint.h>
#include "f11_260.h"
#include "emmintrin.h"
#include "immintrin.h"

residue_wide_t zero_wide = {0};
residue_wide_t one_wide = {
  .limbs = {0, 1},
};
residue_narrow_t zero_narrow = {0};
residue_narrow_t one_narrow = {
  .limbs = {0, 1},
};

#define NVECTORS 3
#define VECTWIDTH 4

__attribute__((__aligned__(32)))
static const int32_t COLLAPSE[8] = { 0, 2, 4, 6, 4, 5, 6, 7 };

// Shrink to 32 bits. Assumes reduction has already occurred, and wide storage
// is being used for vector compatibility.
void narrow(residue_narrow_t *result, const residue_wide_t * __restrict w) {
  __m256i collapse_perm = _mm256_load_si256((__m256i*) COLLAPSE);
  __m128i packed_result;
  #pragma clang loop unroll(full)
  for (int i = 0; i < NVECTORS; ++i) {
    __m256i x = _mm256_load_si256((__m256i*) (&w->limbs[i * VECTWIDTH]));
    packed_result = _mm256_castsi256_si128(
        _mm256_permutevar8x32_epi32(x, collapse_perm));
    _mm_store_si128((__m128i*) &result->limbs[i * VECTWIDTH], packed_result);
  }
}

// Reduce to 10 limbs. For final compression.
void narrow_reduce(
  residue_narrow_reduced_t *result, const residue_narrow_t * __restrict w) {
  residue_narrow_t temp;

  __m128i x = _mm_load_si128((__m128i *) (&w->limbs[0 * VECTWIDTH]));
  __m128i x10 = _mm_broadcastd_epi32(x);
  x = _mm_sub_epi32(x, x10);
  _mm_store_si128((__m128i *) &temp.limbs[0 * VECTWIDTH], x);
  x = _mm_load_si128((__m128i *) (&w->limbs[1 * VECTWIDTH]));
  x = _mm_sub_epi32(x, x10);
  _mm_store_si128((__m128i *) &temp.limbs[1 * VECTWIDTH], x);
  x = _mm_load_si128((__m128i *) (&w->limbs[2 * VECTWIDTH]));
  x = _mm_sub_epi32(x, x10);
  _mm_store_si128((__m128i *) &temp.limbs[2 * VECTWIDTH], x);

  reduce_step_narrow(&temp, &temp);

  // May want to use vpalignr here.
  #pragma clang loop unroll(full)
  for (int i = 0; i < NLIMBS_REDUCED; ++i) {
    result->limbs[i] = temp.limbs[i+1] - temp.limbs[0];
  }
}

// Reduce to unique representative.
// This is expensive. Only used for final signature or DH Key
void narrow_complete(
  residue_narrow_reduced_t *result, const residue_narrow_t * __restrict w) {

  residue_narrow_t temp;
  for (int i = 0; i < NLIMBS; ++i) {
    temp.limbs[i] = w->limbs[i] - w->limbs[0];
  }

  // This may be combined with the final reduction from a multiply.
  reduce_step_narrow(&temp, &temp);

  int gt_mask = 0;
  int lt_mask = 0;
  int32_t limit[NLIMBS];
  for (int i = 1; i < NLIMBS; ++i) {
    temp.limbs[i] = temp.limbs[i] - temp.limbs[0];
    temp.limbs[i] += 1 & gt_mask;
    temp.limbs[i] -= 1 & gt_mask;
    gt_mask = -(temp.limbs[i] > T);
    lt_mask = -(temp.limbs[i] < 0);
    temp.limbs[i] -= (T & gt_mask);
    temp.limbs[i] += (T & lt_mask);
  }
  for (int i = 1; i < NLIMBS - 1; ++i) {
    temp.limbs[i] -= temp.limbs[11];
    limit[i] = T;
  }
  int64_t all_t = -1;
  for (int i = NLIMBS - 3; i >= 1; --i) {
    all_t &= -(temp.limbs[i+1] == T);
    limit[i] -= 1 & (~all_t);
  }
  gt_mask = 0;
  lt_mask = 0;
  for (int i = 1; i < NLIMBS - 1; ++i) {
    temp.limbs[i] += 1 & gt_mask;
    temp.limbs[i] -= 1 & lt_mask;
    gt_mask = -(temp.limbs[i] > limit[i]);
    lt_mask = -(temp.limbs[i] < 0);
    temp.limbs[i] -= (T & gt_mask);
    temp.limbs[i] += (T & lt_mask);
    result->limbs[i-1] = temp.limbs[i];
  }
}

// Reduce to mostly unique representative.
// All coefficients are reduced to 0 <= xi <= t
// Unique up to carries (xi == t) => (xi = 0; x[i+1] += 1);
// This is sufficient to determine if x is even or odd.
// Still pretty expensive. Used in point compression.
void narrow_partial_complete(
  residue_narrow_reduced_t *result, const residue_narrow_t * __restrict w) {

  residue_narrow_t temp;
  for (int i = 0; i < NLIMBS; ++i) {
    temp.limbs[i] = w->limbs[i] - w->limbs[0];
  }

  // This may be combined with the final reduction from a multiply.
  reduce_step_narrow(&temp, &temp);

  int gt_mask = 0;
  int lt_mask = 0;
  for (int i = 1; i < NLIMBS; ++i) {
    temp.limbs[i] = temp.limbs[i] - temp.limbs[0];
    temp.limbs[i] += 1 & gt_mask;
    temp.limbs[i] -= 1 & gt_mask;
    gt_mask = -(temp.limbs[i] > T);
    lt_mask = -(temp.limbs[i] < 0);
    temp.limbs[i] -= (T & gt_mask);
    temp.limbs[i] += (T & lt_mask);
  }
  for (int i = 1; i < NLIMBS; ++i) {
    temp.limbs[i] -= temp.limbs[11];
  }
  gt_mask = 0;
  lt_mask = 0;
  for (int i = 1; i < NLIMBS - 1; ++i) {
    temp.limbs[i] += 1 & gt_mask;
    temp.limbs[i] -= 1 & lt_mask;
    gt_mask = -(temp.limbs[i] > T);
    lt_mask = -(temp.limbs[i] < 0);
    temp.limbs[i] -= (T & gt_mask);
    temp.limbs[i] += (T & lt_mask);
    result->limbs[i-1] = temp.limbs[i];
  }
}

int is_odd(residue_narrow_reduced_t *x) {
  int result = 0;
  for (int i = 0; i < NLIMBS_REDUCED; ++i) {
    result ^= x->limbs[i] & 0x1;
  }
  return result;
}

// Copy a 12x64-bit residue
void copy_wide(
  residue_wide_t *result, const residue_wide_t * __restrict x) {

  for (int i = 0; i < NLIMBS; ++i) {
    result->limbs[i] = x->limbs[i];
  }
}

// Copy a 12x32-bit residue
void copy_narrow(
  residue_narrow_t *result, const residue_narrow_t * __restrict x) {

  for (int i = 0; i < NLIMBS; ++i) {
    result->limbs[i] = x->limbs[i];
  }
}

// Copy a 10x32-bit residue
void copy_narrow_reduced(
  residue_narrow_reduced_t *result,
  const residue_narrow_reduced_t * __restrict x) {

  for (int i = 0; i < NLIMBS_REDUCED; ++i) {
    result->limbs[i] = x->limbs[i];
  }
}

// Produce a 64-bit residue
void widen(
  residue_wide_t *result, const residue_narrow_t * __restrict x) {
  for (int i = 0; i < NLIMBS; ++i) {
    result->limbs[i] = x->limbs[i];
  }
}

// Subtract 2 12x64-bit residues.
void sub_wide(
  residue_wide_t *result, const residue_wide_t * __restrict x,
  const residue_wide_t * __restrict y) {

  for (int i = 0; i < NLIMBS; ++i) {
    result->limbs[i] = x->limbs[i] - y->limbs[i];
  }
}

// negate a 12x64-bit residue.
void negate_wide(residue_wide_t *result, const residue_wide_t *x) {

  __m256i zero = _mm256_setzero_si256();
  #pragma clang loop unroll(full)
  for (int i = 0; i < NVECTORS; ++i) {
    __m256i xv = _mm256_load_si256((__m256i*) (&x->limbs[i * VECTWIDTH]));
    xv = _mm256_sub_epi64(zero, xv);
    _mm256_store_si256((__m256i*) &result->limbs[i * VECTWIDTH], xv);
  }
}

// negate a 12x32-bit residue.
void negate_narrow(
  residue_narrow_t *result, const residue_narrow_t *x) {

  for (int i = 0; i < NLIMBS; ++i) {
    result->limbs[i] = -(x->limbs[i]);
  }
}

// Add 2 12x32-bit residues.
void add_narrow(
  residue_narrow_t *result, const residue_narrow_t * __restrict x,
  const residue_narrow_t * __restrict y) {

  for (int i = 0; i < NLIMBS; ++i) {
    result->limbs[i] = x->limbs[i] + y->limbs[i];
  }
}

// Add 2 12x64-bit residues.
void add_wide(
  residue_wide_t *result, const residue_wide_t * __restrict x,
  const residue_wide_t * __restrict y) {

  for (int i = 0; i < NLIMBS; ++i) {
    result->limbs[i] = x->limbs[i] + y->limbs[i];
  }
}

// Scale a wide residue by 2.
void double_wide(
  residue_wide_t *result, const residue_wide_t *x) {

  for (int i = 0; i < NLIMBS; ++i) {
    result->limbs[i] = x->limbs[i] << 1;
  }
}

#define wrap(x) (((x + (NLIMBS - 1)) % (NLIMBS - 1)) + 1)
// Multiply two wide residues, and produce a wide result. The result is reduced
// to 32 bits, but not narrowed for performance reasons.
void mul_wide(
  residue_wide_t *result, const residue_wide_t *x, const residue_wide_t *y) {

  residue_wide_t temp;

  __m256i sublhs, subrhs, mul; // Temporaries for the actual sub sub mul
  __m128i sublhs128_low7, subrhs128_low7;
  __m128i sublhs128_low3, subrhs128_low3;
  __m128i sublhs128_low10, subrhs128_low10;
  __m128i sublhs128_high, subrhs128_high; // Temporaries for single lane blends
  __m256i accum7, accum3, accum10; // Three accumulators

  // We do 7-10 first.
  __m256i lhs_10 = _mm256_set1_epi32(x->limbs[0]);
  __m256i rhs_10 = _mm256_set1_epi32(y->limbs[0]);

  __m256i lhs_7_10 = _mm256_load_si256((__m256i*) &x->limbs[8]);
  __m256i rhs_7_10 = _mm256_load_si256((__m256i*) &y->limbs[8]);

  __m256i lhs_0 = _mm256_set1_epi32(x->limbs[1]);
  __m256i rhs_0 = _mm256_set1_epi32(y->limbs[1]);

  // 0 * [7 8 9 10]
  sublhs = _mm256_sub_epi32(lhs_0, lhs_7_10);
  subrhs = _mm256_sub_epi32(rhs_7_10, rhs_0);
  accum7 = _mm256_mul_epi32(sublhs, subrhs);

  __m256i lhs_6_9 = _mm256_loadu_si256((__m256i*) &x->limbs[7]);
  __m256i rhs_6_9 = _mm256_loadu_si256((__m256i*) &y->limbs[7]);

  // 10 * [8 9]
  __m128i lhs_8_9 = _mm256_extracti128_si256(lhs_6_9, 1);
  __m128i rhs_8_9 = _mm256_extracti128_si256(rhs_6_9, 1);

  sublhs128_low7 = _mm_sub_epi32(_mm256_castsi256_si128(lhs_10), lhs_8_9);
  subrhs128_low7 = _mm_sub_epi32(rhs_8_9, _mm256_castsi256_si128(rhs_10));

  __m256i lhs_1 = _mm256_set1_epi32(x->limbs[2]);
  __m256i rhs_1 = _mm256_set1_epi32(y->limbs[2]);

  // Second round is interleaved starting here
  __m256i lhs_8 = _mm256_set1_epi32(x->limbs[9]);
  __m256i rhs_8 = _mm256_set1_epi32(y->limbs[9]);

  // 1 * [6 7 8 9]
  sublhs = _mm256_sub_epi32(lhs_1, lhs_6_9);
  subrhs = _mm256_sub_epi32(rhs_6_9, rhs_1);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum7 = _mm256_add_epi64(accum7, mul);

  // 8 * [6 7]
  sublhs128_low3 = _mm_sub_epi32(_mm256_castsi256_si128(lhs_8),
                                 _mm256_castsi256_si128(lhs_6_9));
  subrhs128_low3 = _mm_sub_epi32(_mm256_castsi256_si128(rhs_6_9),
                                 _mm256_castsi256_si128(rhs_8));

  __m256i lhs_5_8 = _mm256_loadu_si256((__m256i*) &x->limbs[6]);
  __m256i rhs_5_8 = _mm256_loadu_si256((__m256i*) &y->limbs[6]);

  __m256i lhs_4 = _mm256_set1_epi32(x->limbs[5]);
  __m256i rhs_4 = _mm256_set1_epi32(y->limbs[5]);

  __m256i lhs_2 = _mm256_set1_epi32(x->limbs[3]);
  __m256i rhs_2 = _mm256_set1_epi32(y->limbs[3]);

  __m256i lhs_9 = _mm256_set1_epi32(x->limbs[10]);
  __m256i rhs_9 = _mm256_set1_epi32(y->limbs[10]);

  // 2 * [5 6 7 8]
  sublhs = _mm256_sub_epi32(lhs_2, lhs_5_8);
  subrhs = _mm256_sub_epi32(rhs_5_8, rhs_2);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum7 = _mm256_add_epi64(accum7, mul);

  // 9 * [5 6 7 8]
  sublhs = _mm256_sub_epi32(lhs_9, lhs_5_8);
  subrhs = _mm256_sub_epi32(rhs_5_8, rhs_9);
  accum3 = _mm256_mul_epi32(sublhs, subrhs);

  __m256i lhs_4_7 = _mm256_loadu_si256((__m256i*) &x->limbs[5]);
  __m256i rhs_4_7 = _mm256_loadu_si256((__m256i*) &y->limbs[5]);

  __m256i lhs_3 = _mm256_set1_epi32(x->limbs[4]);
  __m256i rhs_3 = _mm256_set1_epi32(y->limbs[4]);

  // Third round is interleaved starting here
  __m256i lhs_6 = _mm256_set1_epi32(x->limbs[7]);
  __m256i rhs_6 = _mm256_set1_epi32(y->limbs[7]);

  // 3 * [4 5 6 7]
  sublhs = _mm256_sub_epi32(lhs_3, lhs_4_7);
  subrhs = _mm256_sub_epi32(rhs_4_7, rhs_3);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum7 = _mm256_add_epi64(accum7, mul);

  // 10 * [4 5 6 7]
  sublhs = _mm256_sub_epi32(lhs_10, lhs_4_7);
  subrhs = _mm256_sub_epi32(rhs_4_7, rhs_10);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum3 = _mm256_add_epi64(accum3, mul);

  // 6 * [4 5]
  sublhs128_low10 = _mm_sub_epi32(_mm256_castsi256_si128(lhs_6),
                                  _mm256_castsi256_si128(lhs_4_7));
  subrhs128_low10 = _mm_sub_epi32(_mm256_castsi256_si128(rhs_4_7),
                                  _mm256_castsi256_si128(rhs_6));

  __m256i lhs_3_6 = _mm256_load_si256((__m256i*) &x->limbs[4]);
  __m256i rhs_3_6 = _mm256_load_si256((__m256i*) &y->limbs[4]);

  __m256i lhs_7 = _mm256_set1_epi32(x->limbs[8]);
  __m256i rhs_7 = _mm256_set1_epi32(y->limbs[8]);

  // 0 * [3 4 5 6]
  sublhs = _mm256_sub_epi32(lhs_0, lhs_3_6);
  subrhs = _mm256_sub_epi32(rhs_3_6, rhs_0);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum3 = _mm256_add_epi64(accum3, mul);

  // 7 * [3 4 5 6]
  sublhs = _mm256_sub_epi32(lhs_7, lhs_3_6);
  subrhs = _mm256_sub_epi32(rhs_3_6, rhs_7);
  accum10 = _mm256_mul_epi32(sublhs, subrhs);

  // 4 * [5 6]
  sublhs = _mm256_sub_epi32(lhs_4, lhs_3_6);
  subrhs = _mm256_sub_epi32(rhs_3_6, rhs_4);
  sublhs = _mm256_blend_epi32(
      sublhs, _mm256_castsi128_si256(sublhs128_low7), 0x0f);
  subrhs = _mm256_blend_epi32(
      subrhs, _mm256_castsi128_si256(subrhs128_low7), 0x0f);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum7 = _mm256_add_epi64(accum7, mul);

  // Round 1 drops off here.
  _mm256_store_si256((__m256i*) &temp.limbs[8], accum7);

  __m256i lhs_2_5 = _mm256_loadu_si256((__m256i*) &x->limbs[3]);
  __m256i rhs_2_5 = _mm256_loadu_si256((__m256i*) &y->limbs[3]);

  // 1 * [2 3 4 5]
  sublhs = _mm256_sub_epi32(lhs_1, lhs_2_5);
  subrhs = _mm256_sub_epi32(rhs_2_5, rhs_1);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum3 = _mm256_add_epi64(accum3, mul);

  // 8 * [2 3 4 5]
  sublhs = _mm256_sub_epi32(lhs_8, lhs_2_5);
  subrhs = _mm256_sub_epi32(rhs_2_5, rhs_8);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum10 = _mm256_add_epi64(accum10, mul);

  __m256i lhs_1_4 = _mm256_loadu_si256((__m256i*) &x->limbs[2]);
  __m256i rhs_1_4 = _mm256_loadu_si256((__m256i*) &y->limbs[2]);

  // 9 * [1 2 3 4]
  sublhs = _mm256_sub_epi32(lhs_9, lhs_1_4);
  subrhs = _mm256_sub_epi32(rhs_1_4, rhs_9);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum10 = _mm256_add_epi64(accum10, mul);

  // 0 * [1 2]
  // We do this here, only because we never load [10 0 1 2]
  sublhs128_high = _mm_sub_epi32(_mm256_castsi256_si128(lhs_0),
                                 _mm256_castsi256_si128(lhs_1_4));
  subrhs128_high = _mm_sub_epi32(_mm256_castsi256_si128(rhs_1_4),
                                 _mm256_castsi256_si128(rhs_0));
  sublhs = _mm256_inserti128_si256(
    _mm256_castsi128_si256(sublhs128_low10),
    sublhs128_high, 1);
  subrhs = _mm256_inserti128_si256(
    _mm256_castsi128_si256(subrhs128_low10),
    subrhs128_high, 1);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum10 = _mm256_add_epi64(accum10, mul);

  // 2 * [3 4]
  sublhs = _mm256_sub_epi32(lhs_2, lhs_1_4);
  subrhs = _mm256_sub_epi32(rhs_1_4, rhs_2);
  sublhs = _mm256_blend_epi32(
      sublhs, _mm256_castsi128_si256(sublhs128_low3), 0x0f);
  subrhs = _mm256_blend_epi32(
      subrhs, _mm256_castsi128_si256(subrhs128_low3), 0x0f);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum3 = _mm256_add_epi64(accum3, mul);
  // Round 2 drops off here.
  _mm256_store_si256((__m256i*) &temp.limbs[4], accum3);

  __m256i lhs_0_3 = _mm256_loadu_si256((__m256i*) &x->limbs[1]);
  __m256i rhs_0_3 = _mm256_loadu_si256((__m256i*) &y->limbs[1]);

  // 10 * [0 1 2 3]
  sublhs = _mm256_sub_epi32(lhs_10, lhs_0_3);
  subrhs = _mm256_sub_epi32(rhs_0_3, rhs_10);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum10 = _mm256_add_epi64(accum10, mul);
  // Round 3 drops off here.

  _mm256_store_si256((__m256i*) &temp.limbs[0], accum10);

  reduce_step_wide(&temp, &temp);
  reduce_step_wide(result, &temp);
}

static inline __m256i load_extend_32_64(__m128i *x) {
  return _mm256_cvtepi32_epi64(_mm_load_si128(x));
}

static inline __m256i loadu_extend_32_64(__m128i *x) {
  return _mm256_cvtepi32_epi64(_mm_loadu_si128(x));
}
// Multiply a wide residues by a narrow and produce a wide result. The result is
// reduced to 32 bits, but not narrowed for performance reasons.
void mul_wide_narrow(
  residue_wide_t *result, const residue_wide_t *x, const residue_narrow_t *y) {

  residue_wide_t temp;

  __m256i sublhs, subrhs, mul; // Temporaries for the actual sub sub mul
  __m128i sublhs128_low7, subrhs128_low7;
  __m128i sublhs128_low3, subrhs128_low3;
  __m128i sublhs128_low10, subrhs128_low10;
  __m128i sublhs128_high, subrhs128_high; // Temporaries for single lane blends
  __m256i accum7, accum3, accum10; // Three accumulators

  // We do 7-10 first.
  __m256i lhs_10 = _mm256_set1_epi32(x->limbs[0]);
  __m256i rhs_10 = _mm256_set1_epi32(y->limbs[0]);

  __m256i lhs_7_10 = _mm256_load_si256((__m256i*) &x->limbs[8]);
  __m256i rhs_7_10 = load_extend_32_64((__m128i*) &y->limbs[8]);

  __m256i lhs_0 = _mm256_set1_epi32(x->limbs[1]);
  __m256i rhs_0 = _mm256_set1_epi32(y->limbs[1]);

  // 0 * [7 8 9 10]
  sublhs = _mm256_sub_epi32(lhs_0, lhs_7_10);
  subrhs = _mm256_sub_epi32(rhs_7_10, rhs_0);
  accum7 = _mm256_mul_epi32(sublhs, subrhs);

  __m256i lhs_6_9 = _mm256_loadu_si256((__m256i*) &x->limbs[7]);
  __m256i rhs_6_9 = loadu_extend_32_64((__m128i*) &y->limbs[7]);

  // 10 * [8 9]
  __m128i lhs_8_9 = _mm256_extracti128_si256(lhs_6_9, 1);
  __m128i rhs_8_9 = _mm256_extracti128_si256(rhs_6_9, 1);

  sublhs128_low7 = _mm_sub_epi32(_mm256_castsi256_si128(lhs_10), lhs_8_9);
  subrhs128_low7 = _mm_sub_epi32(rhs_8_9, _mm256_castsi256_si128(rhs_10));

  __m256i lhs_1 = _mm256_set1_epi32(x->limbs[2]);
  __m256i rhs_1 = _mm256_set1_epi32(y->limbs[2]);

  // Second round is interleaved starting here
  __m256i lhs_8 = _mm256_set1_epi32(x->limbs[9]);
  __m256i rhs_8 = _mm256_set1_epi32(y->limbs[9]);

  // 1 * [6 7 8 9]
  sublhs = _mm256_sub_epi32(lhs_1, lhs_6_9);
  subrhs = _mm256_sub_epi32(rhs_6_9, rhs_1);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum7 = _mm256_add_epi64(accum7, mul);

  // 8 * [6 7]
  sublhs128_low3 = _mm_sub_epi32(_mm256_castsi256_si128(lhs_8),
                                 _mm256_castsi256_si128(lhs_6_9));
  subrhs128_low3 = _mm_sub_epi32(_mm256_castsi256_si128(rhs_6_9),
                                 _mm256_castsi256_si128(rhs_8));

  __m256i lhs_5_8 = _mm256_loadu_si256((__m256i*) &x->limbs[6]);
  __m256i rhs_5_8 = loadu_extend_32_64((__m128i*) &y->limbs[6]);

  __m256i lhs_4 = _mm256_set1_epi32(x->limbs[5]);
  __m256i rhs_4 = _mm256_set1_epi32(y->limbs[5]);

  __m256i lhs_2 = _mm256_set1_epi32(x->limbs[3]);
  __m256i rhs_2 = _mm256_set1_epi32(y->limbs[3]);

  __m256i lhs_9 = _mm256_set1_epi32(x->limbs[10]);
  __m256i rhs_9 = _mm256_set1_epi32(y->limbs[10]);

  // 2 * [5 6 7 8]
  sublhs = _mm256_sub_epi32(lhs_2, lhs_5_8);
  subrhs = _mm256_sub_epi32(rhs_5_8, rhs_2);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum7 = _mm256_add_epi64(accum7, mul);

  // 9 * [5 6 7 8]
  sublhs = _mm256_sub_epi32(lhs_9, lhs_5_8);
  subrhs = _mm256_sub_epi32(rhs_5_8, rhs_9);
  accum3 = _mm256_mul_epi32(sublhs, subrhs);

  __m256i lhs_4_7 = _mm256_loadu_si256((__m256i*) &x->limbs[5]);
  __m256i rhs_4_7 = loadu_extend_32_64((__m128i*) &y->limbs[5]);

  __m256i lhs_3 = _mm256_set1_epi32(x->limbs[4]);
  __m256i rhs_3 = _mm256_set1_epi32(y->limbs[4]);

  // Third round is interleaved starting here
  __m256i lhs_6 = _mm256_set1_epi32(x->limbs[7]);
  __m256i rhs_6 = _mm256_set1_epi32(y->limbs[7]);

  // 3 * [4 5 6 7]
  sublhs = _mm256_sub_epi32(lhs_3, lhs_4_7);
  subrhs = _mm256_sub_epi32(rhs_4_7, rhs_3);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum7 = _mm256_add_epi64(accum7, mul);

  // 10 * [4 5 6 7]
  sublhs = _mm256_sub_epi32(lhs_10, lhs_4_7);
  subrhs = _mm256_sub_epi32(rhs_4_7, rhs_10);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum3 = _mm256_add_epi64(accum3, mul);

  // 6 * [4 5]
  sublhs128_low10 = _mm_sub_epi32(_mm256_castsi256_si128(lhs_6),
                                  _mm256_castsi256_si128(lhs_4_7));
  subrhs128_low10 = _mm_sub_epi32(_mm256_castsi256_si128(rhs_4_7),
                                  _mm256_castsi256_si128(rhs_6));

  __m256i lhs_3_6 = _mm256_load_si256((__m256i*) &x->limbs[4]);
  __m256i rhs_3_6 = load_extend_32_64((__m128i*) &y->limbs[4]);

  __m256i lhs_7 = _mm256_set1_epi32(x->limbs[8]);
  __m256i rhs_7 = _mm256_set1_epi32(y->limbs[8]);

  // 0 * [3 4 5 6]
  sublhs = _mm256_sub_epi32(lhs_0, lhs_3_6);
  subrhs = _mm256_sub_epi32(rhs_3_6, rhs_0);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum3 = _mm256_add_epi64(accum3, mul);

  // 7 * [3 4 5 6]
  sublhs = _mm256_sub_epi32(lhs_7, lhs_3_6);
  subrhs = _mm256_sub_epi32(rhs_3_6, rhs_7);
  accum10 = _mm256_mul_epi32(sublhs, subrhs);

  // 4 * [5 6]
  sublhs = _mm256_sub_epi32(lhs_4, lhs_3_6);
  subrhs = _mm256_sub_epi32(rhs_3_6, rhs_4);
  sublhs = _mm256_blend_epi32(
      sublhs, _mm256_castsi128_si256(sublhs128_low7), 0x0f);
  subrhs = _mm256_blend_epi32(
      subrhs, _mm256_castsi128_si256(subrhs128_low7), 0x0f);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum7 = _mm256_add_epi64(accum7, mul);

  // Round 1 drops off here.
  _mm256_store_si256((__m256i*) &temp.limbs[8], accum7);

  __m256i lhs_2_5 = _mm256_loadu_si256((__m256i*) &x->limbs[3]);
  __m256i rhs_2_5 = loadu_extend_32_64((__m128i*) &y->limbs[3]);

  // 1 * [2 3 4 5]
  sublhs = _mm256_sub_epi32(lhs_1, lhs_2_5);
  subrhs = _mm256_sub_epi32(rhs_2_5, rhs_1);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum3 = _mm256_add_epi64(accum3, mul);

  // 8 * [2 3 4 5]
  sublhs = _mm256_sub_epi32(lhs_8, lhs_2_5);
  subrhs = _mm256_sub_epi32(rhs_2_5, rhs_8);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum10 = _mm256_add_epi64(accum10, mul);

  __m256i lhs_1_4 = _mm256_loadu_si256((__m256i*) &x->limbs[2]);
  __m256i rhs_1_4 = loadu_extend_32_64((__m128i*) &y->limbs[2]);

  // 9 * [1 2 3 4]
  sublhs = _mm256_sub_epi32(lhs_9, lhs_1_4);
  subrhs = _mm256_sub_epi32(rhs_1_4, rhs_9);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum10 = _mm256_add_epi64(accum10, mul);

  // 0 * [1 2]
  // We do this here, only because we never load [10 0 1 2]
  sublhs128_high = _mm_sub_epi32(_mm256_castsi256_si128(lhs_0),
                                 _mm256_castsi256_si128(lhs_1_4));
  subrhs128_high = _mm_sub_epi32(_mm256_castsi256_si128(rhs_1_4),
                                 _mm256_castsi256_si128(rhs_0));
  sublhs = _mm256_inserti128_si256(
    _mm256_castsi128_si256(sublhs128_low10),
    sublhs128_high, 1);
  subrhs = _mm256_inserti128_si256(
    _mm256_castsi128_si256(subrhs128_low10),
    subrhs128_high, 1);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum10 = _mm256_add_epi64(accum10, mul);

  // 2 * [3 4]
  sublhs = _mm256_sub_epi32(lhs_2, lhs_1_4);
  subrhs = _mm256_sub_epi32(rhs_1_4, rhs_2);
  sublhs = _mm256_blend_epi32(
      sublhs, _mm256_castsi128_si256(sublhs128_low3), 0x0f);
  subrhs = _mm256_blend_epi32(
      subrhs, _mm256_castsi128_si256(subrhs128_low3), 0x0f);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum3 = _mm256_add_epi64(accum3, mul);
  // Round 2 drops off here.
  _mm256_store_si256((__m256i*) &temp.limbs[4], accum3);

  __m256i lhs_0_3 = _mm256_loadu_si256((__m256i*) &x->limbs[1]);
  __m256i rhs_0_3 = loadu_extend_32_64((__m128i*) &y->limbs[1]);

  // 10 * [0 1 2 3]
  sublhs = _mm256_sub_epi32(lhs_10, lhs_0_3);
  subrhs = _mm256_sub_epi32(rhs_0_3, rhs_10);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum10 = _mm256_add_epi64(accum10, mul);
  // Round 3 drops off here.

  _mm256_store_si256((__m256i*) &temp.limbs[0], accum10);
  reduce_step_wide(&temp, &temp);
  reduce_step_wide(result, &temp);
}

// Multiply two narrow residues and produce a wide result. The result is reduced
// to 32 bits, but not narrowed for performance reasons.
void mul_narrow(
  residue_wide_t *result, const residue_narrow_t *x,
  const residue_narrow_t *y) {

  residue_wide_t temp;

  __m256i sublhs, subrhs, mul; // Temporaries for the actual sub sub mul
  __m128i sublhs128_low7, subrhs128_low7;
  __m128i sublhs128_low3, subrhs128_low3;
  __m128i sublhs128_low10, subrhs128_low10;
  __m128i sublhs128_high, subrhs128_high; // Temporaries for single lane blends
  __m256i accum7, accum3, accum10; // Three accumulators

  // We do 7-10 first.
  __m256i lhs_10 = _mm256_set1_epi32(x->limbs[0]);
  __m256i rhs_10 = _mm256_set1_epi32(y->limbs[0]);

  __m256i lhs_7_10 = load_extend_32_64((__m128i*) &x->limbs[8]);
  __m256i rhs_7_10 = load_extend_32_64((__m128i*) &y->limbs[8]);

  __m256i lhs_0 = _mm256_set1_epi32(x->limbs[1]);
  __m256i rhs_0 = _mm256_set1_epi32(y->limbs[1]);

  // 0 * [7 8 9 10]
  sublhs = _mm256_sub_epi32(lhs_0, lhs_7_10);
  subrhs = _mm256_sub_epi32(rhs_7_10, rhs_0);
  accum7 = _mm256_mul_epi32(sublhs, subrhs);

  __m256i lhs_6_9 = loadu_extend_32_64((__m128i*) &x->limbs[7]);
  __m256i rhs_6_9 = loadu_extend_32_64((__m128i*) &y->limbs[7]);

  // 10 * [8 9]
  __m128i lhs_8_9 = _mm256_extracti128_si256(lhs_6_9, 1);
  __m128i rhs_8_9 = _mm256_extracti128_si256(rhs_6_9, 1);

  sublhs128_low7 = _mm_sub_epi32(_mm256_castsi256_si128(lhs_10), lhs_8_9);
  subrhs128_low7 = _mm_sub_epi32(rhs_8_9, _mm256_castsi256_si128(rhs_10));

  __m256i lhs_1 = _mm256_set1_epi32(x->limbs[2]);
  __m256i rhs_1 = _mm256_set1_epi32(y->limbs[2]);

  // Second round is interleaved starting here
  __m256i lhs_8 = _mm256_set1_epi32(x->limbs[9]);
  __m256i rhs_8 = _mm256_set1_epi32(y->limbs[9]);

  // 1 * [6 7 8 9]
  sublhs = _mm256_sub_epi32(lhs_1, lhs_6_9);
  subrhs = _mm256_sub_epi32(rhs_6_9, rhs_1);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum7 = _mm256_add_epi64(accum7, mul);

  // 8 * [6 7]
  sublhs128_low3 = _mm_sub_epi32(_mm256_castsi256_si128(lhs_8),
                                 _mm256_castsi256_si128(lhs_6_9));
  subrhs128_low3 = _mm_sub_epi32(_mm256_castsi256_si128(rhs_6_9),
                                 _mm256_castsi256_si128(rhs_8));

  __m256i lhs_5_8 = loadu_extend_32_64((__m128i*) &x->limbs[6]);
  __m256i rhs_5_8 = loadu_extend_32_64((__m128i*) &y->limbs[6]);

  __m256i lhs_4 = _mm256_set1_epi32(x->limbs[5]);
  __m256i rhs_4 = _mm256_set1_epi32(y->limbs[5]);

  __m256i lhs_2 = _mm256_set1_epi32(x->limbs[3]);
  __m256i rhs_2 = _mm256_set1_epi32(y->limbs[3]);

  __m256i lhs_9 = _mm256_set1_epi32(x->limbs[10]);
  __m256i rhs_9 = _mm256_set1_epi32(y->limbs[10]);

  // 2 * [5 6 7 8]
  sublhs = _mm256_sub_epi32(lhs_2, lhs_5_8);
  subrhs = _mm256_sub_epi32(rhs_5_8, rhs_2);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum7 = _mm256_add_epi64(accum7, mul);

  // 9 * [5 6 7 8]
  sublhs = _mm256_sub_epi32(lhs_9, lhs_5_8);
  subrhs = _mm256_sub_epi32(rhs_5_8, rhs_9);
  accum3 = _mm256_mul_epi32(sublhs, subrhs);

  __m256i lhs_4_7 = loadu_extend_32_64((__m128i*) &x->limbs[5]);
  __m256i rhs_4_7 = loadu_extend_32_64((__m128i*) &y->limbs[5]);

  __m256i lhs_3 = _mm256_set1_epi32(x->limbs[4]);
  __m256i rhs_3 = _mm256_set1_epi32(y->limbs[4]);

  // Third round is interleaved starting here
  __m256i lhs_6 = _mm256_set1_epi32(x->limbs[7]);
  __m256i rhs_6 = _mm256_set1_epi32(y->limbs[7]);

  // 3 * [4 5 6 7]
  sublhs = _mm256_sub_epi32(lhs_3, lhs_4_7);
  subrhs = _mm256_sub_epi32(rhs_4_7, rhs_3);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum7 = _mm256_add_epi64(accum7, mul);

  // 10 * [4 5 6 7]
  sublhs = _mm256_sub_epi32(lhs_10, lhs_4_7);
  subrhs = _mm256_sub_epi32(rhs_4_7, rhs_10);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum3 = _mm256_add_epi64(accum3, mul);

  // 6 * [4 5]
  sublhs128_low10 = _mm_sub_epi32(_mm256_castsi256_si128(lhs_6),
                                  _mm256_castsi256_si128(lhs_4_7));
  subrhs128_low10 = _mm_sub_epi32(_mm256_castsi256_si128(rhs_4_7),
                                  _mm256_castsi256_si128(rhs_6));

  __m256i lhs_3_6 = load_extend_32_64((__m128i*) &x->limbs[4]);
  __m256i rhs_3_6 = load_extend_32_64((__m128i*) &y->limbs[4]);

  __m256i lhs_7 = _mm256_set1_epi32(x->limbs[8]);
  __m256i rhs_7 = _mm256_set1_epi32(y->limbs[8]);

  // 0 * [3 4 5 6]
  sublhs = _mm256_sub_epi32(lhs_0, lhs_3_6);
  subrhs = _mm256_sub_epi32(rhs_3_6, rhs_0);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum3 = _mm256_add_epi64(accum3, mul);

  // 7 * [3 4 5 6]
  sublhs = _mm256_sub_epi32(lhs_7, lhs_3_6);
  subrhs = _mm256_sub_epi32(rhs_3_6, rhs_7);
  accum10 = _mm256_mul_epi32(sublhs, subrhs);

  // 4 * [5 6]
  sublhs = _mm256_sub_epi32(lhs_4, lhs_3_6);
  subrhs = _mm256_sub_epi32(rhs_3_6, rhs_4);
  sublhs = _mm256_blend_epi32(
      sublhs, _mm256_castsi128_si256(sublhs128_low7), 0x0f);
  subrhs = _mm256_blend_epi32(
      subrhs, _mm256_castsi128_si256(subrhs128_low7), 0x0f);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum7 = _mm256_add_epi64(accum7, mul);

  // Round 1 drops off here.
  _mm256_store_si256((__m256i*) &temp.limbs[8], accum7);

  __m256i lhs_2_5 = loadu_extend_32_64((__m128i*) &x->limbs[3]);
  __m256i rhs_2_5 = loadu_extend_32_64((__m128i*) &y->limbs[3]);

  // 1 * [2 3 4 5]
  sublhs = _mm256_sub_epi32(lhs_1, lhs_2_5);
  subrhs = _mm256_sub_epi32(rhs_2_5, rhs_1);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum3 = _mm256_add_epi64(accum3, mul);

  // 8 * [2 3 4 5]
  sublhs = _mm256_sub_epi32(lhs_8, lhs_2_5);
  subrhs = _mm256_sub_epi32(rhs_2_5, rhs_8);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum10 = _mm256_add_epi64(accum10, mul);

  __m256i lhs_1_4 = loadu_extend_32_64((__m128i*) &x->limbs[2]);
  __m256i rhs_1_4 = loadu_extend_32_64((__m128i*) &y->limbs[2]);

  // 9 * [1 2 3 4]
  sublhs = _mm256_sub_epi32(lhs_9, lhs_1_4);
  subrhs = _mm256_sub_epi32(rhs_1_4, rhs_9);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum10 = _mm256_add_epi64(accum10, mul);

  // 0 * [1 2]
  // We do this here, only because we never load [10 0 1 2]
  sublhs128_high = _mm_sub_epi32(_mm256_castsi256_si128(lhs_0),
                                 _mm256_castsi256_si128(lhs_1_4));
  subrhs128_high = _mm_sub_epi32(_mm256_castsi256_si128(rhs_1_4),
                                 _mm256_castsi256_si128(rhs_0));
  sublhs = _mm256_inserti128_si256(
    _mm256_castsi128_si256(sublhs128_low10),
    sublhs128_high, 1);
  subrhs = _mm256_inserti128_si256(
    _mm256_castsi128_si256(subrhs128_low10),
    subrhs128_high, 1);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum10 = _mm256_add_epi64(accum10, mul);

  // 2 * [3 4]
  sublhs = _mm256_sub_epi32(lhs_2, lhs_1_4);
  subrhs = _mm256_sub_epi32(rhs_1_4, rhs_2);
  sublhs = _mm256_blend_epi32(
      sublhs, _mm256_castsi128_si256(sublhs128_low3), 0x0f);
  subrhs = _mm256_blend_epi32(
      subrhs, _mm256_castsi128_si256(subrhs128_low3), 0x0f);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum3 = _mm256_add_epi64(accum3, mul);
  // Round 2 drops off here.
  _mm256_store_si256((__m256i*) &temp.limbs[4], accum3);

  __m256i lhs_0_3 = loadu_extend_32_64((__m128i*) &x->limbs[1]);
  __m256i rhs_0_3 = loadu_extend_32_64((__m128i*) &y->limbs[1]);

  // 10 * [0 1 2 3]
  sublhs = _mm256_sub_epi32(lhs_10, lhs_0_3);
  subrhs = _mm256_sub_epi32(rhs_0_3, rhs_10);
  mul = _mm256_mul_epi32(sublhs, subrhs);
  accum10 = _mm256_add_epi64(accum10, mul);
  // Round 3 drops off here.

  _mm256_store_si256((__m256i*) &temp.limbs[0], accum10);
  reduce_step_wide(&temp, &temp);
  reduce_step_wide(result, &temp);
}

// Multiply a wide residue by a small constant. The result is reduced to 32
// bits, but not narrowed for performance reasons.
void mul_wide_const(
  residue_wide_t *result, const residue_wide_t *x, int32_t d) {

  residue_wide_t temp;
  for (int i = 0; i < NLIMBS; ++i) {
    temp.limbs[i] = x->limbs[i] * d;
  }
  reduce_step_wide(result, &temp);
}

// Multiply a narrow residue by a small constant. The result is reduced to 32
// bits, but not narrowed for performance reasons.
void mul_narrow_const(
  residue_wide_t *result, const residue_narrow_t *x, int32_t d) {

  residue_wide_t temp;
  for (int i = 0; i < NLIMBS - 1; ++i) {
    temp.limbs[i] = ((uint64_t) x->limbs[i]) * d;
  }
  reduce_step_wide(result, &temp);
}

// Square a wide residue and produce a wide result. The result is reduced to 32
// bits but not narrowed for performance reasons.
void square_wide(
  residue_wide_t *result, const residue_wide_t *x) {

  residue_wide_t temp;

  __m256i sublhs, mul; // Temporaries for the actual sub sub mul
  __m128i sublhs128_low7;
  __m128i sublhs128_low3;
  __m128i sublhs128_low10;
  __m128i sublhs128_high; // Temporaries for single lane blends
  __m256i accum7, accum3, accum10; // Three accumulators
  accum7 = accum3 = accum10 = _mm256_setzero_si256();

  // We do 7-10 first.
  __m256i lhs_7_10 = _mm256_load_si256((__m256i*) &x->limbs[8]);
  __m256i lhs_0 = _mm256_set1_epi32(x->limbs[1]);

  // 0 * [7 8 9 10]
  sublhs = _mm256_sub_epi32(lhs_0, lhs_7_10);
  mul = _mm256_mul_epi32(sublhs, sublhs);
  accum7 = _mm256_sub_epi64(accum7, mul);

  __m256i lhs_6_9 = _mm256_loadu_si256((__m256i*) &x->limbs[7]);

  __m256i lhs_1 = _mm256_set1_epi32(x->limbs[2]);

  // Second round is interleaved starting here
  __m256i lhs_8 = _mm256_set1_epi32(x->limbs[9]);

  // 1 * [6 7 8 9]
  sublhs = _mm256_sub_epi32(lhs_1, lhs_6_9);
  mul = _mm256_mul_epi32(sublhs, sublhs);
  accum7 = _mm256_sub_epi64(accum7, mul);

  // 8 * [6 7]
  sublhs128_low3 = _mm_sub_epi32(_mm256_castsi256_si128(lhs_8),
                                 _mm256_castsi256_si128(lhs_6_9));

  __m256i lhs_5_8 = _mm256_loadu_si256((__m256i*) &x->limbs[6]);
  __m256i lhs_4 = _mm256_set1_epi32(x->limbs[5]);
  __m256i lhs_2 = _mm256_set1_epi32(x->limbs[3]);
  __m256i lhs_9 = _mm256_set1_epi32(x->limbs[10]);

  // 2 * [5 6 7 8]
  sublhs = _mm256_sub_epi32(lhs_2, lhs_5_8);
  mul = _mm256_mul_epi32(sublhs, sublhs);
  accum7 = _mm256_sub_epi64(accum7, mul);

  // 9 * [5 6 7 8]
  sublhs = _mm256_sub_epi32(lhs_9, lhs_5_8);
  mul = _mm256_mul_epi32(sublhs, sublhs);
  accum3 = _mm256_sub_epi64(accum3, mul);

  __m256i lhs_4_7 = _mm256_loadu_si256((__m256i*) &x->limbs[5]);
  __m256i lhs_10 = _mm256_set1_epi32(x->limbs[0]);
  __m256i lhs_3 = _mm256_set1_epi32(x->limbs[4]);

  // Third round is interleaved starting here
  __m256i lhs_6 = _mm256_set1_epi32(x->limbs[7]);

  // 3 * [4 5 6 7]
  sublhs = _mm256_sub_epi32(lhs_3, lhs_4_7);
  mul = _mm256_mul_epi32(sublhs, sublhs);
  accum7 = _mm256_sub_epi64(accum7, mul);

  // 10 * [4 5 6 7]
  sublhs = _mm256_sub_epi32(lhs_10, lhs_4_7);
  mul = _mm256_mul_epi32(sublhs, sublhs);
  accum3 = _mm256_sub_epi64(accum3, mul);

  // 6 * [4 5]
  sublhs128_low10 = _mm_sub_epi32(_mm256_castsi256_si128(lhs_6),
                                  _mm256_castsi256_si128(lhs_4_7));

  __m256i lhs_3_6 = _mm256_load_si256((__m256i*) &x->limbs[4]);
  __m256i lhs_7 = _mm256_set1_epi32(x->limbs[8]);

  // 0 * [3 4 5 6]
  sublhs = _mm256_sub_epi32(lhs_0, lhs_3_6);
  mul = _mm256_mul_epi32(sublhs, sublhs);
  accum3 = _mm256_sub_epi64(accum3, mul);

  // 7 * [3 4 5 6]
  sublhs = _mm256_sub_epi32(lhs_7, lhs_3_6);
  mul = _mm256_mul_epi32(sublhs, sublhs);
  accum10 = _mm256_sub_epi64(accum10, mul);

  // 4 * [5 6], 10 * [8 9]
  __m128i lhs_8_9 = _mm256_extracti128_si256(lhs_6_9, 1);

  sublhs128_low7 = _mm_sub_epi32(_mm256_castsi256_si128(lhs_10), lhs_8_9);

  sublhs = _mm256_sub_epi32(lhs_4, lhs_3_6);
  sublhs = _mm256_blend_epi32(
      sublhs, _mm256_castsi128_si256(sublhs128_low7), 0x0f);
  mul = _mm256_mul_epi32(sublhs, sublhs);
  accum7 = _mm256_sub_epi64(accum7, mul);

  // Round 1 drops off here.
  _mm256_store_si256((__m256i*) &temp.limbs[8], accum7);

  __m256i lhs_2_5 = _mm256_loadu_si256((__m256i*) &x->limbs[3]);

  // 1 * [2 3 4 5]
  sublhs = _mm256_sub_epi32(lhs_1, lhs_2_5);
  mul = _mm256_mul_epi32(sublhs, sublhs);
  accum3 = _mm256_sub_epi64(accum3, mul);

  // 8 * [2 3 4 5]
  sublhs = _mm256_sub_epi32(lhs_8, lhs_2_5);
  mul = _mm256_mul_epi32(sublhs, sublhs);
  accum10 = _mm256_sub_epi64(accum10, mul);

  __m256i lhs_1_4 = _mm256_loadu_si256((__m256i*) &x->limbs[2]);

  // 9 * [1 2 3 4]
  sublhs = _mm256_sub_epi32(lhs_9, lhs_1_4);
  mul = _mm256_mul_epi32(sublhs, sublhs);
  accum10 = _mm256_sub_epi64(accum10, mul);

  // 0 * [1 2]
  // We do this here, only because we never load [10 0 1 2]
  sublhs128_high = _mm_sub_epi32(_mm256_castsi256_si128(lhs_0),
                                 _mm256_castsi256_si128(lhs_1_4));
  sublhs = _mm256_inserti128_si256(
    _mm256_castsi128_si256(sublhs128_low10),
    sublhs128_high, 1);
  mul = _mm256_mul_epi32(sublhs, sublhs);
  accum10 = _mm256_sub_epi64(accum10, mul);

  // 2 * [3 4]
  sublhs = _mm256_sub_epi32(lhs_2, lhs_1_4);
  sublhs = _mm256_blend_epi32(
      sublhs, _mm256_castsi128_si256(sublhs128_low3), 0x0f);
  mul = _mm256_mul_epi32(sublhs, sublhs);
  accum3 = _mm256_sub_epi64(accum3, mul);
  // Round 2 drops off here.
  _mm256_store_si256((__m256i*) &temp.limbs[4], accum3);

  __m256i lhs_0_3 = _mm256_loadu_si256((__m256i*) &x->limbs[1]);

  // 10 * [0 1 2 3]
  sublhs = _mm256_sub_epi32(lhs_10, lhs_0_3);
  mul = _mm256_mul_epi32(sublhs, sublhs);
  accum10 = _mm256_sub_epi64(accum10, mul);
  // Round 3 drops off here.

  _mm256_store_si256((__m256i*) &temp.limbs[0], accum10);
  reduce_step_wide(&temp, &temp);
  reduce_step_wide(result, &temp);
}

// Square a narrow residue and produce a wide result. The result is reduced to
// 32 bits but not narrowed for performance reasons.
void square_narrow(
  residue_wide_t *result, const residue_narrow_t *x) {

  residue_wide_t temp;

  __m256i sublhs, mul; // Temporaries for the actual sub sub mul
  __m128i sublhs128_low7;
  __m128i sublhs128_low3;
  __m128i sublhs128_low10;
  __m128i sublhs128_high; // Temporaries for single lane blends
  __m256i accum7, accum3, accum10; // Three accumulators
  accum7 = accum3 = accum10 = _mm256_setzero_si256();

  // We do 7-10 first.
  __m256i lhs_10 = _mm256_set1_epi32(x->limbs[0]);
  __m256i lhs_7_10 = load_extend_32_64((__m128i*) &x->limbs[8]);
  __m256i lhs_0 = _mm256_set1_epi32(x->limbs[1]);

  // 0 * [7 8 9 10]
  sublhs = _mm256_sub_epi32(lhs_0, lhs_7_10);
  mul = _mm256_mul_epi32(sublhs, sublhs);
  accum7 = _mm256_sub_epi64(accum7, mul);

  __m256i lhs_6_9 = loadu_extend_32_64((__m128i*) &x->limbs[7]);

  // 10 * [8 9]
  __m128i lhs_8_9 = _mm256_extracti128_si256(lhs_6_9, 1);

  sublhs128_low7 = _mm_sub_epi32(_mm256_castsi256_si128(lhs_10), lhs_8_9);

  __m256i lhs_1 = _mm256_set1_epi32(x->limbs[2]);

  // Second round is interleaved starting here
  __m256i lhs_8 = _mm256_set1_epi32(x->limbs[9]);

  // 1 * [6 7 8 9]
  sublhs = _mm256_sub_epi32(lhs_1, lhs_6_9);
  mul = _mm256_mul_epi32(sublhs, sublhs);
  accum7 = _mm256_sub_epi64(accum7, mul);

  // 8 * [6 7]
  sublhs128_low3 = _mm_sub_epi32(_mm256_castsi256_si128(lhs_8),
                                 _mm256_castsi256_si128(lhs_6_9));

  __m256i lhs_5_8 = loadu_extend_32_64((__m128i*) &x->limbs[6]);
  __m256i lhs_4 = _mm256_set1_epi32(x->limbs[5]);
  __m256i lhs_2 = _mm256_set1_epi32(x->limbs[3]);
  __m256i lhs_9 = _mm256_set1_epi32(x->limbs[10]);

  // 2 * [5 6 7 8]
  sublhs = _mm256_sub_epi32(lhs_2, lhs_5_8);
  mul = _mm256_mul_epi32(sublhs, sublhs);
  accum7 = _mm256_sub_epi64(accum7, mul);

  // 9 * [5 6 7 8]
  sublhs = _mm256_sub_epi32(lhs_9, lhs_5_8);
  mul = _mm256_mul_epi32(sublhs, sublhs);
  accum3 = _mm256_sub_epi64(accum3, mul);

  __m256i lhs_4_7 = loadu_extend_32_64((__m128i*) &x->limbs[5]);
  __m256i lhs_3 = _mm256_set1_epi32(x->limbs[4]);

  // Third round is interleaved starting here
  __m256i lhs_6 = _mm256_set1_epi32(x->limbs[7]);

  // 3 * [4 5 6 7]
  sublhs = _mm256_sub_epi32(lhs_3, lhs_4_7);
  mul = _mm256_mul_epi32(sublhs, sublhs);
  accum7 = _mm256_sub_epi64(accum7, mul);

  // 10 * [4 5 6 7]
  sublhs = _mm256_sub_epi32(lhs_10, lhs_4_7);
  mul = _mm256_mul_epi32(sublhs, sublhs);
  accum3 = _mm256_sub_epi64(accum3, mul);

  // 6 * [4 5]
  sublhs128_low10 = _mm_sub_epi32(_mm256_castsi256_si128(lhs_6),
                                  _mm256_castsi256_si128(lhs_4_7));

  __m256i lhs_3_6 = load_extend_32_64((__m128i*) &x->limbs[4]);
  __m256i lhs_7 = _mm256_set1_epi32(x->limbs[8]);

  // 0 * [3 4 5 6]
  sublhs = _mm256_sub_epi32(lhs_0, lhs_3_6);
  mul = _mm256_mul_epi32(sublhs, sublhs);
  accum3 = _mm256_sub_epi64(accum3, mul);

  // 7 * [3 4 5 6]
  sublhs = _mm256_sub_epi32(lhs_7, lhs_3_6);
  mul = _mm256_mul_epi32(sublhs, sublhs);
  accum10 = _mm256_sub_epi64(accum10, mul);

  // 4 * [5 6]
  sublhs = _mm256_sub_epi32(lhs_4, lhs_3_6);
  sublhs = _mm256_blend_epi32(
      sublhs, _mm256_castsi128_si256(sublhs128_low7), 0x0f);
  mul = _mm256_mul_epi32(sublhs, sublhs);
  accum7 = _mm256_sub_epi64(accum7, mul);

  // Round 1 drops off here.
  _mm256_store_si256((__m256i*) &temp.limbs[8], accum7);

  __m256i lhs_2_5 = loadu_extend_32_64((__m128i*) &x->limbs[3]);

  // 1 * [2 3 4 5]
  sublhs = _mm256_sub_epi32(lhs_1, lhs_2_5);
  mul = _mm256_mul_epi32(sublhs, sublhs);
  accum3 = _mm256_sub_epi64(accum3, mul);

  // 8 * [2 3 4 5]
  sublhs = _mm256_sub_epi32(lhs_8, lhs_2_5);
  mul = _mm256_mul_epi32(sublhs, sublhs);
  accum10 = _mm256_sub_epi64(accum10, mul);

  __m256i lhs_1_4 = loadu_extend_32_64((__m128i*) &x->limbs[2]);

  // 9 * [1 2 3 4]
  sublhs = _mm256_sub_epi32(lhs_9, lhs_1_4);
  mul = _mm256_mul_epi32(sublhs, sublhs);
  accum10 = _mm256_sub_epi64(accum10, mul);

  // 0 * [1 2]
  // We do this here, only because we never load [10 0 1 2]
  sublhs128_high = _mm_sub_epi32(_mm256_castsi256_si128(lhs_0),
                                 _mm256_castsi256_si128(lhs_1_4));
  sublhs = _mm256_inserti128_si256(
    _mm256_castsi128_si256(sublhs128_low10),
    sublhs128_high, 1);
  mul = _mm256_mul_epi32(sublhs, sublhs);
  accum10 = _mm256_sub_epi64(accum10, mul);

  // 2 * [3 4]
  sublhs = _mm256_sub_epi32(lhs_2, lhs_1_4);
  sublhs = _mm256_blend_epi32(
      sublhs, _mm256_castsi128_si256(sublhs128_low3), 0x0f);
  mul = _mm256_mul_epi32(sublhs, sublhs);
  accum3 = _mm256_sub_epi64(accum3, mul);
  // Round 2 drops off here.
  _mm256_store_si256((__m256i*) &temp.limbs[4], accum3);

  __m256i lhs_0_3 = loadu_extend_32_64((__m128i*) &x->limbs[1]);

  // 10 * [0 1 2 3]
  sublhs = _mm256_sub_epi32(lhs_10, lhs_0_3);
  mul = _mm256_mul_epi32(sublhs, sublhs);
  accum10 = _mm256_sub_epi64(accum10, mul);
  // Round 3 drops off here.

  _mm256_store_si256((__m256i*) &temp.limbs[0], accum10);
  reduce_step_wide(&temp, &temp);
  reduce_step_wide(result, &temp);
}

// Approximately divide each coefficient by t. Carry the results.
void reduce_step_narrow(
  residue_narrow_t *result, const residue_narrow_t *x) {

  __m256i accum0, error0, shift_error0, carry_rot0;
  __m128i accum8, error8, shift_error8, carry_rot8;

  __m256i mask = _mm256_set1_epi32(0x3ffffff);
  uint32_t slide_permute[8] = { 7, 0, 1, 2, 3, 4, 5, 6 };
  __m256i slide_permute_v = _mm256_loadu_si256((__m256i*) slide_permute);

  accum8 = _mm_load_si128((__m128i*) &x->limbs[8]);
  accum0 = _mm256_loadu_si256((__m256i*) &x->limbs[0]);

  error8 = _mm_srai_epi32(accum8, 26);
  carry_rot8 = _mm_shuffle_epi32(error8, 0x92);
  shift_error8 = _mm_slli_epi32(error8, 4);

  error0 = _mm256_srai_epi32(accum0, 26);
  carry_rot0 = _mm256_permutevar8x32_epi32(error0, slide_permute_v);
  shift_error0 = _mm256_slli_epi32(error0, 4);

  accum8 = _mm_and_si128(accum8, _mm256_castsi256_si128(mask));
  accum0 = _mm256_and_si256(accum0, mask);

  accum8 = _mm_sub_epi32(accum8, error8);
  accum0 = _mm256_sub_epi32(accum0, error0);

  accum8 = _mm_add_epi32(accum8, shift_error8);
  accum0 = _mm256_add_epi32(accum0, shift_error0);

  __m128i merged_carry8 = _mm_blend_epi32(
      carry_rot8, _mm256_castsi256_si128(carry_rot0), 0x01);
  accum8 = _mm_add_epi32(accum8, merged_carry8);
  __m256i merged_carry0 = _mm256_blend_epi32(
    carry_rot0, _mm256_castsi128_si256(carry_rot8), 0x01);
  accum0 = _mm256_add_epi32(accum0, merged_carry0);
  _mm_store_si128((__m128i*) &result->limbs[8], accum8);
  _mm256_storeu_si256((__m256i*) &result->limbs[0], accum0);
}

// Approximately divide each coefficient by t. Carry the results.
void reduce_step_wide(
  residue_wide_t *result, const residue_wide_t *x) {

  __m256i accum0, accum4, accum8;

  __m256i logical_shift;
  __m256i arithmetic_shift;
  __m256i carry_rot0, carry_rot4, carry_rot8;
  __m256i error0, error4, error8;
  __m256i shift_error0, shift_error4, shift_error8;
  __m256i merged_carry;
  __m256i mask = _mm256_set1_epi64x(0x3ffffff);

  accum8 = _mm256_load_si256((__m256i*) &x->limbs[8]);
  logical_shift = _mm256_srli_epi64(accum8, 26);
  arithmetic_shift = _mm256_srai_epi32(accum8, 26);
  accum8 = _mm256_and_si256(accum8, mask);
  error8 = _mm256_blend_epi32(logical_shift, arithmetic_shift, 0xaa);
  accum8 = _mm256_sub_epi64(accum8, error8);
  carry_rot8 = _mm256_permute4x64_epi64(error8, 0x92);
  shift_error8 = _mm256_slli_epi64(error8, 4);
  accum8 = _mm256_add_epi64(accum8, shift_error8);

  accum4 = _mm256_load_si256((__m256i*) &x->limbs[4]);
  logical_shift = _mm256_srli_epi64(accum4, 26);
  arithmetic_shift = _mm256_srai_epi32(accum4, 26);
  accum4 = _mm256_and_si256(accum4, mask);
  error4 = _mm256_blend_epi32(logical_shift, arithmetic_shift, 0xaa);
  accum4 = _mm256_sub_epi64(accum4, error4);
  carry_rot4 = _mm256_permute4x64_epi64(error4, 0x93);
  shift_error4 = _mm256_slli_epi64(error4, 4);
  accum4 = _mm256_add_epi64(accum4, shift_error4);

  merged_carry = _mm256_blend_epi32(carry_rot8, carry_rot4, 0x03);
  accum8 = _mm256_add_epi64(accum8, merged_carry);
  _mm256_store_si256((__m256i*) &result->limbs[8], accum8);

  accum0 = _mm256_load_si256((__m256i*) &x->limbs[0]);
  logical_shift = _mm256_srli_epi64(accum0, 26);
  arithmetic_shift = _mm256_srai_epi32(accum0, 26);
  accum0 = _mm256_and_si256(accum0, mask);
  error0 = _mm256_blend_epi32(logical_shift, arithmetic_shift, 0xaa);
  accum0 = _mm256_sub_epi64(accum0, error0);
  carry_rot0 = _mm256_permute4x64_epi64(error0, 0x93);

  merged_carry = _mm256_blend_epi32(carry_rot4, carry_rot0, 0x03);
  accum4 = _mm256_add_epi64(accum4, merged_carry);
  _mm256_store_si256((__m256i*) &result->limbs[4], accum4);

  shift_error0 = _mm256_slli_epi64(error0, 4);
  accum0 = _mm256_add_epi64(accum0, shift_error0);

  merged_carry = _mm256_blend_epi32(carry_rot0, carry_rot8, 0x03);
  accum0 = _mm256_add_epi64(accum0, merged_carry);
  _mm256_store_si256((__m256i*) &result->limbs[0], accum0);
}

// Takes advantage of the fact that if a residue z *is zero* then after setting
// one coefficient to T/2, all the remaining coefficients should be near to
// T/2. They should therefore resolve all carries in a single step, and all be
// equal to the same value. Some other value may not reduce completely, but this
// is fine, we will know it is not zero.
int equal_wide(const residue_wide_t *x, const residue_wide_t *y) {
  residue_wide_t temp;

  sub_wide(&temp, x, y);
  int64_t delta = -temp.limbs[0] + (T / 2);
  for (int i = 0; i < NLIMBS; ++i) {
    temp.limbs[i] += delta;
  }

  reduce_step_wide(&temp, &temp);

  delta = temp.limbs[0];
  int result = 0;
  for (int i = 1; i < NLIMBS; ++i) {
    result |= (temp.limbs[i] ^ delta);
  }

  return !result;
}

int equal_narrow_reduced(
  const residue_narrow_reduced_t * x, const residue_narrow_reduced_t * y) {

  int result = 0;
  for (int i = 0; i < NLIMBS_REDUCED; ++i) {
    result |= (x->limbs[i] ^ y->limbs[i]);
  }

  return !result;
}

static inline void nsquare_wide(
  residue_wide_t *result, const residue_wide_t *x, int n) {

  square_wide(result, x);
  for (int i = 1; i < n; ++i) {
    square_wide(result, result);
  }
}

static void raise_to_t(
  residue_wide_t *result, const residue_wide_t *x) {
  // zi = z^(2^i - 1), z1 = x
  residue_wide_t z2;
  residue_wide_t z3;
  residue_wide_t z5;
  residue_wide_t z10;
  residue_wide_t z11;
  residue_wide_t z22;
  residue_wide_t result_t;

  square_wide(&z2, x);
  mul_wide(&z2, &z2, x);
  square_wide(&z3, &z2);
  mul_wide(&z3, &z3, x);
  nsquare_wide(&z5, &z3, 2);
  mul_wide(&z5, &z5, &z2);
  nsquare_wide(&z10, &z5, 5);
  mul_wide(&z10, &z10, &z5);
  square_wide(&z11, &z10);
  mul_wide(&z11, &z11, x);
  nsquare_wide(&z22, &z11, 11);
  mul_wide(&z22, &z22, &z11);
  nsquare_wide(&result_t, &z22, 4);
  mul_wide(result, &result_t, x);
}

static void raise_to_phi_t(
  residue_wide_t *result, const residue_wide_t *x, int n) {
  residue_wide_t temp;

  raise_to_t(&temp, x);

  for (int i = 1; i < n; ++i) {
    mul_wide(&temp, &temp, x);
    raise_to_t(&temp, &temp);
  }

  mul_wide(result, &temp, x);
}

static void raise_to_t_minus_1_over_4(
  residue_wide_t *result, const residue_wide_t *x) {
  // zi = z^(2^i - 1), z1 = x
  residue_wide_t z2;
  residue_wide_t z3;
  residue_wide_t z5;
  residue_wide_t z10;
  residue_wide_t z11;
  residue_wide_t z22;

  square_wide(&z2, x);
  mul_wide(&z2, &z2, x);
  square_wide(&z3, &z2);
  mul_wide(&z3, &z3, x);
  nsquare_wide(&z5, &z3, 2);
  mul_wide(&z5, &z5, &z2);
  nsquare_wide(&z10, &z5, 5);
  mul_wide(&z10, &z10, &z5);
  square_wide(&z11, &z10);
  mul_wide(&z11, &z11, x);
  nsquare_wide(&z22, &z11, 11);
  mul_wide(&z22, &z22, &z11);
  nsquare_wide(result, &z22, 2);
}

static void raise_to_p_minus_3_over_4(
  residue_wide_t *result, const residue_wide_t *x) {

  residue_wide_t z4; //z to (t-1)/4
  residue_wide_t z2; //z to (t-1)/2
  residue_wide_t z3_4; //z to (3t+1)/4
  residue_wide_t y_small;
  residue_wide_t y, y_t4_y;
  residue_wide_t raised;

  raise_to_t_minus_1_over_4(&z4, x);
  square_wide(&z2, &z4);
  mul_wide(&z3_4, &z2, &z4);
  mul_wide(&z3_4, &z3_4, x);
  raise_to_t(&raised, &z4);
  mul_wide(&y_small, &z2, &raised);
  raise_to_t(&raised, &y_small);
  mul_wide(&y, &z3_4, &raised);
  raise_to_t(&raised, &y);
  raise_to_t(&raised, &raised);
  raise_to_t(&raised, &raised);
  raise_to_t(&raised, &raised);
  mul_wide(&y_t4_y, &raised, &y);
  raise_to_t(&raised, &y_t4_y);
  raise_to_t(&raised, &raised);
  raise_to_t(&raised, &raised);
  mul_wide(result, &raised, &y_small);
}

int sqrt_inv_wide(
  residue_wide_t *result, const residue_wide_t * __restrict x,
  const residue_wide_t * __restrict y) {
  residue_wide_t xy;
  residue_wide_t y2;
  residue_wide_t xy3;
  residue_wide_t xy3_p_3_over_4;
  residue_wide_t cand2;
  residue_wide_t should_be_x;

  square_wide(&y2, y);
  mul_wide(&xy, x, y);
  mul_wide(&xy3, &xy, &y2);
  raise_to_p_minus_3_over_4(&xy3_p_3_over_4, &xy3);
  mul_wide(result, &xy, &xy3_p_3_over_4);
  square_wide(&cand2, result);
  mul_wide(&should_be_x, y, &cand2);

  return equal_wide(&should_be_x, x);
}

void invert_wide(
  residue_wide_t *result, const residue_wide_t * __restrict x) {

  residue_wide_t x_t_minus_1_over_4;
  residue_wide_t x_t_minus_1;
  residue_wide_t x_t;
  residue_wide_t phi_8_x_t;
  residue_wide_t phi_8_x_t_t;

  raise_to_t_minus_1_over_4(&x_t_minus_1_over_4, x);
  nsquare_wide(&x_t_minus_1, &x_t_minus_1_over_4, 2);
  mul_wide(&x_t, &x_t_minus_1, x);
  raise_to_phi_t(&phi_8_x_t, &x_t, 8);
  raise_to_t(&phi_8_x_t_t, &phi_8_x_t);
  mul_wide(result, &phi_8_x_t_t, &x_t_minus_1);
}

void encode(uint8_t *out, const residue_narrow_reduced_t * __restrict x) {
  uint32_t collect = x->limbs[0];

  int space = 32 - TBITS;
  int i = 1;
  int bits_remaining = TBITS * NLIMBS_REDUCED;
  while (bits_remaining > 0) {
    *out++ = collect & 0xff;
    collect >>= 8;
    space += 8;
    bits_remaining -= 8;
    if (space >= TBITS && i < NLIMBS_REDUCED) {
      collect |= x->limbs[i] << (32 - space);
      space -= TBITS;
      ++i;
    }
  }
}

void decode(residue_narrow_reduced_t *out, const uint8_t *in) {
  uint32_t collect = 0;

  int shift = 0;
  int i = 0;
  int bits_remaining = TBITS * NLIMBS_REDUCED;
  while (bits_remaining > 0) {
    collect |= (*in++) << shift;
    shift += 8;
    bits_remaining -= 8;
    if (shift >= TBITS) {
      if (bits_remaining > 0) {
        out->limbs[i] = collect & TMASK;
        collect >>= 26;
        shift -= 26;
        ++i;
      } else {
        out->limbs[i] = collect;
      }
    }
  }
}
