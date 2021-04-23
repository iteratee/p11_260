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

#define FORCE_BLEND 0

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

static inline __m256i load_extend_32_64(__m128i *x) {
  return _mm256_cvtepi32_epi64(_mm_load_si128(x));
}

static inline __m256i loadu_extend_32_64(__m128i *x) {
  return _mm256_cvtepi32_epi64(_mm_loadu_si128(x));
}

static inline __m512i load512_extend_32_64(__m256i *x) {
  return _mm512_cvtepi32_epi64(_mm256_load_si256(x));
}

static inline __m512i loadu512_extend_32_64(__m256i *x) {
  return _mm512_cvtepi32_epi64(_mm256_loadu_si256(x));
}

static inline __m512i loadu512_mask_extend_32_64(__m256i *x, __mmask8 k) {
  return _mm512_cvtepi32_epi64(
    _mm256_mask_loadu_epi32(_mm256_setzero_si256(), k, x));
}

// Produce a 64-bit residue
void widen(
  residue_wide_t *result, const residue_narrow_t * __restrict x) {
  __m256i wide10 = loadu_extend_32_64((__m128i *) x);
  __m512i wide3 = loadu512_extend_32_64((__m256i *) &x[4]);
  _mm256_store_si256((__m256i*) &result->limbs[0], wide10);
  _mm512_storeu_si512((__m512i*) &result->limbs[4], wide3);
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

#include <stdio.h>
#include <string.h>
// static void print4x64(__m256i x, const char * preamble) {
//   uint64_t x_vals[4];
//   memcpy(x_vals, &x, sizeof(x_vals));
//   printf("%s\n", preamble);
//   for (int i = 0; i < 4; ++i) {
//     printf("%#lx\n", x_vals[i]);
//   }
// }

// static void print8x64(__m512i x, const char * preamble) {
//   uint64_t x_vals[8];
//   memcpy(x_vals, &x, sizeof(x_vals));
//   printf("%s\n", preamble);
//   for (int i = 0; i < 8; ++i) {
//     printf("%#lx\n", x_vals[i]);
//   }
// }

#define wrap(x) (((x + (NLIMBS - 1)) % (NLIMBS - 1)) + 1)
// Multiply two wide residues, and produce a wide result. The result is reduced
// to 32 bits, but not narrowed for performance reasons.
void mul_wide(
  residue_wide_t *result, const residue_wide_t *x, const residue_wide_t *y) {

  residue_wide_t temp;

  // Two accumulators
  __m256i accum10;
  __m512i accum3;

  // Temporaries for sub sub mul
  __m512i sublhs_512, subrhs_512, mul_512;
  __m256i sublhs_256, subrhs_256, mul_256;

  __m512i blend_var_lhs, blend_var_rhs;
  __m512i blend_const_lhs, blend_const_rhs;

  __m256i blend_var_lhs_256, blend_var_rhs_256;
  __m256i blend_const_lhs_256, blend_const_rhs_256;

  __mmask16 low4, low8, low12;

  low4 = _cvtu32_mask16(0xf);
  low8 = _cvtu32_mask16(0xff);
  low12 = _cvtu32_mask16(0xfff);

  // 8 * [6 7] and 2 * [3 4 5 6 7 8]
  __m512i lhs_6_7 = _mm512_castsi128_si512(_mm_loadu_si128((__m128i *) &x->limbs[7]));
  __m512i rhs_6_7 = _mm512_castsi128_si512(_mm_loadu_si128((__m128i *) &y->limbs[7]));
//   print8x64(lhs_6_7, "lhs_6_7");
//   print8x64(rhs_6_7, "rhs_6_7");

  __m512i lhs_1_8 = _mm512_loadu_si512((__m512i *) &x->limbs[2]);
  __m512i rhs_1_8 = _mm512_loadu_si512((__m512i *) &y->limbs[2]);

//   print8x64(lhs_1_8, "lhs_1_8");
//   print8x64(rhs_1_8, "rhs_1_8");

  __m512i lhs_8 = _mm512_set1_epi32(x->limbs[9]);
  __m512i rhs_8 = _mm512_set1_epi32(y->limbs[9]);

//   print8x64(lhs_8, "lhs_8");
//   print8x64(rhs_8, "rhs_8");

  __m512i lhs_2 = _mm512_set1_epi32(x->limbs[3]);
  __m512i rhs_2 = _mm512_set1_epi32(y->limbs[3]);

//   print8x64(lhs_2, "lhs_2");
//   print8x64(rhs_2, "rhs_2");

#if FORCE_BLEND
  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_var_lhs)
      : [mask]  "Yk"((__mmask16) low4),
        [a]     "v"(lhs_6_7),
        [b]     "v"(lhs_1_8));
  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_var_rhs)
      : [mask]  "Yk"((__mmask16) low4),
        [a]     "v"(rhs_6_7),
        [b]     "v"(rhs_1_8));

  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_const_lhs)
      : [mask]  "Yk"((__mmask16) low4),
        [a]     "v"(lhs_8),
        [b]     "v"(lhs_2));
  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_const_rhs)
      : [mask]  "Yk"((__mmask16) low4),
        [a]     "v"(rhs_8),
        [b]     "v"(rhs_2));
#else
  blend_var_lhs = _mm512_mask_blend_epi32(low4, lhs_6_7, lhs_1_8);
  blend_var_rhs = _mm512_mask_blend_epi32(low4, rhs_6_7, rhs_1_8);

  blend_const_lhs = _mm512_mask_blend_epi32(low4, lhs_8, lhs_2);
  blend_const_rhs = _mm512_mask_blend_epi32(low4, rhs_8, rhs_2);
#endif

//   print8x64(blend_var_lhs, "blend_var_lhs");
//   print8x64(blend_var_rhs, "blend_var_rhs");
//   print8x64(blend_const_lhs, "blend_const_lhs");
//   print8x64(blend_const_rhs, "blend_const_rhs");

  sublhs_512 = _mm512_sub_epi32(blend_const_lhs, blend_var_lhs);
  subrhs_512 = _mm512_sub_epi32(blend_var_rhs, blend_const_rhs);
//   print8x64(sublhs_512, "sublhs_512");
//   print8x64(subrhs_512, "subrhs_512");

  accum3 = _mm512_mul_epi32(sublhs_512, subrhs_512);
//   print8x64(accum3, "accum3");

  // 9 * [5 6 7 8] and 3 * [4 5 6 7]
  __m512i lhs_5_8 = _mm512_castsi256_si512(_mm256_loadu_si256((__m256i *) &x->limbs[6]));
  __m512i rhs_5_8 = _mm512_castsi256_si512(_mm256_loadu_si256((__m256i *) &y->limbs[6]));

//   print8x64(lhs_5_8, "lhs_5_8");
//   print8x64(rhs_5_8, "rhs_5_8");

  __m512i lhs_0_7 = _mm512_loadu_si512((__m512i *) &x->limbs[1]);
  __m512i rhs_0_7 = _mm512_loadu_si512((__m512i *) &y->limbs[1]);

//   print8x64(lhs_0_7, "lhs_0_7");
//   print8x64(rhs_0_7, "rhs_0_7");

  __m512i lhs_9 = _mm512_set1_epi32(x->limbs[10]);
  __m512i rhs_9 = _mm512_set1_epi32(y->limbs[10]);

//   print8x64(lhs_9, "lhs_9");
//   print8x64(rhs_9, "rhs_9");

  __m512i lhs_3 = _mm512_set1_epi32(x->limbs[4]);
  __m512i rhs_3 = _mm512_set1_epi32(y->limbs[4]);

//   print8x64(lhs_3, "lhs_3");
//   print8x64(rhs_3, "rhs_3");

#if FORCE_BLEND
  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_var_lhs)
      : [mask]  "Yk"((__mmask16) low8),
        [a]     "v"(lhs_5_8),
        [b]     "v"(lhs_0_7));
  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_var_rhs)
      : [mask]  "Yk"((__mmask16) low8),
        [a]     "v"(rhs_5_8),
        [b]     "v"(rhs_0_7));

  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_const_lhs)
      : [mask]  "Yk"((__mmask16) low8),
        [a]     "v"(lhs_9),
        [b]     "v"(lhs_3));
  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_const_rhs)
      : [mask]  "Yk"((__mmask16) low8),
        [a]     "v"(rhs_9),
        [b]     "v"(rhs_3));
#else
  blend_var_lhs = _mm512_mask_blend_epi32(low8, lhs_5_8, lhs_0_7);
  blend_var_rhs = _mm512_mask_blend_epi32(low8, rhs_5_8, rhs_0_7);

  blend_const_lhs = _mm512_mask_blend_epi32(low8, lhs_9, lhs_3);
  blend_const_rhs = _mm512_mask_blend_epi32(low8, rhs_9, rhs_3);
#endif

//   print8x64(blend_var_lhs, "blend_var_lhs");
//   print8x64(blend_var_rhs, "blend_var_rhs");
//   print8x64(blend_const_lhs, "blend_const_lhs");
//   print8x64(blend_const_rhs, "blend_const_rhs");

  sublhs_512 = _mm512_sub_epi32(blend_const_lhs, blend_var_lhs);
  subrhs_512 = _mm512_sub_epi32(blend_var_rhs, blend_const_rhs);
//   print8x64(sublhs_512, "sublhs_512");
//   print8x64(subrhs_512, "subrhs_512");

  mul_512 = _mm512_mul_epi32(sublhs_512, subrhs_512);
//   print8x64(mul_512, "mul_512");
  accum3 = _mm512_add_epi64(accum3, mul_512);
//   print8x64(accum3, "accum3");

  // 9 * [1 2 3 4] Done here because we have loaded everything already.
  sublhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(lhs_9),
                                _mm512_castsi512_si256(lhs_1_8));
  subrhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(rhs_1_8),
                                _mm512_castsi512_si256(rhs_9));
//   print4x64(sublhs_256, "sublhs_256");
//   print4x64(subrhs_256, "subrhs_256");

  accum10 = _mm256_mul_epi32(sublhs_256, subrhs_256);
//   print4x64(accum10, "accum10");

  // 10 * [4 5 6 7 8 9] and 4 * [5 6]
  __m512i lhs_4_9 = _mm512_mask_loadu_epi32(_mm512_setzero(), low12, &x->limbs[5]);
  __m512i rhs_4_9 = _mm512_mask_loadu_epi32(_mm512_setzero(), low12, &y->limbs[5]);

//   print8x64(lhs_4_9, "lhs_4_9");
//   print8x64(rhs_4_9, "rhs_4_9");

  __m512i lhs_10_6 = _mm512_loadu_si512((__m512i *) &x->limbs[0]);
  __m512i rhs_10_6 = _mm512_loadu_si512((__m512i *) &y->limbs[0]);

//   print8x64(lhs_10_6, "lhs_10_6");
//   print8x64(rhs_10_6, "rhs_10_6");

  __m512i lhs_10 = _mm512_set1_epi32(x->limbs[0]);
  __m512i rhs_10 = _mm512_set1_epi32(y->limbs[0]);

//   print8x64(lhs_10, "lhs_10");
//   print8x64(rhs_10, "rhs_10");

  __m512i lhs_4 = _mm512_set1_epi32(x->limbs[5]);
  __m512i rhs_4 = _mm512_set1_epi32(y->limbs[5]);

//   print8x64(lhs_4, "lhs_4");
//   print8x64(rhs_4, "rhs_4");

#if FORCE_BLEND
  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_var_lhs)
      : [mask]  "Yk"((__mmask16) low12),
        [a]     "v"(lhs_4_9),
        [b]     "v"(lhs_10_6));
  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_var_rhs)
      : [mask]  "Yk"((__mmask16) low12),
        [a]     "v"(rhs_4_9),
        [b]     "v"(rhs_10_6));

  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_const_lhs)
      : [mask]  "Yk"((__mmask16) low12),
        [a]     "v"(lhs_10),
        [b]     "v"(lhs_4));
  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_const_rhs)
      : [mask]  "Yk"((__mmask16) low12),
        [a]     "v"(rhs_10),
        [b]     "v"(rhs_4));
#else
  blend_var_lhs = _mm512_mask_blend_epi32(low12, lhs_4_9, lhs_10_6);
  blend_var_rhs = _mm512_mask_blend_epi32(low12, rhs_4_9, rhs_10_6);

  blend_const_lhs = _mm512_mask_blend_epi32(low12, lhs_10, lhs_4);
  blend_const_rhs = _mm512_mask_blend_epi32(low12, rhs_10, rhs_4);
#endif

//   print8x64(blend_var_lhs, "blend_var_lhs");
//   print8x64(blend_var_rhs, "blend_var_rhs");
//   print8x64(blend_const_lhs, "blend_const_lhs");
//   print8x64(blend_const_rhs, "blend_const_rhs");

  sublhs_512 = _mm512_sub_epi32(blend_const_lhs, blend_var_lhs);
  subrhs_512 = _mm512_sub_epi32(blend_var_rhs, blend_const_rhs);

//   print8x64(sublhs_512, "sublhs_512");
//   print8x64(subrhs_512, "subrhs_512");

  mul_512 = _mm512_mul_epi32(sublhs_512, subrhs_512);
//   print8x64(mul_512, "mul_512");
  accum3 = _mm512_add_epi64(accum3, mul_512);
//   print8x64(accum3, "accum3");

  // 10 * [0 1 2 3] Done here because we have loaded everything already.
  sublhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(lhs_10),
                                _mm512_castsi512_si256(lhs_0_7));
  subrhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(rhs_0_7),
                                _mm512_castsi512_si256(rhs_10));
//   print4x64(sublhs_256, "sublhs_256");
//   print4x64(subrhs_256, "subrhs_256");

  mul_256 = _mm256_mul_epi32(sublhs_256, subrhs_256);
//   print4x64(mul_256, "mul_256");
  accum10 = _mm256_add_epi64(accum10, mul_256);
//   print4x64(accum10, "accum10");

  // 0 * [3 4 5 6 7 8 9 10]
  __m512i lhs_3_10 = _mm512_loadu_si512((__m512i *) &x->limbs[4]);
  __m512i rhs_3_10 = _mm512_loadu_si512((__m512i *) &y->limbs[4]);

//   print8x64(lhs_3_10, "lhs_3_10");
//   print8x64(rhs_3_10, "rhs_3_10");

  __m512i lhs_0 = _mm512_set1_epi32(x->limbs[1]);
  __m512i rhs_0 = _mm512_set1_epi32(y->limbs[1]);

//   print8x64(lhs_0, "lhs_0");
//   print8x64(rhs_0, "rhs_0");

  sublhs_512 = _mm512_sub_epi32(lhs_0, lhs_3_10);
  subrhs_512 = _mm512_sub_epi32(rhs_3_10, rhs_0);
//   print8x64(sublhs_512, "sublhs_512");
//   print8x64(subrhs_512, "subrhs_512");

  mul_512 = _mm512_mul_epi32(sublhs_512, subrhs_512);
//   print8x64(mul_512, "mul_512");
  accum3 = _mm512_add_epi64(accum3, mul_512);
//   print8x64(accum3, "accum3");

  // 6 * [4 5] and 0 * [1 2]
  __m256i lhs_6 = _mm256_set1_epi32(x->limbs[7]);
  __m256i rhs_6 = _mm256_set1_epi32(y->limbs[7]);

//   print4x64(lhs_6, "lhs_6");
//   print4x64(rhs_6, "rhs_6");

  blend_var_lhs_256 = _mm256_blend_epi32(_mm512_castsi512_si256(lhs_10_6),
                                         _mm512_castsi512_si256(lhs_4_9),
                                         0x0f);
  blend_var_rhs_256 = _mm256_blend_epi32(_mm512_castsi512_si256(rhs_10_6),
                                         _mm512_castsi512_si256(rhs_4_9),
                                         0x0f);

  blend_const_lhs_256 = _mm256_blend_epi32(_mm512_castsi512_si256(lhs_0),
                                           lhs_6,
                                           0x0f);
  blend_const_rhs_256 = _mm256_blend_epi32(_mm512_castsi512_si256(rhs_0),
                                           rhs_6,
                                           0x0f);
//   print4x64(blend_var_lhs_256, "blend_var_lhs_256");
//   print4x64(blend_var_rhs_256, "blend_var_rhs_256");
//   print4x64(blend_const_lhs_256, "blend_const_lhs_256");
//   print4x64(blend_const_rhs_256, "blend_const_rhs_256");

  sublhs_256 = _mm256_sub_epi32(blend_const_lhs_256, blend_var_lhs_256);
  subrhs_256 = _mm256_sub_epi32(blend_var_rhs_256, blend_const_rhs_256);
//   print4x64(sublhs_256, "sublhs_256");
//   print4x64(subrhs_256, "subrhs_256");

  mul_256 = _mm256_mul_epi32(sublhs_256, subrhs_256);
//   print4x64(mul_256, "mul_256");
  accum10 = _mm256_add_epi64(accum10, mul_256);
//   print4x64(accum10, "accum10");

  // 7 * [3 4 5 6]
  __m256i lhs_7 = _mm256_set1_epi32(x->limbs[8]);
  __m256i rhs_7 = _mm256_set1_epi32(y->limbs[8]);

//   print4x64(lhs_7, "lhs_7");
//   print4x64(rhs_7, "rhs_7");

  sublhs_256 = _mm256_sub_epi32(lhs_7, _mm512_castsi512_si256(lhs_3_10));
  subrhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(rhs_3_10), rhs_7);
//   print4x64(sublhs_256, "sublhs_256");
//   print4x64(subrhs_256, "subrhs_256");
  mul_256 = _mm256_mul_epi32(sublhs_256, subrhs_256);
//   print4x64(mul_256, "mul_256");
  accum10 = _mm256_add_epi64(accum10, mul_256);
//   print4x64(accum10, "accum10");

  // 1 * [2 3 4 5 6 7 8 9]
  __m512i lhs_2_9 = _mm512_loadu_si512((__m512i *) &x->limbs[3]);
  __m512i rhs_2_9 = _mm512_loadu_si512((__m512i *) &y->limbs[3]);

//   print8x64(lhs_2_9, "lhs_2_9");
//   print8x64(rhs_2_9, "rhs_2_9");

  __m512i lhs_1 = _mm512_set1_epi32(x->limbs[2]);
  __m512i rhs_1 = _mm512_set1_epi32(y->limbs[2]);

//   print8x64(lhs_1, "lhs_1");
//   print8x64(rhs_1, "rhs_1");

  sublhs_512 = _mm512_sub_epi32(lhs_1, lhs_2_9);
  subrhs_512 = _mm512_sub_epi32(rhs_2_9, rhs_1);
//   print8x64(sublhs_512, "sublhs_512");
//   print8x64(subrhs_512, "subrhs_512");
  mul_512 = _mm512_mul_epi32(sublhs_512, subrhs_512);
//   print8x64(mul_512, "mul_512");
  accum3 = _mm512_add_epi64(accum3, mul_512);
//   print8x64(accum3, "accum3");

  _mm512_storeu_si512((__m512i*) &temp.limbs[4], accum3);

  // 8 * [2 3 4 5]
  sublhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(lhs_8),
                                _mm512_castsi512_si256(lhs_2_9));
  subrhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(rhs_2_9),
                                _mm512_castsi512_si256(rhs_8));
//   print4x64(sublhs_256, "sublhs_256");
//   print4x64(subrhs_256, "subrhs_256");
  mul_256 = _mm256_mul_epi32(sublhs_256, subrhs_256);
//   print4x64(mul_256, "mul_256");
  accum10 = _mm256_add_epi64(accum10, mul_256);
//   print4x64(accum10, "accum10");
  _mm256_store_si256((__m256i*) &temp.limbs[0], accum10);

  reduce_step_wide(&temp, &temp);
  reduce_step_wide(result, &temp);
}

void mul_wide_narrow_slow(
  residue_wide_t *result, const residue_wide_t *x, const residue_narrow_t *y) {
  residue_wide_t y_temp;
  widen(&y_temp, y);
  mul_wide(result, x, &y_temp);
}
// Multiply a wide residues by a narrow and produce a wide result. The result is
// reduced to 32 bits, but not narrowed for performance reasons.
void mul_wide_narrow(
  residue_wide_t *result, const residue_wide_t *x, const residue_narrow_t *y) {

  residue_wide_t temp;

  // Two accumulators
  __m256i accum10;
  __m512i accum3;

  // Temporaries for sub sub mul
  __m512i sublhs_512, subrhs_512, mul_512;
  __m256i sublhs_256, subrhs_256, mul_256;

  __m512i blend_var_lhs, blend_var_rhs;
  __m512i blend_const_lhs, blend_const_rhs;

  __m256i blend_var_lhs_256, blend_var_rhs_256;
  __m256i blend_const_lhs_256, blend_const_rhs_256;

  __mmask16 low4, low8, low12;

  low4 = _cvtu32_mask16(0xf);
  low8 = _cvtu32_mask16(0xff);
  low12 = _cvtu32_mask16(0xfff);

  // 8 * [6 7] and 2 * [3 4 5 6 7 8]
  __m512i lhs_6_7 = _mm512_castsi128_si512(_mm_loadu_si128((__m128i *) &x->limbs[7]));
  __m512i rhs_6_7 = _mm512_castsi256_si512(loadu_extend_32_64((__m128i *) &y->limbs[7]));
//   print8x64(lhs_6_7, "lhs_6_7");
//   print8x64(rhs_6_7, "rhs_6_7");

  __m512i lhs_1_8 = _mm512_loadu_si512((__m512i *) &x->limbs[2]);
  __m512i rhs_1_8 = loadu512_extend_32_64((__m256i *) &y->limbs[2]);

//   print8x64(lhs_1_8, "lhs_1_8");
//   print8x64(rhs_1_8, "rhs_1_8");

  __m512i lhs_8 = _mm512_set1_epi32(x->limbs[9]);
  __m512i rhs_8 = _mm512_set1_epi32(y->limbs[9]);

//   print8x64(lhs_8, "lhs_8");
//   print8x64(rhs_8, "rhs_8");

  __m512i lhs_2 = _mm512_set1_epi32(x->limbs[3]);
  __m512i rhs_2 = _mm512_set1_epi32(y->limbs[3]);

//   print8x64(lhs_2, "lhs_2");
//   print8x64(rhs_2, "rhs_2");

#if FORCE_BLEND
  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_var_lhs)
      : [mask]  "Yk"((__mmask16) low4),
        [a]     "v"(lhs_6_7),
        [b]     "v"(lhs_1_8));
  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_var_rhs)
      : [mask]  "Yk"((__mmask16) low4),
        [a]     "v"(rhs_6_7),
        [b]     "v"(rhs_1_8));

  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_const_lhs)
      : [mask]  "Yk"((__mmask16) low4),
        [a]     "v"(lhs_8),
        [b]     "v"(lhs_2));
  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_const_rhs)
      : [mask]  "Yk"((__mmask16) low4),
        [a]     "v"(rhs_8),
        [b]     "v"(rhs_2));
#else
  blend_var_lhs = _mm512_mask_blend_epi32(low4, lhs_6_7, lhs_1_8);
  blend_var_rhs = _mm512_mask_blend_epi32(low4, rhs_6_7, rhs_1_8);

  blend_const_lhs = _mm512_mask_blend_epi32(low4, lhs_8, lhs_2);
  blend_const_rhs = _mm512_mask_blend_epi32(low4, rhs_8, rhs_2);
#endif

//   print8x64(blend_var_lhs, "blend_var_lhs");
//   print8x64(blend_var_rhs, "blend_var_rhs");
//   print8x64(blend_const_lhs, "blend_const_lhs");
//   print8x64(blend_const_rhs, "blend_const_rhs");

  sublhs_512 = _mm512_sub_epi32(blend_const_lhs, blend_var_lhs);
  subrhs_512 = _mm512_sub_epi32(blend_var_rhs, blend_const_rhs);
//   print8x64(sublhs_512, "sublhs_512");
//   print8x64(subrhs_512, "subrhs_512");

  accum3 = _mm512_mul_epi32(sublhs_512, subrhs_512);
//   print8x64(accum3, "accum3");

  // 9 * [5 6 7 8] and 3 * [4 5 6 7]
  __m512i lhs_5_8 = _mm512_castsi256_si512(_mm256_loadu_si256((__m256i *) &x->limbs[6]));
  __m512i rhs_5_8 = _mm512_castsi256_si512(loadu_extend_32_64((__m128i *) &y->limbs[6]));

//   print8x64(lhs_5_8, "lhs_5_8");
//   print8x64(rhs_5_8, "rhs_5_8");

  __m512i lhs_0_7 = _mm512_loadu_si512((__m512i *) &x->limbs[1]);
  __m512i rhs_0_7 = loadu512_extend_32_64((__m256i *) &y->limbs[1]);

//   print8x64(lhs_0_7, "lhs_0_7");
//   print8x64(rhs_0_7, "rhs_0_7");

  __m512i lhs_9 = _mm512_set1_epi32(x->limbs[10]);
  __m512i rhs_9 = _mm512_set1_epi32(y->limbs[10]);

//   print8x64(lhs_9, "lhs_9");
//   print8x64(rhs_9, "rhs_9");

  __m512i lhs_3 = _mm512_set1_epi32(x->limbs[4]);
  __m512i rhs_3 = _mm512_set1_epi32(y->limbs[4]);

//   print8x64(lhs_3, "lhs_3");
//   print8x64(rhs_3, "rhs_3");

#if FORCE_BLEND
  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_var_lhs)
      : [mask]  "Yk"((__mmask16) low8),
        [a]     "v"(lhs_5_8),
        [b]     "v"(lhs_0_7));
  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_var_rhs)
      : [mask]  "Yk"((__mmask16) low8),
        [a]     "v"(rhs_5_8),
        [b]     "v"(rhs_0_7));

  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_const_lhs)
      : [mask]  "Yk"((__mmask16) low8),
        [a]     "v"(lhs_9),
        [b]     "v"(lhs_3));
  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_const_rhs)
      : [mask]  "Yk"((__mmask16) low8),
        [a]     "v"(rhs_9),
        [b]     "v"(rhs_3));
#else
  blend_var_lhs = _mm512_mask_blend_epi32(low8, lhs_5_8, lhs_0_7);
  blend_var_rhs = _mm512_mask_blend_epi32(low8, rhs_5_8, rhs_0_7);

  blend_const_lhs = _mm512_mask_blend_epi32(low8, lhs_9, lhs_3);
  blend_const_rhs = _mm512_mask_blend_epi32(low8, rhs_9, rhs_3);
#endif

//   print8x64(blend_var_lhs, "blend_var_lhs");
//   print8x64(blend_var_rhs, "blend_var_rhs");
//   print8x64(blend_const_lhs, "blend_const_lhs");
//   print8x64(blend_const_rhs, "blend_const_rhs");

  sublhs_512 = _mm512_sub_epi32(blend_const_lhs, blend_var_lhs);
  subrhs_512 = _mm512_sub_epi32(blend_var_rhs, blend_const_rhs);
//   print8x64(sublhs_512, "sublhs_512");
//   print8x64(subrhs_512, "subrhs_512");

  mul_512 = _mm512_mul_epi32(sublhs_512, subrhs_512);
//   print8x64(mul_512, "mul_512");
  accum3 = _mm512_add_epi64(accum3, mul_512);
//   print8x64(accum3, "accum3");

  // 9 * [1 2 3 4] Done here because we have loaded everything already.
  sublhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(lhs_9),
                                _mm512_castsi512_si256(lhs_1_8));
  subrhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(rhs_1_8),
                                _mm512_castsi512_si256(rhs_9));
//   print4x64(sublhs_256, "sublhs_256");
//   print4x64(subrhs_256, "subrhs_256");

  accum10 = _mm256_mul_epi32(sublhs_256, subrhs_256);
//   print4x64(accum10, "accum10");

  // 10 * [4 5 6 7 8 9] and 4 * [5 6]
  __m512i lhs_4_9 = _mm512_mask_loadu_epi32(_mm512_setzero(), low12, &x->limbs[5]);
  __m512i rhs_4_9 = loadu512_mask_extend_32_64((__m256i *) &y->limbs[5], _cvtu32_mask8(0x3f));

//   print8x64(lhs_4_9, "lhs_4_9");
//   print8x64(rhs_4_9, "rhs_4_9");

  __m512i lhs_10_6 = _mm512_loadu_si512((__m512i *) &x->limbs[0]);
  __m512i rhs_10_6 = loadu512_extend_32_64((__m256i *) &y->limbs[0]);

//   print8x64(lhs_10_6, "lhs_10_6");
//   print8x64(rhs_10_6, "rhs_10_6");

  __m512i lhs_10 = _mm512_set1_epi32(x->limbs[0]);
  __m512i rhs_10 = _mm512_set1_epi32(y->limbs[0]);

//   print8x64(lhs_10, "lhs_10");
//   print8x64(rhs_10, "rhs_10");

  __m512i lhs_4 = _mm512_set1_epi32(x->limbs[5]);
  __m512i rhs_4 = _mm512_set1_epi32(y->limbs[5]);

//   print8x64(lhs_4, "lhs_4");
//   print8x64(rhs_4, "rhs_4");

#if FORCE_BLEND
  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_var_lhs)
      : [mask]  "Yk"((__mmask16) low12),
        [a]     "v"(lhs_4_9),
        [b]     "v"(lhs_10_6));
  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_var_rhs)
      : [mask]  "Yk"((__mmask16) low12),
        [a]     "v"(rhs_4_9),
        [b]     "v"(rhs_10_6));

  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_const_lhs)
      : [mask]  "Yk"((__mmask16) low12),
        [a]     "v"(lhs_10),
        [b]     "v"(lhs_4));
  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_const_rhs)
      : [mask]  "Yk"((__mmask16) low12),
        [a]     "v"(rhs_10),
        [b]     "v"(rhs_4));
#else
  blend_var_lhs = _mm512_mask_blend_epi32(low12, lhs_4_9, lhs_10_6);
  blend_var_rhs = _mm512_mask_blend_epi32(low12, rhs_4_9, rhs_10_6);

  blend_const_lhs = _mm512_mask_blend_epi32(low12, lhs_10, lhs_4);
  blend_const_rhs = _mm512_mask_blend_epi32(low12, rhs_10, rhs_4);
#endif

//   print8x64(blend_var_lhs, "blend_var_lhs");
//   print8x64(blend_var_rhs, "blend_var_rhs");
//   print8x64(blend_const_lhs, "blend_const_lhs");
//   print8x64(blend_const_rhs, "blend_const_rhs");

  sublhs_512 = _mm512_sub_epi32(blend_const_lhs, blend_var_lhs);
  subrhs_512 = _mm512_sub_epi32(blend_var_rhs, blend_const_rhs);

//   print8x64(sublhs_512, "sublhs_512");
//   print8x64(subrhs_512, "subrhs_512");

  mul_512 = _mm512_mul_epi32(sublhs_512, subrhs_512);
//   print8x64(mul_512, "mul_512");
  accum3 = _mm512_add_epi64(accum3, mul_512);
//   print8x64(accum3, "accum3");

  // 10 * [0 1 2 3] Done here because we have loaded everything already.
  sublhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(lhs_10),
                                _mm512_castsi512_si256(lhs_0_7));
  subrhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(rhs_0_7),
                                _mm512_castsi512_si256(rhs_10));
//   print4x64(sublhs_256, "sublhs_256");
//   print4x64(subrhs_256, "subrhs_256");

  mul_256 = _mm256_mul_epi32(sublhs_256, subrhs_256);
//   print4x64(mul_256, "mul_256");
  accum10 = _mm256_add_epi64(accum10, mul_256);
//   print4x64(accum10, "accum10");

  // 0 * [3 4 5 6 7 8 9 10]
  __m512i lhs_3_10 = _mm512_loadu_si512((__m512i *) &x->limbs[4]);
  __m512i rhs_3_10 = loadu512_extend_32_64((__m256i *) &y->limbs[4]);

//   print8x64(lhs_3_10, "lhs_3_10");
//   print8x64(rhs_3_10, "rhs_3_10");

  __m512i lhs_0 = _mm512_set1_epi32(x->limbs[1]);
  __m512i rhs_0 = _mm512_set1_epi32(y->limbs[1]);

//   print8x64(lhs_0, "lhs_0");
//   print8x64(rhs_0, "rhs_0");

  sublhs_512 = _mm512_sub_epi32(lhs_0, lhs_3_10);
  subrhs_512 = _mm512_sub_epi32(rhs_3_10, rhs_0);
//   print8x64(sublhs_512, "sublhs_512");
//   print8x64(subrhs_512, "subrhs_512");

  mul_512 = _mm512_mul_epi32(sublhs_512, subrhs_512);
//   print8x64(mul_512, "mul_512");
  accum3 = _mm512_add_epi64(accum3, mul_512);
//   print8x64(accum3, "accum3");

  // 6 * [4 5] and 0 * [1 2]
  __m256i lhs_6 = _mm256_set1_epi32(x->limbs[7]);
  __m256i rhs_6 = _mm256_set1_epi32(y->limbs[7]);

//   print4x64(lhs_6, "lhs_6");
//   print4x64(rhs_6, "rhs_6");

  blend_var_lhs_256 = _mm256_blend_epi32(_mm512_castsi512_si256(lhs_10_6),
                                         _mm512_castsi512_si256(lhs_4_9),
                                         0x0f);
  blend_var_rhs_256 = _mm256_blend_epi32(_mm512_castsi512_si256(rhs_10_6),
                                         _mm512_castsi512_si256(rhs_4_9),
                                         0x0f);

  blend_const_lhs_256 = _mm256_blend_epi32(_mm512_castsi512_si256(lhs_0),
                                           lhs_6,
                                           0x0f);
  blend_const_rhs_256 = _mm256_blend_epi32(_mm512_castsi512_si256(rhs_0),
                                           rhs_6,
                                           0x0f);
//   print4x64(blend_var_lhs_256, "blend_var_lhs_256");
//   print4x64(blend_var_rhs_256, "blend_var_rhs_256");
//   print4x64(blend_const_lhs_256, "blend_const_lhs_256");
//   print4x64(blend_const_rhs_256, "blend_const_rhs_256");

  sublhs_256 = _mm256_sub_epi32(blend_const_lhs_256, blend_var_lhs_256);
  subrhs_256 = _mm256_sub_epi32(blend_var_rhs_256, blend_const_rhs_256);
//   print4x64(sublhs_256, "sublhs_256");
//   print4x64(subrhs_256, "subrhs_256");

  mul_256 = _mm256_mul_epi32(sublhs_256, subrhs_256);
//   print4x64(mul_256, "mul_256");
  accum10 = _mm256_add_epi64(accum10, mul_256);
//   print4x64(accum10, "accum10");

  // 7 * [3 4 5 6]
  __m256i lhs_7 = _mm256_set1_epi32(x->limbs[8]);
  __m256i rhs_7 = _mm256_set1_epi32(y->limbs[8]);

//   print4x64(lhs_7, "lhs_7");
//   print4x64(rhs_7, "rhs_7");

  sublhs_256 = _mm256_sub_epi32(lhs_7, _mm512_castsi512_si256(lhs_3_10));
  subrhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(rhs_3_10), rhs_7);
//   print4x64(sublhs_256, "sublhs_256");
//   print4x64(subrhs_256, "subrhs_256");
  mul_256 = _mm256_mul_epi32(sublhs_256, subrhs_256);
//   print4x64(mul_256, "mul_256");
  accum10 = _mm256_add_epi64(accum10, mul_256);
//   print4x64(accum10, "accum10");

  // 1 * [2 3 4 5 6 7 8 9]
  __m512i lhs_2_9 = _mm512_loadu_si512((__m512i *) &x->limbs[3]);
  __m512i rhs_2_9 = loadu512_extend_32_64((__m256i *) &y->limbs[3]);

//   print8x64(lhs_2_9, "lhs_2_9");
//   print8x64(rhs_2_9, "rhs_2_9");

  __m512i lhs_1 = _mm512_set1_epi32(x->limbs[2]);
  __m512i rhs_1 = _mm512_set1_epi32(y->limbs[2]);

//   print8x64(lhs_1, "lhs_1");
//   print8x64(rhs_1, "rhs_1");

  sublhs_512 = _mm512_sub_epi32(lhs_1, lhs_2_9);
  subrhs_512 = _mm512_sub_epi32(rhs_2_9, rhs_1);
//   print8x64(sublhs_512, "sublhs_512");
//   print8x64(subrhs_512, "subrhs_512");
  mul_512 = _mm512_mul_epi32(sublhs_512, subrhs_512);
//   print8x64(mul_512, "mul_512");
  accum3 = _mm512_add_epi64(accum3, mul_512);
//   print8x64(accum3, "accum3");

  _mm512_storeu_si512((__m512i*) &temp.limbs[4], accum3);

  // 8 * [2 3 4 5]
  sublhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(lhs_8),
                                _mm512_castsi512_si256(lhs_2_9));
  subrhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(rhs_2_9),
                                _mm512_castsi512_si256(rhs_8));
//   print4x64(sublhs_256, "sublhs_256");
//   print4x64(subrhs_256, "subrhs_256");
  mul_256 = _mm256_mul_epi32(sublhs_256, subrhs_256);
//   print4x64(mul_256, "mul_256");
  accum10 = _mm256_add_epi64(accum10, mul_256);
//   print4x64(accum10, "accum10");
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

  // Two accumulators
  __m256i accum10;
  __m512i accum3;

  // Temporaries for sub sub mul
  __m512i sublhs_512, subrhs_512, mul_512;
  __m256i sublhs_256, subrhs_256, mul_256;

  __m512i blend_var_lhs, blend_var_rhs;
  __m512i blend_const_lhs, blend_const_rhs;

  __m256i blend_var_lhs_256, blend_var_rhs_256;
  __m256i blend_const_lhs_256, blend_const_rhs_256;

  __mmask16 low4, low8, low12;

  low4 = _cvtu32_mask16(0xf);
  low8 = _cvtu32_mask16(0xff);
  low12 = _cvtu32_mask16(0xfff);

  // 8 * [6 7] and 2 * [3 4 5 6 7 8]
  __m512i lhs_6_7 = _mm512_castsi256_si512(loadu_extend_32_64((__m128i *) &x->limbs[7]));
  __m512i rhs_6_7 = _mm512_castsi256_si512(loadu_extend_32_64((__m128i *) &y->limbs[7]));
//   print8x64(lhs_6_7, "lhs_6_7");
//   print8x64(rhs_6_7, "rhs_6_7");

  __m512i lhs_1_8 = loadu512_extend_32_64((__m256i *) &x->limbs[2]);
  __m512i rhs_1_8 = loadu512_extend_32_64((__m256i *) &y->limbs[2]);

//   print8x64(lhs_1_8, "lhs_1_8");
//   print8x64(rhs_1_8, "rhs_1_8");

  __m512i lhs_8 = _mm512_set1_epi32(x->limbs[9]);
  __m512i rhs_8 = _mm512_set1_epi32(y->limbs[9]);

//   print8x64(lhs_8, "lhs_8");
//   print8x64(rhs_8, "rhs_8");

  __m512i lhs_2 = _mm512_set1_epi32(x->limbs[3]);
  __m512i rhs_2 = _mm512_set1_epi32(y->limbs[3]);

//   print8x64(lhs_2, "lhs_2");
//   print8x64(rhs_2, "rhs_2");

#if FORCE_BLEND
  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_var_lhs)
      : [mask]  "Yk"((__mmask16) low4),
        [a]     "v"(lhs_6_7),
        [b]     "v"(lhs_1_8));
  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_var_rhs)
      : [mask]  "Yk"((__mmask16) low4),
        [a]     "v"(rhs_6_7),
        [b]     "v"(rhs_1_8));

  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_const_lhs)
      : [mask]  "Yk"((__mmask16) low4),
        [a]     "v"(lhs_8),
        [b]     "v"(lhs_2));
  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_const_rhs)
      : [mask]  "Yk"((__mmask16) low4),
        [a]     "v"(rhs_8),
        [b]     "v"(rhs_2));
#else
  blend_var_lhs = _mm512_mask_blend_epi32(low4, lhs_6_7, lhs_1_8);
  blend_var_rhs = _mm512_mask_blend_epi32(low4, rhs_6_7, rhs_1_8);

  blend_const_lhs = _mm512_mask_blend_epi32(low4, lhs_8, lhs_2);
  blend_const_rhs = _mm512_mask_blend_epi32(low4, rhs_8, rhs_2);
#endif

//   print8x64(blend_var_lhs, "blend_var_lhs");
//   print8x64(blend_var_rhs, "blend_var_rhs");
//   print8x64(blend_const_lhs, "blend_const_lhs");
//   print8x64(blend_const_rhs, "blend_const_rhs");

  sublhs_512 = _mm512_sub_epi32(blend_const_lhs, blend_var_lhs);
  subrhs_512 = _mm512_sub_epi32(blend_var_rhs, blend_const_rhs);
//   print8x64(sublhs_512, "sublhs_512");
//   print8x64(subrhs_512, "subrhs_512");

  accum3 = _mm512_mul_epi32(sublhs_512, subrhs_512);
//   print8x64(accum3, "accum3");

  // 9 * [5 6 7 8] and 3 * [4 5 6 7]
  __m512i lhs_5_8 = _mm512_castsi256_si512(loadu_extend_32_64((__m128i *) &x->limbs[6]));
  __m512i rhs_5_8 = _mm512_castsi256_si512(loadu_extend_32_64((__m128i *) &y->limbs[6]));

//   print8x64(lhs_5_8, "lhs_5_8");
//   print8x64(rhs_5_8, "rhs_5_8");

  __m512i lhs_0_7 = loadu512_extend_32_64((__m256i *) &x->limbs[1]);
  __m512i rhs_0_7 = loadu512_extend_32_64((__m256i *) &y->limbs[1]);

//   print8x64(lhs_0_7, "lhs_0_7");
//   print8x64(rhs_0_7, "rhs_0_7");

  __m512i lhs_9 = _mm512_set1_epi32(x->limbs[10]);
  __m512i rhs_9 = _mm512_set1_epi32(y->limbs[10]);

//   print8x64(lhs_9, "lhs_9");
//   print8x64(rhs_9, "rhs_9");

  __m512i lhs_3 = _mm512_set1_epi32(x->limbs[4]);
  __m512i rhs_3 = _mm512_set1_epi32(y->limbs[4]);

//   print8x64(lhs_3, "lhs_3");
//   print8x64(rhs_3, "rhs_3");

#if FORCE_BLEND
  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_var_lhs)
      : [mask]  "Yk"((__mmask16) low8),
        [a]     "v"(lhs_5_8),
        [b]     "v"(lhs_0_7));
  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_var_rhs)
      : [mask]  "Yk"((__mmask16) low8),
        [a]     "v"(rhs_5_8),
        [b]     "v"(rhs_0_7));

  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_const_lhs)
      : [mask]  "Yk"((__mmask16) low8),
        [a]     "v"(lhs_9),
        [b]     "v"(lhs_3));
  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_const_rhs)
      : [mask]  "Yk"((__mmask16) low8),
        [a]     "v"(rhs_9),
        [b]     "v"(rhs_3));
#else
  blend_var_lhs = _mm512_mask_blend_epi32(low8, lhs_5_8, lhs_0_7);
  blend_var_rhs = _mm512_mask_blend_epi32(low8, rhs_5_8, rhs_0_7);

  blend_const_lhs = _mm512_mask_blend_epi32(low8, lhs_9, lhs_3);
  blend_const_rhs = _mm512_mask_blend_epi32(low8, rhs_9, rhs_3);
#endif

//   print8x64(blend_var_lhs, "blend_var_lhs");
//   print8x64(blend_var_rhs, "blend_var_rhs");
//   print8x64(blend_const_lhs, "blend_const_lhs");
//   print8x64(blend_const_rhs, "blend_const_rhs");

  sublhs_512 = _mm512_sub_epi32(blend_const_lhs, blend_var_lhs);
  subrhs_512 = _mm512_sub_epi32(blend_var_rhs, blend_const_rhs);
//   print8x64(sublhs_512, "sublhs_512");
//   print8x64(subrhs_512, "subrhs_512");

  mul_512 = _mm512_mul_epi32(sublhs_512, subrhs_512);
//   print8x64(mul_512, "mul_512");
  accum3 = _mm512_add_epi64(accum3, mul_512);
//   print8x64(accum3, "accum3");

  // 9 * [1 2 3 4] Done here because we have loaded everything already.
  sublhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(lhs_9),
                                _mm512_castsi512_si256(lhs_1_8));
  subrhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(rhs_1_8),
                                _mm512_castsi512_si256(rhs_9));
//   print4x64(sublhs_256, "sublhs_256");
//   print4x64(subrhs_256, "subrhs_256");

  accum10 = _mm256_mul_epi32(sublhs_256, subrhs_256);
//   print4x64(accum10, "accum10");

  // 10 * [4 5 6 7 8 9] and 4 * [5 6]
  __m512i lhs_4_9 = loadu512_mask_extend_32_64((__m256i *) &x->limbs[5], _cvtu32_mask8(0x3f));
  __m512i rhs_4_9 = loadu512_mask_extend_32_64((__m256i *) &y->limbs[5], _cvtu32_mask8(0x3f));

//   print8x64(lhs_4_9, "lhs_4_9");
//   print8x64(rhs_4_9, "rhs_4_9");

  __m512i lhs_10_6 = loadu512_extend_32_64((__m256i *) &x->limbs[0]);
  __m512i rhs_10_6 = loadu512_extend_32_64((__m256i *) &y->limbs[0]);

//   print8x64(lhs_10_6, "lhs_10_6");
//   print8x64(rhs_10_6, "rhs_10_6");

  __m512i lhs_10 = _mm512_set1_epi32(x->limbs[0]);
  __m512i rhs_10 = _mm512_set1_epi32(y->limbs[0]);

//   print8x64(lhs_10, "lhs_10");
//   print8x64(rhs_10, "rhs_10");

  __m512i lhs_4 = _mm512_set1_epi32(x->limbs[5]);
  __m512i rhs_4 = _mm512_set1_epi32(y->limbs[5]);

//   print8x64(lhs_4, "lhs_4");
//   print8x64(rhs_4, "rhs_4");

#if FORCE_BLEND
  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_var_lhs)
      : [mask]  "Yk"((__mmask16) low12),
        [a]     "v"(lhs_4_9),
        [b]     "v"(lhs_10_6));
  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_var_rhs)
      : [mask]  "Yk"((__mmask16) low12),
        [a]     "v"(rhs_4_9),
        [b]     "v"(rhs_10_6));

  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_const_lhs)
      : [mask]  "Yk"((__mmask16) low12),
        [a]     "v"(lhs_10),
        [b]     "v"(lhs_4));
  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_const_rhs)
      : [mask]  "Yk"((__mmask16) low12),
        [a]     "v"(rhs_10),
        [b]     "v"(rhs_4));
#else
  blend_var_lhs = _mm512_mask_blend_epi32(low12, lhs_4_9, lhs_10_6);
  blend_var_rhs = _mm512_mask_blend_epi32(low12, rhs_4_9, rhs_10_6);

  blend_const_lhs = _mm512_mask_blend_epi32(low12, lhs_10, lhs_4);
  blend_const_rhs = _mm512_mask_blend_epi32(low12, rhs_10, rhs_4);
#endif

//   print8x64(blend_var_lhs, "blend_var_lhs");
//   print8x64(blend_var_rhs, "blend_var_rhs");
//   print8x64(blend_const_lhs, "blend_const_lhs");
//   print8x64(blend_const_rhs, "blend_const_rhs");

  sublhs_512 = _mm512_sub_epi32(blend_const_lhs, blend_var_lhs);
  subrhs_512 = _mm512_sub_epi32(blend_var_rhs, blend_const_rhs);

//   print8x64(sublhs_512, "sublhs_512");
//   print8x64(subrhs_512, "subrhs_512");

  mul_512 = _mm512_mul_epi32(sublhs_512, subrhs_512);
//   print8x64(mul_512, "mul_512");
  accum3 = _mm512_add_epi64(accum3, mul_512);
//   print8x64(accum3, "accum3");

  // 10 * [0 1 2 3] Done here because we have loaded everything already.
  sublhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(lhs_10),
                                _mm512_castsi512_si256(lhs_0_7));
  subrhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(rhs_0_7),
                                _mm512_castsi512_si256(rhs_10));
//   print4x64(sublhs_256, "sublhs_256");
//   print4x64(subrhs_256, "subrhs_256");

  mul_256 = _mm256_mul_epi32(sublhs_256, subrhs_256);
//   print4x64(mul_256, "mul_256");
  accum10 = _mm256_add_epi64(accum10, mul_256);
//   print4x64(accum10, "accum10");

  // 0 * [3 4 5 6 7 8 9 10]
  __m512i lhs_3_10 = loadu512_extend_32_64((__m256i *) &x->limbs[4]);
  __m512i rhs_3_10 = loadu512_extend_32_64((__m256i *) &y->limbs[4]);

//   print8x64(lhs_3_10, "lhs_3_10");
//   print8x64(rhs_3_10, "rhs_3_10");

  __m512i lhs_0 = _mm512_set1_epi32(x->limbs[1]);
  __m512i rhs_0 = _mm512_set1_epi32(y->limbs[1]);

//   print8x64(lhs_0, "lhs_0");
//   print8x64(rhs_0, "rhs_0");

  sublhs_512 = _mm512_sub_epi32(lhs_0, lhs_3_10);
  subrhs_512 = _mm512_sub_epi32(rhs_3_10, rhs_0);
//   print8x64(sublhs_512, "sublhs_512");
//   print8x64(subrhs_512, "subrhs_512");

  mul_512 = _mm512_mul_epi32(sublhs_512, subrhs_512);
//   print8x64(mul_512, "mul_512");
  accum3 = _mm512_add_epi64(accum3, mul_512);
//   print8x64(accum3, "accum3");

  // 6 * [4 5] and 0 * [1 2]
  __m256i lhs_6 = _mm256_set1_epi32(x->limbs[7]);
  __m256i rhs_6 = _mm256_set1_epi32(y->limbs[7]);

//   print4x64(lhs_6, "lhs_6");
//   print4x64(rhs_6, "rhs_6");

  blend_var_lhs_256 = _mm256_blend_epi32(_mm512_castsi512_si256(lhs_10_6),
                                         _mm512_castsi512_si256(lhs_4_9),
                                         0x0f);
  blend_var_rhs_256 = _mm256_blend_epi32(_mm512_castsi512_si256(rhs_10_6),
                                         _mm512_castsi512_si256(rhs_4_9),
                                         0x0f);

  blend_const_lhs_256 = _mm256_blend_epi32(_mm512_castsi512_si256(lhs_0),
                                           lhs_6,
                                           0x0f);
  blend_const_rhs_256 = _mm256_blend_epi32(_mm512_castsi512_si256(rhs_0),
                                           rhs_6,
                                           0x0f);
//   print4x64(blend_var_lhs_256, "blend_var_lhs_256");
//   print4x64(blend_var_rhs_256, "blend_var_rhs_256");
//   print4x64(blend_const_lhs_256, "blend_const_lhs_256");
//   print4x64(blend_const_rhs_256, "blend_const_rhs_256");

  sublhs_256 = _mm256_sub_epi32(blend_const_lhs_256, blend_var_lhs_256);
  subrhs_256 = _mm256_sub_epi32(blend_var_rhs_256, blend_const_rhs_256);
//   print4x64(sublhs_256, "sublhs_256");
//   print4x64(subrhs_256, "subrhs_256");

  mul_256 = _mm256_mul_epi32(sublhs_256, subrhs_256);
//   print4x64(mul_256, "mul_256");
  accum10 = _mm256_add_epi64(accum10, mul_256);
//   print4x64(accum10, "accum10");

  // 7 * [3 4 5 6]
  __m256i lhs_7 = _mm256_set1_epi32(x->limbs[8]);
  __m256i rhs_7 = _mm256_set1_epi32(y->limbs[8]);

//   print4x64(lhs_7, "lhs_7");
//   print4x64(rhs_7, "rhs_7");

  sublhs_256 = _mm256_sub_epi32(lhs_7, _mm512_castsi512_si256(lhs_3_10));
  subrhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(rhs_3_10), rhs_7);
//   print4x64(sublhs_256, "sublhs_256");
//   print4x64(subrhs_256, "subrhs_256");
  mul_256 = _mm256_mul_epi32(sublhs_256, subrhs_256);
//   print4x64(mul_256, "mul_256");
  accum10 = _mm256_add_epi64(accum10, mul_256);
//   print4x64(accum10, "accum10");

  // 1 * [2 3 4 5 6 7 8 9]
  __m512i lhs_2_9 = loadu512_extend_32_64((__m256i *) &x->limbs[3]);
  __m512i rhs_2_9 = loadu512_extend_32_64((__m256i *) &y->limbs[3]);

//   print8x64(lhs_2_9, "lhs_2_9");
//   print8x64(rhs_2_9, "rhs_2_9");

  __m512i lhs_1 = _mm512_set1_epi32(x->limbs[2]);
  __m512i rhs_1 = _mm512_set1_epi32(y->limbs[2]);

//   print8x64(lhs_1, "lhs_1");
//   print8x64(rhs_1, "rhs_1");

  sublhs_512 = _mm512_sub_epi32(lhs_1, lhs_2_9);
  subrhs_512 = _mm512_sub_epi32(rhs_2_9, rhs_1);
//   print8x64(sublhs_512, "sublhs_512");
//   print8x64(subrhs_512, "subrhs_512");
  mul_512 = _mm512_mul_epi32(sublhs_512, subrhs_512);
//   print8x64(mul_512, "mul_512");
  accum3 = _mm512_add_epi64(accum3, mul_512);
//   print8x64(accum3, "accum3");

  _mm512_storeu_si512((__m512i*) &temp.limbs[4], accum3);

  // 8 * [2 3 4 5]
  sublhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(lhs_8),
                                _mm512_castsi512_si256(lhs_2_9));
  subrhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(rhs_2_9),
                                _mm512_castsi512_si256(rhs_8));
//   print4x64(sublhs_256, "sublhs_256");
//   print4x64(subrhs_256, "subrhs_256");
  mul_256 = _mm256_mul_epi32(sublhs_256, subrhs_256);
//   print4x64(mul_256, "mul_256");
  accum10 = _mm256_add_epi64(accum10, mul_256);
//   print4x64(accum10, "accum10");
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

  // Two accumulators
  __m512i accum3 = _mm512_setzero();
  __m256i accum10 = _mm256_setzero_si256();

  // Temporaries for sub sub mul
  __m512i sublhs_512, mul_512;
  __m256i sublhs_256, mul_256;

  __m512i blend_var_lhs;
  __m512i blend_const_lhs;

  __m256i blend_var_lhs_256;
  __m256i blend_const_lhs_256;

  __mmask16 low4, low8, low12;

  low4 = _cvtu32_mask16(0xf);
  low8 = _cvtu32_mask16(0xff);
  low12 = _cvtu32_mask16(0xfff);

  // 8 * [6 7] and 2 * [3 4 5 6 7 8]
  __m512i lhs_6_7 = _mm512_castsi128_si512(_mm_loadu_si128((__m128i *) &x->limbs[7]));
//   print8x64(lhs_6_7, "lhs_6_7");

  __m512i lhs_1_8 = _mm512_loadu_si512((__m512i *) &x->limbs[2]);

//   print8x64(lhs_1_8, "lhs_1_8");

  __m512i lhs_8 = _mm512_set1_epi32(x->limbs[9]);

//   print8x64(lhs_8, "lhs_8");

  __m512i lhs_2 = _mm512_set1_epi32(x->limbs[3]);

//   print8x64(lhs_2, "lhs_2");

#if FORCE_BLEND
  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_var_lhs)
      : [mask]  "Yk"((__mmask16) low4),
        [a]     "v"(lhs_6_7),
        [b]     "v"(lhs_1_8));

  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_const_lhs)
      : [mask]  "Yk"((__mmask16) low4),
        [a]     "v"(lhs_8),
        [b]     "v"(lhs_2));
#else
  blend_var_lhs = _mm512_mask_blend_epi32(low4, lhs_6_7, lhs_1_8);

  blend_const_lhs = _mm512_mask_blend_epi32(low4, lhs_8, lhs_2);
#endif

//   print8x64(blend_var_lhs, "blend_var_lhs");
//   print8x64(blend_const_lhs, "blend_const_lhs");

  sublhs_512 = _mm512_sub_epi32(blend_const_lhs, blend_var_lhs);
//   print8x64(sublhs_512, "sublhs_512");

  mul_512 = _mm512_mul_epi32(sublhs_512, sublhs_512);
  accum3 = _mm512_sub_epi64(accum3, mul_512);
//   print8x64(accum3, "accum3");

  // 9 * [5 6 7 8] and 3 * [4 5 6 7]
  __m512i lhs_5_8 = _mm512_castsi256_si512(_mm256_loadu_si256((__m256i *) &x->limbs[6]));

//   print8x64(lhs_5_8, "lhs_5_8");

  __m512i lhs_0_7 = _mm512_loadu_si512((__m512i *) &x->limbs[1]);

//   print8x64(lhs_0_7, "lhs_0_7");

  __m512i lhs_9 = _mm512_set1_epi32(x->limbs[10]);

//   print8x64(lhs_9, "lhs_9");

  __m512i lhs_3 = _mm512_set1_epi32(x->limbs[4]);

//   print8x64(lhs_3, "lhs_3");

#if FORCE_BLEND
  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_var_lhs)
      : [mask]  "Yk"((__mmask16) low8),
        [a]     "v"(lhs_5_8),
        [b]     "v"(lhs_0_7));

  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_const_lhs)
      : [mask]  "Yk"((__mmask16) low8),
        [a]     "v"(lhs_9),
        [b]     "v"(lhs_3));
#else
  blend_var_lhs = _mm512_mask_blend_epi32(low8, lhs_5_8, lhs_0_7);

  blend_const_lhs = _mm512_mask_blend_epi32(low8, lhs_9, lhs_3);
#endif

//   print8x64(blend_var_lhs, "blend_var_lhs");
//   print8x64(blend_const_lhs, "blend_const_lhs");

  sublhs_512 = _mm512_sub_epi32(blend_const_lhs, blend_var_lhs);
//   print8x64(sublhs_512, "sublhs_512");

  mul_512 = _mm512_mul_epi32(sublhs_512, sublhs_512);
//   print8x64(mul_512, "mul_512");
  accum3 = _mm512_sub_epi64(accum3, mul_512);
//   print8x64(accum3, "accum3");

  // 9 * [1 2 3 4] Done here because we have loaded everything already.
  sublhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(lhs_9),
                                _mm512_castsi512_si256(lhs_1_8));
//   print4x64(sublhs_256, "sublhs_256");

  mul_256 = _mm256_mul_epi32(sublhs_256, sublhs_256);
  accum10 = _mm256_sub_epi64(accum10, mul_256);
//   print4x64(accum10, "accum10");

  // 10 * [4 5 6 7 8 9] and 4 * [5 6]
  __m512i lhs_4_9 = _mm512_mask_loadu_epi32(_mm512_setzero(), low12, &x->limbs[5]);

//   print8x64(lhs_4_9, "lhs_4_9");

  __m512i lhs_10_6 = _mm512_loadu_si512((__m512i *) &x->limbs[0]);

//   print8x64(lhs_10_6, "lhs_10_6");

  __m512i lhs_10 = _mm512_set1_epi32(x->limbs[0]);

//   print8x64(lhs_10, "lhs_10");

  __m512i lhs_4 = _mm512_set1_epi32(x->limbs[5]);

//   print8x64(lhs_4, "lhs_4");

#if FORCE_BLEND
  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_var_lhs)
      : [mask]  "Yk"((__mmask16) low12),
        [a]     "v"(lhs_4_9),
        [b]     "v"(lhs_10_6));

  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_const_lhs)
      : [mask]  "Yk"((__mmask16) low12),
        [a]     "v"(lhs_10),
        [b]     "v"(lhs_4));
#else
  blend_var_lhs = _mm512_mask_blend_epi32(low12, lhs_4_9, lhs_10_6);

  blend_const_lhs = _mm512_mask_blend_epi32(low12, lhs_10, lhs_4);
#endif

//   print8x64(blend_var_lhs, "blend_var_lhs");
//   print8x64(blend_const_lhs, "blend_const_lhs");

  sublhs_512 = _mm512_sub_epi32(blend_const_lhs, blend_var_lhs);

//   print8x64(sublhs_512, "sublhs_512");

  mul_512 = _mm512_mul_epi32(sublhs_512, sublhs_512);
//   print8x64(mul_512, "mul_512");
  accum3 = _mm512_sub_epi64(accum3, mul_512);
//   print8x64(accum3, "accum3");

  // 10 * [0 1 2 3] Done here because we have loaded everything already.
  sublhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(lhs_10),
                                _mm512_castsi512_si256(lhs_0_7));
//   print4x64(sublhs_256, "sublhs_256");

  mul_256 = _mm256_mul_epi32(sublhs_256, sublhs_256);
//   print4x64(mul_256, "mul_256");
  accum10 = _mm256_sub_epi64(accum10, mul_256);
//   print4x64(accum10, "accum10");

  // 0 * [3 4 5 6 7 8 9 10]
  __m512i lhs_3_10 = _mm512_loadu_si512((__m512i *) &x->limbs[4]);

//   print8x64(lhs_3_10, "lhs_3_10");

  __m512i lhs_0 = _mm512_set1_epi32(x->limbs[1]);

//   print8x64(lhs_0, "lhs_0");

  sublhs_512 = _mm512_sub_epi32(lhs_0, lhs_3_10);
//   print8x64(sublhs_512, "sublhs_512");

  mul_512 = _mm512_mul_epi32(sublhs_512, sublhs_512);
//   print8x64(mul_512, "mul_512");
  accum3 = _mm512_sub_epi64(accum3, mul_512);
//   print8x64(accum3, "accum3");

  // 6 * [4 5] and 0 * [1 2]
  __m256i lhs_6 = _mm256_set1_epi32(x->limbs[7]);

//   print4x64(lhs_6, "lhs_6");

  blend_var_lhs_256 = _mm256_blend_epi32(_mm512_castsi512_si256(lhs_10_6),
                                         _mm512_castsi512_si256(lhs_4_9),
                                         0x0f);

  blend_const_lhs_256 = _mm256_blend_epi32(_mm512_castsi512_si256(lhs_0),
                                           lhs_6,
                                           0x0f);
//   print4x64(blend_var_lhs_256, "blend_var_lhs_256");
//   print4x64(blend_const_lhs_256, "blend_const_lhs_256");

  sublhs_256 = _mm256_sub_epi32(blend_const_lhs_256, blend_var_lhs_256);
//   print4x64(sublhs_256, "sublhs_256");

  mul_256 = _mm256_mul_epi32(sublhs_256, sublhs_256);
//   print4x64(mul_256, "mul_256");
  accum10 = _mm256_sub_epi64(accum10, mul_256);
//   print4x64(accum10, "accum10");

  // 7 * [3 4 5 6]
  __m256i lhs_7 = _mm256_set1_epi32(x->limbs[8]);

//   print4x64(lhs_7, "lhs_7");

  sublhs_256 = _mm256_sub_epi32(lhs_7, _mm512_castsi512_si256(lhs_3_10));
//   print4x64(sublhs_256, "sublhs_256");
  mul_256 = _mm256_mul_epi32(sublhs_256, sublhs_256);
//   print4x64(mul_256, "mul_256");
  accum10 = _mm256_sub_epi64(accum10, mul_256);
//   print4x64(accum10, "accum10");

  // 1 * [2 3 4 5 6 7 8 9]
  __m512i lhs_2_9 = _mm512_loadu_si512((__m512i *) &x->limbs[3]);

//   print8x64(lhs_2_9, "lhs_2_9");

  __m512i lhs_1 = _mm512_set1_epi32(x->limbs[2]);

//   print8x64(lhs_1, "lhs_1");

  sublhs_512 = _mm512_sub_epi32(lhs_1, lhs_2_9);
//   print8x64(sublhs_512, "sublhs_512");
  mul_512 = _mm512_mul_epi32(sublhs_512, sublhs_512);
//   print8x64(mul_512, "mul_512");
  accum3 = _mm512_sub_epi64(accum3, mul_512);
//   print8x64(accum3, "accum3");

  _mm512_storeu_si512((__m512i*) &temp.limbs[4], accum3);

  // 8 * [2 3 4 5]
  sublhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(lhs_8),
                                _mm512_castsi512_si256(lhs_2_9));
//   print4x64(sublhs_256, "sublhs_256");
  mul_256 = _mm256_mul_epi32(sublhs_256, sublhs_256);
//   print4x64(mul_256, "mul_256");
  accum10 = _mm256_sub_epi64(accum10, mul_256);
//   print4x64(accum10, "accum10");
  _mm256_store_si256((__m256i*) &temp.limbs[0], accum10);
  reduce_step_wide(&temp, &temp);
  reduce_step_wide(result, &temp);
}

// Square a narrow residue and produce a wide result. The result is reduced to
// 32 bits but not narrowed for performance reasons.
void square_narrow(
  residue_wide_t *result, const residue_narrow_t *x) {

  residue_wide_t temp;

  // Two accumulators
  __m512i accum3 = _mm512_setzero();
  __m256i accum10 = _mm256_setzero_si256();

  // Temporaries for sub sub mul
  __m512i sublhs_512, mul_512;
  __m256i sublhs_256, mul_256;

  __m512i blend_var_lhs;
  __m512i blend_const_lhs;

  __m256i blend_var_lhs_256;
  __m256i blend_const_lhs_256;

  __mmask16 low4, low8, low12;

  low4 = _cvtu32_mask16(0xf);
  low8 = _cvtu32_mask16(0xff);
  low12 = _cvtu32_mask16(0xfff);

  // 8 * [6 7] and 2 * [3 4 5 6 7 8]
  __m512i lhs_6_7 = _mm512_castsi256_si512(loadu_extend_32_64((__m128i *) &x->limbs[7]));
//   print8x64(lhs_6_7, "lhs_6_7");

  __m512i lhs_1_8 = loadu512_extend_32_64((__m256i *) &x->limbs[2]);

//   print8x64(lhs_1_8, "lhs_1_8");

  __m512i lhs_8 = _mm512_set1_epi32(x->limbs[9]);

//   print8x64(lhs_8, "lhs_8");

  __m512i lhs_2 = _mm512_set1_epi32(x->limbs[3]);

//   print8x64(lhs_2, "lhs_2");

#if FORCE_BLEND
  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_var_lhs)
      : [mask]  "Yk"((__mmask16) low4),
        [a]     "v"(lhs_6_7),
        [b]     "v"(lhs_1_8));

  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_const_lhs)
      : [mask]  "Yk"((__mmask16) low4),
        [a]     "v"(lhs_8),
        [b]     "v"(lhs_2));
#else
  blend_var_lhs = _mm512_mask_blend_epi32(low4, lhs_6_7, lhs_1_8);

  blend_const_lhs = _mm512_mask_blend_epi32(low4, lhs_8, lhs_2);
#endif

//   print8x64(blend_var_lhs, "blend_var_lhs");
//   print8x64(blend_const_lhs, "blend_const_lhs");

  sublhs_512 = _mm512_sub_epi32(blend_const_lhs, blend_var_lhs);
//   print8x64(sublhs_512, "sublhs_512");

  mul_512 = _mm512_mul_epi32(sublhs_512, sublhs_512);
  accum3 = _mm512_sub_epi64(accum3, mul_512);
//   print8x64(accum3, "accum3");

  // 9 * [5 6 7 8] and 3 * [4 5 6 7]
  __m512i lhs_5_8 = _mm512_castsi256_si512(loadu_extend_32_64((__m128i *) &x->limbs[6]));

//   print8x64(lhs_5_8, "lhs_5_8");

  __m512i lhs_0_7 = loadu512_extend_32_64((__m256i *) &x->limbs[1]);

//   print8x64(lhs_0_7, "lhs_0_7");

  __m512i lhs_9 = _mm512_set1_epi32(x->limbs[10]);

//   print8x64(lhs_9, "lhs_9");

  __m512i lhs_3 = _mm512_set1_epi32(x->limbs[4]);

//   print8x64(lhs_3, "lhs_3");

#if FORCE_BLEND
  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_var_lhs)
      : [mask]  "Yk"((__mmask16) low8),
        [a]     "v"(lhs_5_8),
        [b]     "v"(lhs_0_7));

  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_const_lhs)
      : [mask]  "Yk"((__mmask16) low8),
        [a]     "v"(lhs_9),
        [b]     "v"(lhs_3));
#else
  blend_var_lhs = _mm512_mask_blend_epi32(low8, lhs_5_8, lhs_0_7);

  blend_const_lhs = _mm512_mask_blend_epi32(low8, lhs_9, lhs_3);
#endif

//   print8x64(blend_var_lhs, "blend_var_lhs");
//   print8x64(blend_const_lhs, "blend_const_lhs");

  sublhs_512 = _mm512_sub_epi32(blend_const_lhs, blend_var_lhs);
//   print8x64(sublhs_512, "sublhs_512");

  mul_512 = _mm512_mul_epi32(sublhs_512, sublhs_512);
//   print8x64(mul_512, "mul_512");
  accum3 = _mm512_sub_epi64(accum3, mul_512);
//   print8x64(accum3, "accum3");

  // 9 * [1 2 3 4] Done here because we have loaded everything already.
  sublhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(lhs_9),
                                _mm512_castsi512_si256(lhs_1_8));
//   print4x64(sublhs_256, "sublhs_256");

  mul_256 = _mm256_mul_epi32(sublhs_256, sublhs_256);
  accum10 = _mm256_sub_epi64(accum10, mul_256);
//   print4x64(accum10, "accum10");

  // 10 * [4 5 6 7 8 9] and 4 * [5 6]
  __m512i lhs_4_9 = loadu512_mask_extend_32_64((__m256i *) &x->limbs[5], _cvtu32_mask8(0x3f));

//   print8x64(lhs_4_9, "lhs_4_9");

  __m512i lhs_10_6 = loadu512_extend_32_64((__m256i *) &x->limbs[0]);

//   print8x64(lhs_10_6, "lhs_10_6");

  __m512i lhs_10 = _mm512_set1_epi32(x->limbs[0]);

//   print8x64(lhs_10, "lhs_10");

  __m512i lhs_4 = _mm512_set1_epi32(x->limbs[5]);

//   print8x64(lhs_4, "lhs_4");

#if FORCE_BLEND
  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_var_lhs)
      : [mask]  "Yk"((__mmask16) low12),
        [a]     "v"(lhs_4_9),
        [b]     "v"(lhs_10_6));

  __asm__(
      "vpblendmd %[a], %[b], %[blend] %{%[mask]%}"
      : [blend] "=v"(blend_const_lhs)
      : [mask]  "Yk"((__mmask16) low12),
        [a]     "v"(lhs_10),
        [b]     "v"(lhs_4));
#else
  blend_var_lhs = _mm512_mask_blend_epi32(low12, lhs_4_9, lhs_10_6);

  blend_const_lhs = _mm512_mask_blend_epi32(low12, lhs_10, lhs_4);
#endif

//   print8x64(blend_var_lhs, "blend_var_lhs");
//   print8x64(blend_const_lhs, "blend_const_lhs");

  sublhs_512 = _mm512_sub_epi32(blend_const_lhs, blend_var_lhs);

//   print8x64(sublhs_512, "sublhs_512");

  mul_512 = _mm512_mul_epi32(sublhs_512, sublhs_512);
//   print8x64(mul_512, "mul_512");
  accum3 = _mm512_sub_epi64(accum3, mul_512);
//   print8x64(accum3, "accum3");

  // 10 * [0 1 2 3] Done here because we have loaded everything already.
  sublhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(lhs_10),
                                _mm512_castsi512_si256(lhs_0_7));
//   print4x64(sublhs_256, "sublhs_256");

  mul_256 = _mm256_mul_epi32(sublhs_256, sublhs_256);
//   print4x64(mul_256, "mul_256");
  accum10 = _mm256_sub_epi64(accum10, mul_256);
//   print4x64(accum10, "accum10");

  // 0 * [3 4 5 6 7 8 9 10]
  __m512i lhs_3_10 = loadu512_extend_32_64((__m256i *) &x->limbs[4]);

//   print8x64(lhs_3_10, "lhs_3_10");

  __m512i lhs_0 = _mm512_set1_epi32(x->limbs[1]);

//   print8x64(lhs_0, "lhs_0");

  sublhs_512 = _mm512_sub_epi32(lhs_0, lhs_3_10);
//   print8x64(sublhs_512, "sublhs_512");

  mul_512 = _mm512_mul_epi32(sublhs_512, sublhs_512);
//   print8x64(mul_512, "mul_512");
  accum3 = _mm512_sub_epi64(accum3, mul_512);
//   print8x64(accum3, "accum3");

  // 6 * [4 5] and 0 * [1 2]
  __m256i lhs_6 = _mm256_set1_epi32(x->limbs[7]);

//   print4x64(lhs_6, "lhs_6");

  blend_var_lhs_256 = _mm256_blend_epi32(_mm512_castsi512_si256(lhs_10_6),
                                         _mm512_castsi512_si256(lhs_4_9),
                                         0x0f);

  blend_const_lhs_256 = _mm256_blend_epi32(_mm512_castsi512_si256(lhs_0),
                                           lhs_6,
                                           0x0f);
//   print4x64(blend_var_lhs_256, "blend_var_lhs_256");
//   print4x64(blend_const_lhs_256, "blend_const_lhs_256");

  sublhs_256 = _mm256_sub_epi32(blend_const_lhs_256, blend_var_lhs_256);
//   print4x64(sublhs_256, "sublhs_256");

  mul_256 = _mm256_mul_epi32(sublhs_256, sublhs_256);
//   print4x64(mul_256, "mul_256");
  accum10 = _mm256_sub_epi64(accum10, mul_256);
//   print4x64(accum10, "accum10");

  // 7 * [3 4 5 6]
  __m256i lhs_7 = _mm256_set1_epi32(x->limbs[8]);

//   print4x64(lhs_7, "lhs_7");

  sublhs_256 = _mm256_sub_epi32(lhs_7, _mm512_castsi512_si256(lhs_3_10));
//   print4x64(sublhs_256, "sublhs_256");
  mul_256 = _mm256_mul_epi32(sublhs_256, sublhs_256);
//   print4x64(mul_256, "mul_256");
  accum10 = _mm256_sub_epi64(accum10, mul_256);
//   print4x64(accum10, "accum10");

  // 1 * [2 3 4 5 6 7 8 9]
  __m512i lhs_2_9 = loadu512_extend_32_64((__m256i *) &x->limbs[3]);

//   print8x64(lhs_2_9, "lhs_2_9");

  __m512i lhs_1 = _mm512_set1_epi32(x->limbs[2]);

//   print8x64(lhs_1, "lhs_1");

  sublhs_512 = _mm512_sub_epi32(lhs_1, lhs_2_9);
//   print8x64(sublhs_512, "sublhs_512");
  mul_512 = _mm512_mul_epi32(sublhs_512, sublhs_512);
//   print8x64(mul_512, "mul_512");
  accum3 = _mm512_sub_epi64(accum3, mul_512);
//   print8x64(accum3, "accum3");

  _mm512_storeu_si512((__m512i*) &temp.limbs[4], accum3);

  // 8 * [2 3 4 5]
  sublhs_256 = _mm256_sub_epi32(_mm512_castsi512_si256(lhs_8),
                                _mm512_castsi512_si256(lhs_2_9));
//   print4x64(sublhs_256, "sublhs_256");
  mul_256 = _mm256_mul_epi32(sublhs_256, sublhs_256);
//   print4x64(mul_256, "mul_256");
  accum10 = _mm256_sub_epi64(accum10, mul_256);
//   print4x64(accum10, "accum10");
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

  __attribute__((__aligned__(64)))
  static const int64_t carry_permute_3[8] = {
    0xb, 0x0, 0x1, 0x2, 0x3, 0x4, 0x5, 0x6,
  };

  static const int64_t carry_permute_10[8] = {
    0x6, 0x8, 0x9, 0xa, 0xc, 0xd, 0xe, 0xf,
  };

  __m512i accum3;
  __m256i accum10;

  __m512i carry_permute_3_512 = _mm512_load_si512(carry_permute_3);
  __m512i carry_permute_10_512 = _mm512_load_si512(carry_permute_10);
  __m512i mask = _mm512_set1_epi64(0x3ffffff);

  accum3 = _mm512_loadu_si512((__m512i*) &x->limbs[4]);
//   print8x64(accum3, "accum3");
  __m512i shift3 = _mm512_srai_epi64(accum3, 26);
//   print8x64(shift3, "shift3");
  accum3 = _mm512_and_si512(accum3, mask);
//   print8x64(accum3, "accum3");
  accum3 = _mm512_sub_epi64(accum3, shift3);
//   print8x64(accum3, "accum3");
  __m512i shift_error3 = _mm512_slli_epi64(shift3, 4);
//   print8x64(shift_error3, "shift_error3");
  accum3 = _mm512_add_epi64(accum3, shift_error3);
//   print8x64(accum3, "accum3");

  accum10 = _mm256_load_si256((__m256i*) &x->limbs[0]);
//   print4x64(accum10, "accum10");
  __m256i shift10 = _mm256_srai_epi64(accum10, 26);
//   print4x64(shift10, "shift10");
  accum10 = _mm256_and_si256(accum10, _mm512_castsi512_si256(mask));
//   print4x64(accum10, "accum10");
  accum10 = _mm256_sub_epi64(accum10, shift10);
//   print4x64(accum10, "accum10");
  __m256i shift_error10 = _mm256_slli_epi64(shift10, 4);
//   print4x64(shift_error10, "shift_error10");
  accum10 = _mm256_add_epi64(accum10, shift_error10);
//   print4x64(accum10, "accum10");

  __m512i carry3 = _mm512_permutex2var_epi64(
      shift3, carry_permute_3_512, _mm512_castsi256_si512(shift10));
//   print8x64(carry3, "carry3");
  __m256i carry10 = _mm512_castsi512_si256(
      _mm512_permutex2var_epi64(
          shift3, carry_permute_10_512, _mm512_castsi256_si512(shift10)));
//   print4x64(carry10, "carry10");

  accum3 = _mm512_add_epi64(accum3, carry3);
//   print8x64(accum3, "accum3");
  _mm512_storeu_si512((__m512i*) &result->limbs[4], accum3);
  accum10 = _mm256_add_epi64(accum10, carry10);
//   print4x64(accum10, "accum10");
  _mm256_store_si256((__m256i*) &result->limbs[0], accum10);
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

static void raise_to_t2(
  residue_wide_t *result, const residue_wide_t *x) {
  // t^2 = 0xfffff880000e1
  // zi = z^(2^i - 1), z1 = x
  residue_wide_t z2;
  residue_wide_t z3;
  residue_wide_t z5;
  residue_wide_t z10;
  residue_wide_t z20;
  residue_wide_t result_t;

  square_wide(&z2, x);
  mul_wide(&z2, &z2, x);
  square_wide(&z3, &z2);
  mul_wide(&z3, &z3, x);
  nsquare_wide(&z5, &z3, 2);
  mul_wide(&z5, &z5, &z2);
  nsquare_wide(&z10, &z5, 5);
  mul_wide(&z10, &z10, &z5);
  nsquare_wide(&z20, &z10, 10);
  mul_wide(&z20, &z20, &z10);
  square_wide(&result_t, &z20);
  mul_wide(&result_t, &result_t, x);
  nsquare_wide(&result_t, &result_t, 4);
  mul_wide(&result_t, &result_t, x);
  // 22 = 3 for zeros in 8, 16 for zeros in 0000, 3 to make room for e.
  nsquare_wide(&result_t, &result_t, 22);
  mul_wide(&result_t, &result_t, &z3);
  nsquare_wide(&result_t, &result_t, 5);
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
  // x^2 (trades a multiply for a square)
  residue_wide_t x2;
  // rho_k = x^((t^k - 1)/(t - 1))
  // rho_1 = x
  residue_wide_t rho_2, rho_4, rho_8, rho_9;
  residue_wide_t result_t;

  raise_to_t_minus_1_over_4(&x_t_minus_1_over_4, x);
  nsquare_wide(&x_t_minus_1, &x_t_minus_1_over_4, 2);
  square_wide(&x2, x);
  mul_wide(&rho_2, &x_t_minus_1, &x2);
  raise_to_t2(&rho_4, &rho_2);
  mul_wide(&rho_4, &rho_4, &rho_2);
  raise_to_t2(&rho_8, &rho_4);
  raise_to_t2(&rho_8, &rho_8);
  mul_wide(&rho_8, &rho_8, &rho_4);
  raise_to_t(&rho_9, &rho_8);
  mul_wide(&rho_9, &rho_9, x);
  raise_to_t2(&result_t, &rho_9);
  mul_wide(result, &result_t, &x_t_minus_1);
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
