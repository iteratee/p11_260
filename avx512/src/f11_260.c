#include <stdint.h>
#include "f11_260.h"
#include "mul_inline.h"
#include "emmintrin.h"
#include "immintrin.h"

residue_narrow_t zero_narrow = {0};
residue_narrow_t one_narrow = {
  .limbs = {1},
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

// Reduce to 10 limbs. Useful for debugging.
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
    result->limbs[i] = temp.limbs[i] - temp.limbs[10];
  }
}

// Reduce to unique representative.
// This is expensive. Only used for final signature or DH Key
void narrow_complete(
  residue_narrow_reduced_t *result, const residue_narrow_t * __restrict w) {

  residue_narrow_t temp;
  for (int i = 0; i < NLIMBS; ++i) {
    temp.limbs[i] = w->limbs[i] - w->limbs[10];
  }

  // This may be combined with the final reduction from a multiply.
  reduce_step_narrow(&temp, &temp);

  int gt_mask = 0;
  int lt_mask = 0;
  int32_t limit[NLIMBS];
  for (int i = 0; i < NLIMBS; ++i) {
    temp.limbs[i] = temp.limbs[i] - temp.limbs[10];
    temp.limbs[i] += 1 & gt_mask;
    temp.limbs[i] -= 1 & lt_mask;
    gt_mask = -(temp.limbs[i] > T);
    lt_mask = -(temp.limbs[i] < 0);
    temp.limbs[i] -= (T & gt_mask);
    temp.limbs[i] += (T & lt_mask);
  }
  for (int i = 0; i < NLIMBS_REDUCED; ++i) {
    temp.limbs[i] -= temp.limbs[10];
    limit[i] = T;
  }
  int64_t all_t = -1;
  for (int i = NLIMBS_REDUCED - 2; i >= 0; --i) {
    all_t &= -(temp.limbs[i+1] == T);
    limit[i] -= 1 & (~all_t);
  }
  gt_mask = 0;
  lt_mask = 0;
  for (int i = 0; i < NLIMBS_REDUCED; ++i) {
    temp.limbs[i] += 1 & gt_mask;
    temp.limbs[i] -= 1 & lt_mask;
    gt_mask = -(temp.limbs[i] > limit[i]);
    lt_mask = -(temp.limbs[i] < 0);
    temp.limbs[i] -= (T & gt_mask);
    temp.limbs[i] += (T & lt_mask);
    result->limbs[i] = temp.limbs[i];
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
    temp.limbs[i] = w->limbs[i] - w->limbs[10];
  }

  // This may be combined with the final reduction from a multiply.
  reduce_step_narrow(&temp, &temp);

  int gt_mask = 0;
  int lt_mask = 0;
  for (int i = 0; i < NLIMBS; ++i) {
    temp.limbs[i] = temp.limbs[i] - temp.limbs[10];
    temp.limbs[i] += 1 & gt_mask;
    temp.limbs[i] -= 1 & lt_mask;
    gt_mask = -(temp.limbs[i] > T);
    lt_mask = -(temp.limbs[i] < 0);
    temp.limbs[i] -= (T & gt_mask);
    temp.limbs[i] += (T & lt_mask);
  }
  for (int i = 0; i < NLIMBS - 1; ++i) {
    temp.limbs[i] -= temp.limbs[10];
  }
  gt_mask = 0;
  lt_mask = 0;
  for (int i = 0; i < NLIMBS_REDUCED; ++i) {
    temp.limbs[i] += 1 & gt_mask;
    temp.limbs[i] -= 1 & lt_mask;
    gt_mask = -(temp.limbs[i] > T);
    lt_mask = -(temp.limbs[i] < 0);
    temp.limbs[i] -= (T & gt_mask);
    temp.limbs[i] += (T & lt_mask);
    result->limbs[i] = temp.limbs[i];
  }
}

int is_odd(residue_narrow_reduced_t *x) {
  int result = 0;
  for (int i = 0; i < NLIMBS_REDUCED; ++i) {
    result ^= x->limbs[i] & 0x1;
  }
  return result;
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

void reduce_step_narrow(
  residue_narrow_t *result, const residue_narrow_t *x) {
  return reduce_step_narrow_i(result, x);
}

void reduce_step_wide(
  residue_wide_t *result, const residue_wide_t *x) {
  return reduce_step_wide_i(result, x);
}

void mul_narrow(
  residue_narrow_t *result, const residue_narrow_t *x,
  const residue_narrow_t *y) {
  return mul_narrow_i(result, x, y);
}

void square_narrow(
  residue_narrow_t *result, const residue_narrow_t *x) {
  return square_narrow_i(result, x);
}

// Produce a 64-bit residue
void widen(
  residue_wide_t *result, const residue_narrow_t * __restrict x) {
  __m256i wide10 = loadu_extend_32_64((__m128i *) x);
  __m512i wide3 = loadu512_extend_32_64((__m256i *) &x[4]);
  _mm256_store_si256((__m256i*) &result->limbs[0], wide10);
  _mm512_storeu_si512((__m512i*) &result->limbs[4], wide3);
}

// Subtract 2 12x32-bit residues.
void sub_narrow(
  residue_narrow_t *result, const residue_narrow_t * __restrict x,
  const residue_narrow_t * __restrict y) {

  __m512i lhs = _mm512_load_si512((__m512i*) &x->limbs[0]);
  __m512i rhs = _mm512_load_si512((__m512i*) &y->limbs[0]);
  __m512i sub = _mm512_sub_epi32(lhs, rhs);
  _mm512_store_si512((__m512i*) &result->limbs[0], sub);
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

  __m512i lhs = _mm512_load_si512((__m512i*) &x->limbs[0]);
  __m512i zero = _mm512_setzero();
  __m512i neg = _mm512_sub_epi32(zero, lhs);
  _mm512_store_si512((__m512i*) &result->limbs[0], neg);
}

// Add 2 12x32-bit residues.
void add_narrow(
  residue_narrow_t *result, const residue_narrow_t * __restrict x,
  const residue_narrow_t * __restrict y) {

  __m512i lhs = _mm512_load_si512((__m512i*) &x->limbs[0]);
  __m512i rhs = _mm512_load_si512((__m512i*) &y->limbs[0]);
  __m512i add = _mm512_add_epi32(lhs, rhs);
  _mm512_store_si512((__m512i*) &result->limbs[0], add);
}

// Scale a narrow residue by 2.
void double_narrow(
  residue_narrow_t *result, const residue_narrow_t *x) {

  __m512i lhs = _mm512_load_si512((__m512i*) &x->limbs[0]);
  __m512i dub = _mm512_slli_epi32(lhs, 1);
  _mm512_store_si512((__m512i*) &result->limbs[0], dub);
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

// static void print16x32(__m512i x, const char * preamble) {
//   uint32_t x_vals[16];
//   memcpy(x_vals, &x, sizeof(x_vals));
//   printf("%s\n", preamble);
//   for (int i = 0; i < 16; ++i) {
//     printf("%#x\n", x_vals[i]);
//   }
// }

// The swaps below trade 32 bit words within 128 bit lanes
// in low endian order: 01 00 11 10
// in big endian order 10 11 00 01 = 0xb1

static const int32_t permute_final_result[16] = {
  0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30,
};

// Multiply a narrow residue by a small constant. The result is reduced to 32
// bits.
void mul_narrow_const(
  residue_narrow_t *result, const residue_narrow_t *x, int32_t d) {

  residue_wide_t temp;
  for (int i = 0; i < NLIMBS; ++i) {
    temp.limbs[i] = ((uint64_t) x->limbs[i]) * d;
  }
  reduce_step_wide(&temp, &temp);
  __m512i _permute_final_result = _mm512_load_si512(permute_final_result);
  __m512i accum0 = _mm512_load_si512((__m512i*) &temp.limbs[0]);
  __m256i accum8 = _mm256_load_si256((__m256i*) &temp.limbs[8]);
  __m512i final_result = _mm512_permutex2var_epi32(
    accum0, _permute_final_result, _mm512_castsi256_si512(accum8));
  _mm512_store_si512((__m512i*) &result->limbs[0], final_result);
}

// Takes advantage of the fact that if a residue z *is zero* then after setting
// one coefficient to T/2, all the remaining coefficients should be near to
// T/2. They should therefore resolve all carries in a single step, and all be
// equal to the same value. Some other value may not reduce completely, but this
// is fine, we will know it is not zero.
int equal_narrow(const residue_narrow_t *x, const residue_narrow_t *y) {
  residue_narrow_t temp;

  sub_narrow(&temp, x, y);
  int32_t delta = -temp.limbs[0] + (T / 2);
  for (int i = 0; i < NLIMBS; ++i) {
    temp.limbs[i] += delta;
  }

  reduce_step_narrow(&temp, &temp);

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

static inline void nsquare_narrow(
  residue_narrow_t *result, const residue_narrow_t *x, int n) {

  square_narrow_i(result, x);
  for (int i = 1; i < n; ++i) {
    square_narrow_i(result, result);
  }
}

static void raise_to_t(
  residue_narrow_t *result, const residue_narrow_t *x) {
  // zi = z^(2^i - 1), z1 = x
  residue_narrow_t z2;
  residue_narrow_t z3;
  residue_narrow_t z5;
  residue_narrow_t z10;
  residue_narrow_t z11;
  residue_narrow_t z22;
  residue_narrow_t result_t;

  square_narrow_i(&z2, x);
  mul_narrow_i(&z2, &z2, x);
  square_narrow_i(&z3, &z2);
  mul_narrow_i(&z3, &z3, x);
  nsquare_narrow(&z5, &z3, 2);
  mul_narrow_i(&z5, &z5, &z2);
  nsquare_narrow(&z10, &z5, 5);
  mul_narrow_i(&z10, &z10, &z5);
  square_narrow_i(&z11, &z10);
  mul_narrow_i(&z11, &z11, x);
  nsquare_narrow(&z22, &z11, 11);
  mul_narrow_i(&z22, &z22, &z11);
  nsquare_narrow(&result_t, &z22, 4);
  mul_narrow_i(result, &result_t, x);
}

static void raise_to_t2(
  residue_narrow_t *result, const residue_narrow_t *x) {
  // t^2 = 0xfffff880000e1
  // zi = z^(2^i - 1), z1 = x
  residue_narrow_t z2;
  residue_narrow_t z3;
  residue_narrow_t z5;
  residue_narrow_t z10;
  residue_narrow_t z20;
  residue_narrow_t result_t;

  square_narrow_i(&z2, x);
  mul_narrow_i(&z2, &z2, x);
  square_narrow_i(&z3, &z2);
  mul_narrow_i(&z3, &z3, x);
  nsquare_narrow(&z5, &z3, 2);
  mul_narrow_i(&z5, &z5, &z2);
  nsquare_narrow(&z10, &z5, 5);
  mul_narrow_i(&z10, &z10, &z5);
  nsquare_narrow(&z20, &z10, 10);
  mul_narrow_i(&z20, &z20, &z10);
  square_narrow_i(&result_t, &z20);
  mul_narrow_i(&result_t, &result_t, x);
  nsquare_narrow(&result_t, &result_t, 4);
  mul_narrow_i(&result_t, &result_t, x);
  // 22 = 3 for zeros in 8, 16 for zeros in 0000, 3 to make room for e.
  nsquare_narrow(&result_t, &result_t, 22);
  mul_narrow_i(&result_t, &result_t, &z3);
  nsquare_narrow(&result_t, &result_t, 5);
  mul_narrow_i(result, &result_t, x);
}

static void raise_to_phi_t(
  residue_narrow_t *result, const residue_narrow_t *x, int n) {
  residue_narrow_t temp;

  raise_to_t(&temp, x);

  for (int i = 1; i < n; ++i) {
    mul_narrow(&temp, &temp, x);
    raise_to_t(&temp, &temp);
  }

  mul_narrow(result, &temp, x);
}

static void raise_to_t_minus_1_over_4(
  residue_narrow_t *result, const residue_narrow_t *x) {
  // zi = z^(2^i - 1), z1 = x
  residue_narrow_t z2;
  residue_narrow_t z3;
  residue_narrow_t z5;
  residue_narrow_t z10;
  residue_narrow_t z11;
  residue_narrow_t z22;

  square_narrow_i(&z2, x);
  mul_narrow_i(&z2, &z2, x);
  square_narrow_i(&z3, &z2);
  mul_narrow_i(&z3, &z3, x);
  nsquare_narrow(&z5, &z3, 2);
  mul_narrow_i(&z5, &z5, &z2);
  nsquare_narrow(&z10, &z5, 5);
  mul_narrow_i(&z10, &z10, &z5);
  square_narrow_i(&z11, &z10);
  mul_narrow_i(&z11, &z11, x);
  nsquare_narrow(&z22, &z11, 11);
  mul_narrow_i(&z22, &z22, &z11);
  nsquare_narrow(result, &z22, 2);
}

static void raise_to_p_minus_3_over_4(
  residue_narrow_t *result, const residue_narrow_t *x) {

  residue_narrow_t z4; //z to (t-1)/4
  residue_narrow_t z2; //z to (t-1)/2
  residue_narrow_t z3_4; //z to (3t+1)/4
  residue_narrow_t y_small;
  residue_narrow_t y, y_t4_y;
  residue_narrow_t raised;

  raise_to_t_minus_1_over_4(&z4, x);
  square_narrow(&z2, &z4);
  mul_narrow(&z3_4, &z2, &z4);
  mul_narrow(&z3_4, &z3_4, x);
  raise_to_t(&raised, &z4);
  mul_narrow(&y_small, &z2, &raised);
  raise_to_t(&raised, &y_small);
  mul_narrow(&y, &z3_4, &raised);
  raise_to_t(&raised, &y);
  raise_to_t(&raised, &raised);
  raise_to_t(&raised, &raised);
  raise_to_t(&raised, &raised);
  mul_narrow(&y_t4_y, &raised, &y);
  raise_to_t(&raised, &y_t4_y);
  raise_to_t(&raised, &raised);
  raise_to_t(&raised, &raised);
  mul_narrow(result, &raised, &y_small);
}

int sqrt_inv_narrow(
  residue_narrow_t *result, const residue_narrow_t * __restrict x,
  const residue_narrow_t * __restrict y) {
  residue_narrow_t xy;
  residue_narrow_t y2;
  residue_narrow_t xy3;
  residue_narrow_t xy3_p_3_over_4;
  residue_narrow_t cand2;
  residue_narrow_t should_be_x;

  square_narrow(&y2, y);
  mul_narrow(&xy, x, y);
  mul_narrow(&xy3, &xy, &y2);
  raise_to_p_minus_3_over_4(&xy3_p_3_over_4, &xy3);
  mul_narrow(result, &xy, &xy3_p_3_over_4);
  square_narrow(&cand2, result);
  mul_narrow(&should_be_x, y, &cand2);

  return equal_narrow(&should_be_x, x);
}

void invert_narrow(
  residue_narrow_t *result, const residue_narrow_t * __restrict x) {

  residue_narrow_t x_t_minus_1_over_4;
  residue_narrow_t x_t_minus_1;
  // x^2 (trades a multiply for a square)
  residue_narrow_t x2;
  // rho_k = x^((t^k - 1)/(t - 1))
  // rho_1 = x
  residue_narrow_t rho_2, rho_4, rho_8, rho_9;
  residue_narrow_t result_t;

  raise_to_t_minus_1_over_4(&x_t_minus_1_over_4, x);
  nsquare_narrow(&x_t_minus_1, &x_t_minus_1_over_4, 2);
  square_narrow_i(&x2, x);
  mul_narrow_i(&rho_2, &x_t_minus_1, &x2);
  raise_to_t2(&rho_4, &rho_2);
  mul_narrow_i(&rho_4, &rho_4, &rho_2);
  raise_to_t2(&rho_8, &rho_4);
  raise_to_t2(&rho_8, &rho_8);
  mul_narrow_i(&rho_8, &rho_8, &rho_4);
  raise_to_t(&rho_9, &rho_8);
  mul_narrow_i(&rho_9, &rho_9, x);
  raise_to_t2(&result_t, &rho_9);
  mul_narrow_i(result, &result_t, &x_t_minus_1);
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
