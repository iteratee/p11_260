#include <stdint.h>
#include "f11_260.h"
#include "curve.h"

#include "emmintrin.h"
#include "immintrin.h"

static inline void mask_copy_narrow(
  int32_t mask, residue_narrow_t *result,
  residue_narrow_t *x) {

  #pragma clang loop unroll(full)
  for (int i = 0; i < NLIMBS_REDUCED; ++i) {
    result->limbs[i] |= x->limbs[i] & mask;
  }
}

// Due to padding, it's one 512 bit read per coordinate.
void constant_time_extended_narrow_lookup(
  extended_pt_readd_narrow_t *result, int i, int n,
  extended_pt_readd_narrow_t *table) {

  __m512i accum[4];
  __m512i big_i = _mm512_set1_epi32(i);
  __m512i big_one = _mm512_set1_epi32(1);
  #pragma clang loop unroll(full)
  for (int j = 0; j < 4; ++j) {
    accum[j] = _mm512_setzero();
  }
  for (int j = 0; j < n; ++j) {
    __mmask16 mask = _mm512_cmpeq_epi32_mask(big_i, _mm512_setzero());
    #pragma clang loop unroll(full)
    for (int k = 0; k < 4; ++k) {
      __m512i temp = _mm512_load_si512(((__m512i*) &table[j]) + k);
      accum[k] = _mm512_mask_or_epi32(accum[k], mask, accum[k], temp);
    }
    big_i = _mm512_sub_epi32(big_i, big_one);
  }
  for (int j = 0; j < 4; ++j) {
    _mm512_store_si512(((__m512i*) result) + j, accum[j]);
  }
}

// Due to padding, it's one 512 bit read per coordinate.
void constant_time_extended_affine_narrow_lookup(
  extended_affine_pt_readd_narrow_t *result, int i, int n,
  extended_affine_pt_readd_narrow_t *table) {

  __m512i accum[3];
  __m512i big_i = _mm512_set1_epi32(i);
  __m512i big_one = _mm512_set1_epi32(1);
  #pragma clang loop unroll(full)
  for (int j = 0; j < 3; ++j) {
    accum[j] = _mm512_setzero();
  }
  for (int j = 0; j < n; ++j) {
    __mmask16 mask = _mm512_cmpeq_epi32_mask(big_i, _mm512_setzero());
    #pragma clang loop unroll(full)
    for (int k = 0; k < 3; ++k) {
      __m512i temp = _mm512_load_si512(((__m512i*) &table[j]) + k);
      accum[k] = _mm512_mask_or_epi32(accum[k], mask, accum[k], temp);
    }
    big_i = _mm512_sub_epi32(big_i, big_one);
  }
  for (int j = 0; j < 3; ++j) {
    _mm512_store_si512(((__m512i*) result) + j, accum[j]);
  }
}

void constant_time_cond_extended_negate(
  extended_pt_readd_narrow_t *x, int32_t mask32) {
  __m512i zero = _mm512_setzero();
  __m512i mask = _mm512_set1_epi32(mask32);
  __m512i not_mask = _mm512_set1_epi32(~mask32);

  #pragma clang loop unroll(full)
  for (int i = 0; i < 2; ++i) {
    __m512i temp = _mm512_load_si512(((__m512i*) x) + i);
    __m512i neg_temp = _mm512_sub_epi32(zero, temp);
    temp = _mm512_and_si512(not_mask, temp);
    neg_temp = _mm512_and_si512(mask, neg_temp);
    temp = _mm512_or_si512(temp, neg_temp);
    _mm512_store_si512(((__m512i*) x) + i, temp);
  }
}

void constant_time_cond_extended_affine_negate(
  extended_affine_pt_readd_narrow_t *x, int32_t mask32) {
  __m512i zero = _mm512_setzero();
  __m512i mask = _mm512_set1_epi32(mask32);
  __m512i not_mask = _mm512_set1_epi32(~mask32);

  #pragma clang loop unroll(full)
  for (int i = 0; i < 2; ++i) {
    __m512i temp = _mm512_load_si512(((__m512i*) x) + i);
    __m512i neg_temp = _mm512_sub_epi32(zero, temp);
    temp = _mm512_and_si512(not_mask, temp);
    neg_temp = _mm512_and_si512(mask, neg_temp);
    temp = _mm512_or_si512(temp, neg_temp);
    _mm512_store_si512(((__m512i*) x) + i, temp);
  }
}
