#define _DEFAULT_SOURCE
#include <string.h>
#include <stdint.h>
#include "f11_260.h"
#include "scalar.h"

// Plenty of inspiration for this file was taken from Mike Hamburg's
// Ed448 code.

// Constants:
__attribute__((__aligned__(32)))
const scalar_t l_bits = {
  .limbs = {0x28ad9c41, 0xe6dcf7e8, 0x34b804af, 0x5af91169,
            0x5cf68f2f, 0x125277f4, 0x9c1bf9f, 0xffff6b00, 0x3,},
};

__attribute__((__aligned__(32)))
const scalar_t signed_bits_set_adjustment = {
  .limbs = {0x5d498efb, 0x648c205f, 0x2d1fed40, 0x941bba5b,
            0x8c25c342, 0xb6b6202e, 0xd8f90183, 0x000253ff, 0x0,},
};

__attribute__((__aligned__(32)))
const scalar_t SCALAR_MONT_R2 = {
  .limbs = {0x30ba45c7, 0xf3422093, 0x054bbbf6, 0x017ab264,
            0x914ee18b, 0x250f1097, 0xf6bc1224, 0x5e97c70e, 0x2,},
};

const uint32_t SCALAR_MONT_N_PRIME = 0xb3138c3f;

__attribute__((__aligned__(32)))
const scalar_t SCALAR_MONT_R2_HASH = {
  .limbs = {
            0x202dd8e7, 0xcb1bf7be, 0xd219daf6, 0xb85aba0a,
            0xdc8da05f, 0xbd23bfce, 0xb7642c95, 0xbb13e4ad, 0x0,},
};

void divide_by_2_mod_l(
  scalar_t *result, const scalar_t *x) {

  uint32_t mask = -(x->limbs[0] & 1);

  uint64_t chain = 0;
  for (int i = 0; i < SCALAR_LIMBS; ++i) {
    chain = (chain + x->limbs[i]) + (mask & l_bits.limbs[i]);
    result->limbs[i] = chain;
    chain >>= SCALAR_LIMB_BITS;
  }

  int i;
  for (i = 0; i < SCALAR_LIMBS - 1; ++i) {
    result->limbs[i] = result->limbs[i] >> 1 |
      (result->limbs[i+1] << (SCALAR_LIMB_BITS - 1));
  }
  result->limbs[i] >>= 1;
}

void add_mod_l(
  scalar_t *result, const scalar_t *x,
  const scalar_t * __restrict y) {

  uint64_t chain = 0;
  int i;
  for (i = 0; i < SCALAR_LIMBS; ++i) {
    chain = (chain + x->limbs[i]) + y->limbs[i];
    result->limbs[i] = chain;
    chain >>= SCALAR_LIMB_BITS;
  }

  sub_mod_l(result, result, &l_bits);
}

void sub_mod_l(
  scalar_t *result, const scalar_t *x,
  const scalar_t *y) {
  sub_mod_l_accum(result, x->limbs, y);
}

// x is a pointer and not a scalar_t so that this function can be used to reduce
// accumulators after multiplication.
void sub_mod_l_accum(
  scalar_t *result, const uint32_t *x,
  const scalar_t *y) {

  int64_t chain = 0;
  int i;
  for (i = 0; i < SCALAR_LIMBS; ++i) {
    chain = (chain + x[i]) - y->limbs[i];
    result->limbs[i] = chain;
    chain >>= SCALAR_LIMB_BITS;
  }

  //Should be 0 or -1 (to function as a mask)
  int32_t borrow = chain;

  chain = 0;
  for (i = 0; i < SCALAR_LIMBS; ++i) {
    chain = (chain + result->limbs[i]) + (l_bits.limbs[i] & borrow);
    result->limbs[i] = chain;
    chain >>= SCALAR_LIMB_BITS;
  }
}

void convert_to_sabs(
  scalar_t *result, const scalar_t *x) {
  add_mod_l(result, x, &signed_bits_set_adjustment);
  divide_by_2_mod_l(result, result);
}

#include <stdio.h>
void mont_reduce_hash_mod_l(
  scalar_t *result, const scalar_hash_t * __restrict x) {
  uint32_t accum[HASH_LIMBS];

  for (int i = 0; i < HASH_LIMBS; ++i) {
    accum[i] = x->limbs[i];
  }

  uint64_t chain = 0;
  for (int i = 0; i <= HASH_LIMBS - SCALAR_LIMBS; ++i) {
    uint32_t q = accum[0] * SCALAR_MONT_N_PRIME;
    for (int j = 0; j < SCALAR_LIMBS; ++j) {
      chain += accum[j] + ((uint64_t) q) * l_bits.limbs[j];
      if (j > 0) {
        accum[j - 1] = chain;
      }
      chain >>= SCALAR_LIMB_BITS;
    }
    int j;
    for (j = SCALAR_LIMBS; j < HASH_LIMBS - i; ++j) {
      chain += accum[j];
      accum[j - 1] = chain;
      chain >>= SCALAR_LIMB_BITS;
    }
    accum[j - 1] = chain;
  }

  for (int i = 0; i < SCALAR_LIMBS; ++i) {
    result->limbs[i] = accum[i];
  }
  explicit_bzero(accum, sizeof(accum));
}

void reduce_hash_mod_l(scalar_t *result, const scalar_hash_t * __restrict x) {
  mont_reduce_hash_mod_l(result, x);
  mont_mult_mod_l(result, result, &SCALAR_MONT_R2_HASH);
}

void mont_mult_mod_l(scalar_t *result, const scalar_t *x,
                     const scalar_t *y) {
  uint32_t accum[SCALAR_LIMBS + 1] = {0};

  for (int i = 0; i < SCALAR_LIMBS; ++i) {
    uint32_t x_limb = x->limbs[i];

    uint64_t chain = 0;
    int j;
    for (j = 0; j < SCALAR_LIMBS; ++j) {
      chain += accum[j] + ((uint64_t) y->limbs[j]) * x_limb;
      accum[j] = chain;
      chain >>= SCALAR_LIMB_BITS;
    }

    // 2 bit value
    accum[j] = chain;

    uint32_t q = accum[0] * SCALAR_MONT_N_PRIME;
    chain = 0;
    for (int j = 0; j < SCALAR_LIMBS; ++j) {
      chain += accum[j] + ((uint64_t) l_bits.limbs[j]) * q;
      if (j > 0) {
        accum[j - 1] = chain;
      }
      chain >>= SCALAR_LIMB_BITS;
    }

    // chain is a 2-bit value with a possible carry.
    // result is a 3 bit value
    chain += accum[j];
    accum[j - 1] = chain;
  }

  sub_mod_l_accum(result, accum, &l_bits);
  explicit_bzero(accum, sizeof(accum));
}

void mult_mod_l(scalar_t *result, const scalar_t * __restrict x,
                const scalar_t * __restrict y) {
  scalar_t temp;
  mont_mult_mod_l(&temp, x, y);
  mont_mult_mod_l(result, &temp, &SCALAR_MONT_R2);
  explicit_bzero(&temp, sizeof(temp));
}
