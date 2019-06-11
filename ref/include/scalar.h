#ifndef SCALAR_H
#define SCALAR_H
#include <stdint.h>
#include "f11_260.h"

// Reduced to 10 limbs. For smaller tables.
typedef struct scalar {
  uint32_t limbs[9];
} scalar_t;

typedef struct scalar_hash {
  uint32_t limbs[16];
} scalar_hash_t;

// const int SCALAR_LIMBS = 9;
#define HASH_LIMBS 16
#define SCALAR_LIMBS 9
#define SCALAR_BITS 258
#define SCALAR_BYTES 33
#define SCALAR_LIMB_BITS 32
#define SCALAR_LAST_LIMB_BITS 2
#define SCALAR_LAST_LIMB_MASK 0x3

// Constants
// A scalar representing l, the order of the prime subgroup.
const scalar_t l_bits;
// For converting to SABS representation
const scalar_t signed_bits_set_adjustment;
// l * N' is congruent to -1 mod 2^32
const uint32_t SCALAR_MONT_N_PRIME;
// (2 ^ 32)^18 mod l. Used to convert to montgomery domain.
// Or to fix the result of a single multiply via a 2nd multiply.
const scalar_t SCALAR_MONT_R2;
// (2 ^ 32)^25 mod l.
// Or to fix the result of a single hash by multiply via a 2nd multiply.
const scalar_t SCALAR_MONT_R2_HASH;

// Functions for manipulating scalars. May need more for ECDSA.

// This is used to convert to SABS representation.
void divide_by_2_mod_l(scalar_t *result, const scalar_t * __restrict x);

void add_mod_l(scalar_t *result, const scalar_t * __restrict x,
               const scalar_t * __restrict y);

void sub_mod_l(scalar_t *result, const scalar_t * __restrict x,
               const scalar_t * __restrict y);

void sub_mod_l_accum(scalar_t *result, const uint32_t * __restrict x,
                     const scalar_t * __restrict y);

void mont_mult_mod_l(scalar_t *result, const scalar_t * __restrict x,
                     const scalar_t * __restrict y);

void mult_mod_l(scalar_t *result, const scalar_t * __restrict x,
                const scalar_t * __restrict y);

void mont_reduce_hash_mod_l(
  scalar_t *result, const scalar_hash_t * __restrict x);
void reduce_hash_mod_l(scalar_t *result, const scalar_hash_t * __restrict x);

void convert_to_sabs(scalar_t *result, const scalar_t * __restrict x);
#endif
