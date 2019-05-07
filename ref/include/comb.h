#ifndef COMB_H
#define COMB_H

#include "curve.h"
#include "scalar.h"

#define COMB_TABLE_SIZE 16
#define COMB_TEETH 5
#define COMB_COUNT 4
#define COMB_SEPARATION 13
#define COMB_LOOKUP_MASK 0xf

// A single comb table.
typedef struct sabs_single_comb {
  extended_affine_pt_readd_narrow_reduced_t table[COMB_TABLE_SIZE];
} sabs_single_comb_t;

// A single wide comb table. Used in computing a narrow comb table.
typedef struct sabs_single_comb_wide {
  projective_pt_wide_t table[COMB_TABLE_SIZE];
} sabs_single_comb_wide_t;

// A comb set. There is a precomputed comb set for the base point, but for
// verifications of several signatures from the same key, it would be
// advantageous to precompute a comb.
typedef struct sabs_comb_set {
  sabs_single_comb_t combs[COMB_COUNT];
} sabs_comb_set_t;

// An unreduced comb set. Used just to separate the logic of comb computation
// from comb reduction.
typedef struct sabs_comb_set_wide {
  sabs_single_comb_wide_t combs[COMB_COUNT];
} sabs_comb_set_wide_t;

// used for computing the entries in the comb table.
typedef struct teeth_set {
  extended_pt_readd_wide_t teeth[COMB_TEETH - 1];
} teeth_set_t;

sabs_comb_set_t base_comb;

void compute_comb_set(
  sabs_comb_set_t *result, const affine_pt_narrow_reduced_t *base_pt);

void reduce_comb_set(sabs_comb_set_t *result, sabs_comb_set_wide_t *source);

void scalar_comb_multiply(
  projective_pt_wide_t *result, const sabs_comb_set_t * __restrict comb,
  const scalar_t * __restrict n);
#endif
