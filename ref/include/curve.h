#ifndef CURVE_H
#define CURVE_H
#include "f11_260.h"
#include "scalar.h"

typedef struct affine_pt_narrow_reduced {
  residue_narrow_reduced_t x;
  residue_narrow_reduced_t y;
} affine_pt_narrow_reduced_t;

// For use in single use windows
typedef struct extended_pt_readd_narrow_reduced {
  // x and dt are negated together, so they are contiguous in the structure
  residue_narrow_reduced_t x;
  residue_narrow_reduced_t dt;
  residue_narrow_reduced_t y;
  residue_narrow_reduced_t z;
} extended_pt_readd_narrow_reduced_t;

typedef struct extended_pt_readd_narrow {
  residue_narrow_t x;
  residue_narrow_t dt;
  residue_narrow_t y;
  residue_narrow_t z;
} extended_pt_readd_narrow_t;

typedef struct extended_pt_readd_wide {
  residue_wide_t x;
  residue_wide_t dt;
  residue_wide_t y;
  residue_wide_t z;
} extended_pt_readd_wide_t;

// For use in pre-computed tables.
typedef struct extended_affine_pt_readd_narrow_reduced {
  residue_narrow_reduced_t x;
  residue_narrow_reduced_t dt;
  residue_narrow_reduced_t y;
  // This makes the structure 32 dwords, much nicer for using the vector unit
  // for a constant time table load.
  uint32_t pad[2];
} extended_affine_pt_readd_narrow_reduced_t;

typedef struct extended_affine_pt_readd_narrow {
  residue_narrow_t x;
  residue_narrow_t dt;
  residue_narrow_t y;
} extended_affine_pt_readd_narrow_t;

// For use in doubling.
typedef struct projective_pt_wide {
  residue_wide_t x;
  residue_wide_t y;
  residue_wide_t z;
} projective_pt_wide_t;

// For use in addition.
typedef struct extended_pt_wide {
  residue_wide_t x;
  residue_wide_t y;
  residue_wide_t t;
  residue_wide_t z;
} extended_pt_wide_t;

#define D (-49142)

const affine_pt_narrow_reduced_t B;

void copy_projective_pt_wide(
  projective_pt_wide_t *result, const projective_pt_wide_t *source);

void copy_extended_pt_wide(
  extended_pt_wide_t *result, const extended_pt_wide_t *source);

void copy_extended_pt_readd_wide(
  extended_pt_readd_wide_t *result, const extended_pt_readd_wide_t *source);

void negate_extended_pt_readd_wide(
  extended_pt_readd_wide_t *result,
  const extended_pt_readd_wide_t *source);

void affine_narrow_reduced_to_extended(
  extended_pt_wide_t *result,
  const affine_pt_narrow_reduced_t * __restrict x);

void affine_to_projective(
  projective_pt_wide_t *result,
  const affine_pt_narrow_reduced_t * __restrict x);

void affine_to_readd_wide(
  extended_pt_readd_wide_t *result,
  const affine_pt_narrow_reduced_t * __restrict x);

void extended_to_readd_wide_neg(
  extended_pt_readd_wide_t *result,
  const extended_pt_wide_t * __restrict x);

void affine_to_readd_narrow_reduced(
  extended_pt_readd_narrow_reduced_t *result,
  const affine_pt_narrow_reduced_t * __restrict x);

void projective_to_readd_wide(
  extended_pt_readd_wide_t *result, projective_pt_wide_t * __restrict x);

void extended_to_projective_wide(
  projective_pt_wide_t *result, const extended_pt_wide_t * __restrict x);

void readd_to_projective(
  projective_pt_wide_t *result,
  const extended_pt_readd_narrow_t * __restrict x);

void affine_readd_to_extended(
  extended_pt_wide_t *result,
  const extended_affine_pt_readd_narrow_t * __restrict x);

void affine_double(
  projective_pt_wide_t *result,
  const affine_pt_narrow_reduced_t * __restrict x);

void affine_double_extended(
  extended_pt_wide_t *result, const affine_pt_narrow_reduced_t * __restrict x);

void projective_double(
  projective_pt_wide_t *result, const projective_pt_wide_t *x);

void projective_double_extended(
  extended_pt_wide_t *result, const projective_pt_wide_t * __restrict x);

void extended_double_extended(
  extended_pt_wide_t *result, const extended_pt_wide_t *x);

void extended_add(
  projective_pt_wide_t *result, const extended_pt_wide_t * __restrict x,
  const extended_pt_wide_t * __restrict y);

void extended_add_extended(
  extended_pt_wide_t *result, const extended_pt_wide_t * __restrict x,
  const extended_pt_wide_t * __restrict y);

void extended_readd_narrow(
  projective_pt_wide_t *result, const extended_pt_wide_t * __restrict x,
  const extended_pt_readd_narrow_t * __restrict y);

void extended_readd_narrow_extended(
  extended_pt_wide_t *result, const extended_pt_wide_t * __restrict x,
  const extended_pt_readd_narrow_t * __restrict y);

void extended_readd_affine_narrow_extended(
  extended_pt_wide_t *result, const extended_pt_wide_t * __restrict x,
  const extended_affine_pt_readd_narrow_t * __restrict y);

void extended_add_extended(
  extended_pt_wide_t *result, const extended_pt_wide_t * __restrict x,
  const extended_pt_wide_t * __restrict y);

void extended_readd_readd_narrow_reduce(
  extended_pt_readd_narrow_reduced_t *result,
  const extended_pt_wide_t * __restrict x,
  const extended_pt_readd_narrow_reduced_t * __restrict y);

void extended_readd_wide_extended(
  extended_pt_wide_t *result,
  const extended_pt_wide_t *x1,
  const extended_pt_readd_wide_t * __restrict x2);

void scalar_multiply(
  projective_pt_wide_t *result, const affine_pt_narrow_reduced_t * __restrict x,
  const scalar_t * __restrict n);
#endif
