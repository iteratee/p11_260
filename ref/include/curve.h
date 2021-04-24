#ifndef CURVE_H
#define CURVE_H
#include "f11_260.h"
#include "scalar.h"

typedef struct affine_pt_narrow {
  residue_narrow_t x;
  residue_narrow_t y;
} affine_pt_narrow_t;

typedef struct extended_pt_readd_narrow {
  __attribute__((__aligned__(64)))
  residue_narrow_t x;
  residue_narrow_t dt;
  residue_narrow_t y;
  residue_narrow_t z;
} extended_pt_readd_narrow_t;

typedef struct extended_affine_pt_readd_narrow {
  __attribute__((__aligned__(64)))
  residue_narrow_t x;
  residue_narrow_t dt;
  residue_narrow_t y;
} extended_affine_pt_readd_narrow_t;

// For use in doubling.
typedef struct projective_pt_narrow {
  residue_narrow_t x;
  residue_narrow_t y;
  residue_narrow_t z;
} projective_pt_narrow_t;

// For use in addition.
typedef struct extended_pt_narrow {
  residue_narrow_t x;
  residue_narrow_t y;
  residue_narrow_t t;
  residue_narrow_t z;
} extended_pt_narrow_t;

#define D (-49142)

__attribute__((__aligned__(32)))
const affine_pt_narrow_t B;

void copy_projective_pt_narrow(
  projective_pt_narrow_t *result, const projective_pt_narrow_t *source);

void copy_extended_pt_narrow(
  extended_pt_narrow_t *result, const extended_pt_narrow_t *source);

void copy_extended_pt_readd_narrow(
  extended_pt_readd_narrow_t *result, const extended_pt_readd_narrow_t *source);

void copy_extended_pt_readd_narrow(
  extended_pt_readd_narrow_t *result, const extended_pt_readd_narrow_t *source);

void copy_extended_affine_pt_readd_narrow(
  extended_affine_pt_readd_narrow_t *result,
  const extended_affine_pt_readd_narrow_t *source);

void negate_extended_pt_readd_narrow(
  extended_pt_readd_narrow_t *result,
  const extended_pt_readd_narrow_t *source);

void negate_extended_affine_pt_readd_narrow(
  extended_affine_pt_readd_narrow_t *result,
  const extended_affine_pt_readd_narrow_t *source);

void affine_narrow_to_extended(
  extended_pt_narrow_t *result,
  const affine_pt_narrow_t * __restrict x);

void affine_to_projective(
  projective_pt_narrow_t *result,
  const affine_pt_narrow_t * __restrict x);

void affine_to_readd_narrow(
  extended_pt_readd_narrow_t *result,
  const affine_pt_narrow_t * __restrict x);

void extended_to_readd_narrow_neg(
  extended_pt_readd_narrow_t *result,
  const extended_pt_narrow_t * __restrict x);

void affine_to_readd_narrow(
  extended_pt_readd_narrow_t *result,
  const affine_pt_narrow_t * __restrict x);

void projective_to_extended_narrow(
  extended_pt_narrow_t *result, projective_pt_narrow_t * __restrict x);

void extended_to_projective_narrow(
  projective_pt_narrow_t *result, const extended_pt_narrow_t * __restrict x);

void readd_to_projective(
  projective_pt_narrow_t *result,
  const extended_pt_readd_narrow_t * __restrict x);

void affine_readd_to_extended(
  extended_pt_narrow_t *result,
  const extended_affine_pt_readd_narrow_t * __restrict x);

void negate_extended_affine_pt_readd_narrow(
  extended_affine_pt_readd_narrow_t *result,
  const extended_affine_pt_readd_narrow_t *source);

void affine_double(
  projective_pt_narrow_t *result,
  const affine_pt_narrow_t * __restrict x);

void affine_double_extended(
  extended_pt_narrow_t *result, const affine_pt_narrow_t * __restrict x);

void projective_double(
  projective_pt_narrow_t *result, const projective_pt_narrow_t *x);

void projective_double_extended(
  extended_pt_narrow_t *result, const projective_pt_narrow_t * __restrict x);

void extended_double_extended(
  extended_pt_narrow_t *result, const extended_pt_narrow_t *x);

void projective_add(
  projective_pt_narrow_t *result, const projective_pt_narrow_t * __restrict x1,
  const projective_pt_narrow_t * __restrict x2);

void extended_add(
  projective_pt_narrow_t *result, const extended_pt_narrow_t * __restrict x,
  const extended_pt_narrow_t * __restrict y);

void extended_add_extended(
  extended_pt_narrow_t *result, const extended_pt_narrow_t * __restrict x,
  const extended_pt_narrow_t * __restrict y);

void extended_readd_narrow(
  projective_pt_narrow_t *result, const extended_pt_narrow_t * __restrict x,
  const extended_pt_readd_narrow_t * __restrict y);

void extended_readd_narrow_extended(
  extended_pt_narrow_t *result, const extended_pt_narrow_t * __restrict x,
  const extended_pt_readd_narrow_t * __restrict y);

void extended_readd_affine_narrow_extended(
  extended_pt_narrow_t *result, const extended_pt_narrow_t * __restrict x,
  const extended_affine_pt_readd_narrow_t * __restrict y);

void extended_add_extended(
  extended_pt_narrow_t *result, const extended_pt_narrow_t * __restrict x,
  const extended_pt_narrow_t * __restrict y);

void extended_readd_readd_narrow(
  extended_pt_readd_narrow_t *result,
  const extended_pt_narrow_t * __restrict x,
  const extended_pt_readd_narrow_t * __restrict y);

void extended_readd_narrow_extended(
  extended_pt_narrow_t *result,
  const extended_pt_narrow_t *x1,
  const extended_pt_readd_narrow_t * __restrict x2);

void scalar_multiply(
  projective_pt_narrow_t *result, const affine_pt_narrow_t * __restrict x,
  const scalar_t * __restrict n);

void scalar_multiply_unsafe(
  projective_pt_narrow_t *result, const affine_pt_narrow_t * __restrict x,
  const scalar_t * __restrict n);

int point_decompress(
  affine_pt_narrow_t *result, residue_narrow_reduced_t *y, int low_bit);
#endif
