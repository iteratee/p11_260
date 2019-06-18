#include <stdint.h>
#include "f11_260.h"
#include "curve.h"

static inline void mask_copy_narrow(
  int32_t mask, residue_narrow_t *result,
  residue_narrow_t *x) {

  #pragma clang loop unroll(full)
  for (int i = 0; i < NLIMBS; ++i) {
    result->limbs[i] |= x->limbs[i] & mask;
  }
}

void constant_time_extended_narrow_lookup(
  extended_pt_readd_narrow_t *result, int i, int n,
  extended_pt_readd_narrow_t *table) {

  #pragma clang loop unroll(full)
  for (int i = 0; i < NLIMBS; ++i) {
    result->x.limbs[i] = 0;
    result->dt.limbs[i] = 0;
    result->y.limbs[i] = 0;
    result->z.limbs[i] = 0;
  }
  for (int j = 0; j < n; ++j) {
    int32_t mask = -(i == j);

    mask_copy_narrow(mask, &result->x, &table[j].x);
    mask_copy_narrow(mask, &result->dt, &table[j].dt);
    mask_copy_narrow(mask, &result->y, &table[j].y);
    mask_copy_narrow(mask, &result->z, &table[j].z);
  }
}

void constant_time_extended_affine_narrow_lookup(
  extended_affine_pt_readd_narrow_t *result, int i, int n,
  extended_affine_pt_readd_narrow_t *table) {

  #pragma clang loop unroll(full)
  for (int i = 0; i < NLIMBS; ++i) {
    result->x.limbs[i] = 0;
    result->dt.limbs[i] = 0;
    result->y.limbs[i] = 0;
  }

  for (int j = 0; j < n; ++j) {
    int32_t mask = -(i == j);

    mask_copy_narrow(mask, &result->x, &table[j].x);
    mask_copy_narrow(mask, &result->dt, &table[j].dt);
    mask_copy_narrow(mask, &result->y, &table[j].y);
  }
}

void constant_time_cond_extended_negate(
  extended_pt_readd_narrow_t *x, int32_t mask) {
  #pragma clang loop unroll(full)
  for (int i = 0; i < NLIMBS; ++i) {
    x->x.limbs[i] = (x->x.limbs[i] & ~mask) | ((-x->x.limbs[i]) & mask);
  }
  #pragma clang loop unroll(full)
  for (int i = 0; i < NLIMBS; ++i) {
    x->dt.limbs[i] = (x->dt.limbs[i] & ~mask) | ((-x->dt.limbs[i]) & mask);
  }
}

void constant_time_cond_extended_affine_negate(
  extended_affine_pt_readd_narrow_t *x, int32_t mask) {
  for (int i = 0; i < NLIMBS; ++i) {
    x->x.limbs[i] = (x->x.limbs[i] & ~mask) | ((-x->x.limbs[i]) & mask);
  }
  for (int i = 0; i < NLIMBS; ++i) {
    x->dt.limbs[i] = (x->dt.limbs[i] & ~mask) | ((-x->dt.limbs[i]) & mask);
  }
}
