#include <stdint.h>
#include "f11_260.h"
#include "curve.h"

#include <stdio.h>
static void print_residue_narrow(const residue_narrow_t *x) {
  printf("[");
  for (int i = 0; i < NLIMBS; ++i) {
    printf(" %#x,", x->limbs[i]);
    // printf("x[%d]: %d\n", i, x[i]);
  }
  printf(" ]");
  printf("\n");
}

static inline void mask_copy_narrow_reduced(
  int32_t mask, residue_narrow_t *result,
  residue_narrow_reduced_t *x) {

  #pragma clang loop unroll(full)
  for (int i = 0; i < NLIMBS_REDUCED; ++i) {
    result->limbs[i+1] |= x->limbs[i] & mask;
  }
}

void constant_time_extended_narrow_reduced_lookup(
  extended_pt_readd_narrow_t *result, int i, int n,
  extended_pt_readd_narrow_reduced_t *table) {

  #pragma clang loop unroll(full)
  for (int i = 0; i < NLIMBS; ++i) {
    result->x.limbs[i] = 0;
    result->dt.limbs[i] = 0;
    result->y.limbs[i] = 0;
    result->z.limbs[i] = 0;
  }
  for (int j = 0; j < n; ++j) {
    int32_t mask = -(i == j);

    mask_copy_narrow_reduced(mask, &result->x, &table[j].x);
    mask_copy_narrow_reduced(mask, &result->dt, &table[j].dt);
    mask_copy_narrow_reduced(mask, &result->y, &table[j].y);
    mask_copy_narrow_reduced(mask, &result->z, &table[j].z);
  }
}

void constant_time_extended_affine_narrow_reduced_lookup(
  extended_affine_pt_readd_narrow_t *result, int i, int n,
  extended_affine_pt_readd_narrow_reduced_t *table) {

  #pragma clang loop unroll(full)
  for (int i = 0; i < NLIMBS; ++i) {
    result->x.limbs[i] = 0;
    result->dt.limbs[i] = 0;
    result->y.limbs[i] = 0;
  }

  for (int j = 0; j < n; ++j) {
    int32_t mask = -(i == j);

    mask_copy_narrow_reduced(mask, &result->x, &table[j].x);
    mask_copy_narrow_reduced(mask, &result->dt, &table[j].dt);
    mask_copy_narrow_reduced(mask, &result->y, &table[j].y);
  }
}

void constant_time_cond_extended_negate(
  extended_pt_readd_narrow_t *x, int32_t mask) {
  #pragma clang loop unroll(full)
  for (int i = 0; i < NLIMBS_REDUCED; ++i) {
    x->x.limbs[i + 1] = (x->x.limbs[i + 1] & ~mask) | ((-x->x.limbs[i+1]) & mask);
  }
  #pragma clang loop unroll(full)
  for (int i = 0; i < NLIMBS_REDUCED; ++i) {
    x->dt.limbs[i + 1] = (x->dt.limbs[i + 1] & ~mask) | ((-x->dt.limbs[i+1]) & mask);
  }
}

void constant_time_cond_extended_affine_negate(
  extended_affine_pt_readd_narrow_t *x, int32_t mask) {
  for (int i = 0; i < NLIMBS_REDUCED; ++i) {
    x->x.limbs[i + 1] = (x->x.limbs[i + 1] & ~mask) | ((-x->x.limbs[i+1]) & mask);
  }
  for (int i = 0; i < NLIMBS_REDUCED; ++i) {
    x->dt.limbs[i + 1] = (x->dt.limbs[i + 1] & ~mask) | ((-x->dt.limbs[i+1]) & mask);
  }
}
