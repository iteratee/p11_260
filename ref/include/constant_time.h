#ifndef CONSTANT_TIME_H
#define CONSTANT_TIME_H
#include <stdint.h>
#include "f11_260.h"
#include "curve.h"

void constant_time_extended_narrow_reduced_lookup(
  extended_pt_readd_narrow_t *result, int i, int n,
  const extended_pt_readd_narrow_reduced_t *table);

void constant_time_extended_affine_narrow_reduced_lookup(
  extended_affine_pt_readd_narrow_t *result, int i, int n,
  const extended_affine_pt_readd_narrow_reduced_t *table);

void constant_time_cond_extended_negate(
  extended_pt_readd_narrow_t *x, int32_t mask);

void constant_time_cond_extended_affine_negate(
  extended_affine_pt_readd_narrow_t *x, int32_t mask);
#endif
