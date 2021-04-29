#include <stdint.h>
#include "f11_260.h"

residue_narrow_t zero_narrow = {0};
residue_narrow_t one_narrow = {
  .limbs = {1},
};

// Shrink to 32 bits. Assumes reduction has already occurred, and wide storage
// is being used for vector compatibility.
void narrow(residue_narrow_t *result, const residue_wide_t * __restrict w) {
  for (int i = 0; i < NLIMBS; ++i) {
    result->limbs[i] = w->limbs[i];
  }
}

// Reduce to 10 limbs. Useful for debugging
void narrow_reduce(
  residue_narrow_reduced_t *result, const residue_narrow_t * __restrict w) {
  residue_narrow_t temp;
  for (int i = 0; i < NLIMBS; ++i) {
    temp.limbs[i] = w->limbs[i] - w->limbs[10];
  }

  reduce_step_narrow(&temp, &temp);

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
  for (int i = 0; i < NLIMBS - 1; ++i) {
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

// Subtract 2 12x32-bit residues.
void sub_narrow(
  residue_narrow_t *result, const residue_narrow_t * __restrict x,
  const residue_narrow_t * __restrict y) {

  for (int i = 0; i < NLIMBS; ++i) {
    result->limbs[i] = x->limbs[i] - y->limbs[i];
  }
}

// negate a 12x32-bit residue.
void negate_narrow(
  residue_narrow_t *result, const residue_narrow_t *x) {

  for (int i = 0; i < NLIMBS; ++i) {
    result->limbs[i] = -(x->limbs[i]);
  }
}

// Add 2 12x32-bit residues.
void add_narrow(
  residue_narrow_t *result, const residue_narrow_t * __restrict x,
  const residue_narrow_t * __restrict y) {

  for (int i = 0; i < NLIMBS; ++i) {
    result->limbs[i] = x->limbs[i] + y->limbs[i];
  }
}

// Scale a narrow residue by 2.
void double_narrow(
  residue_narrow_t *result, const residue_narrow_t *x) {

  for (int i = 0; i < NLIMBS; ++i) {
    result->limbs[i] = x->limbs[i] << 1;
  }
}

// Scale a wide residue by 2.
void double_wide(
  residue_wide_t *result, const residue_wide_t *x) {

  for (int i = 0; i < NLIMBS; ++i) {
    result->limbs[i] = x->limbs[i] << 1;
  }
}

#define wrap(x) (((x + NLIMBS) % NLIMBS))
// Multiply two wide residues, and produce a wide result. The result is reduced
// to 32 bits, but not narrowed for performance reasons.
void mul_wide(
  residue_wide_t *result, const residue_wide_t *x, const residue_wide_t *y) {

  residue_wide_t temp;
  for (int i = 0; i < NLIMBS; ++i) {
    temp.limbs[i] = 0;
    int i_2 = (i + (-(i & 1) & NLIMBS)) >> 1;
    for (int j = 1; j <= NLIMBS / 2; ++ j) {
      temp.limbs[i] +=
        (x->limbs[wrap(i_2 + j)] - x->limbs[wrap(i_2 - j)]) *
        (y->limbs[wrap(i_2 - j)] - y->limbs[wrap(i_2 + j)]);
    }
  }
  reduce_step_wide(&temp, &temp);
  reduce_step_wide(result, &temp);
}

// Multiply a wide residues by a narrow and produce a wide result. The result is
// reduced to 32 bits, but not narrowed for performance reasons.
void mul_wide_narrow(
  residue_wide_t *result, const residue_wide_t *x, const residue_narrow_t *y) {

  residue_wide_t temp;
  for (int i = 0; i < NLIMBS; ++i) {
    temp.limbs[i] = 0;
    int i_2 = (i + (-(i & 1) & NLIMBS)) >> 1;
    for (int j = 1; j <= NLIMBS / 2; ++j) {
      temp.limbs[i] +=
        (x->limbs[wrap(i_2 + j)] - x->limbs[wrap(i_2 - j)]) *
        ((int64_t) (y->limbs[wrap(i_2 - j)] - y->limbs[wrap(i_2 + j)]));
    }
  }
  reduce_step_wide(&temp, &temp);
  reduce_step_wide(result, &temp);
}

// Multiply two narrow residues and produce a narrow result.
void mul_narrow(
  residue_narrow_t *result, const residue_narrow_t *x,
  const residue_narrow_t *y) {

  residue_wide_t temp;
  for (int i = 0; i < NLIMBS; ++i) {
    temp.limbs[i] = 0;
    int i_2 = (i + (-(i & 1) & NLIMBS)) >> 1;
    for (int j = 1; j <= NLIMBS / 2; ++ j) {
      temp.limbs[i] +=
        ((int64_t) (x->limbs[wrap(i_2 + j)] - x->limbs[wrap(i_2 - j)])) *
        ((int64_t) (y->limbs[wrap(i_2 - j)] - y->limbs[wrap(i_2 + j)]));
    }
  }
  reduce_step_wide(&temp, &temp);
  reduce_step_wide(&temp, &temp);
  narrow(result, &temp);
}

// Multiply a narrow residue by a small constant. The result is reduced to 32
// bits, but not narrowed for performance reasons.
void mul_narrow_const(
  residue_narrow_t *result, const residue_narrow_t *x, int32_t d) {

  residue_wide_t temp;
  for (int i = 0; i < NLIMBS; ++i) {
    temp.limbs[i] = ((uint64_t) x->limbs[i]) * d;
  }
  reduce_step_wide(&temp, &temp);
  narrow(result, &temp);
}


// Square a narrow residue and produce a wide result. The result is reduced to
// 32 bits but not narrowed for performance reasons.
void square_narrow(
  residue_narrow_t *result, const residue_narrow_t *x) {

  residue_wide_t temp;
  for (int i = 0; i < NLIMBS; ++i) {
    temp.limbs[i] = 0;
    int i_2 = (i + (-(i & 1) & NLIMBS)) >> 1;
    for (int j = 1; j <= NLIMBS / 2; ++ j) {
      temp.limbs[i] -=
        ((int64_t) (x->limbs[wrap(i_2 + j)] - x->limbs[wrap(i_2 - j)])) *
        ((int64_t) (x->limbs[wrap(i_2 + j)] - x->limbs[wrap(i_2 - j)]));
    }
  }
  reduce_step_wide(&temp, &temp);
  reduce_step_wide(&temp, &temp);
  narrow(result, &temp);
}

// Approximately divide each coefficient by t. Carry the results.
void reduce_step_narrow(
  residue_narrow_t *result, const residue_narrow_t *x) {

  int32_t carries[NLIMBS];

  for (int i = 0; i < NLIMBS; ++i) {
    carries[i] = x->limbs[i] >> TBITS;
    result->limbs[i] = (x->limbs[i] & TMASK) +
      (carries[i] << T_CBITS) - carries[i];
  }

  for (int i = 1; i < NLIMBS; ++i) {
    result->limbs[i] += carries[i - 1];
  }
  result->limbs[0] += carries[NLIMBS - 1];
}

// Approximately divide each coefficient by t. Carry the results.
void reduce_step_wide(
  residue_wide_t *result, const residue_wide_t *x) {

  int64_t carries[NLIMBS];

  for (int i = 0; i < NLIMBS; ++i) {
    carries[i] = x->limbs[i] >> TBITS;
    result->limbs[i] = (x->limbs[i] & TMASK) +
      (carries[i] << T_CBITS) - carries[i];
  }

  for (int i = 1; i < NLIMBS; ++i) {
    result->limbs[i] += carries[i - 1];
  }
  result->limbs[0] += carries[NLIMBS - 1];
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

  square_narrow(result, x);
  for (int i = 1; i < n; ++i) {
    square_narrow(result, result);
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

  square_narrow(&z2, x);
  mul_narrow(&z2, &z2, x);
  square_narrow(&z3, &z2);
  mul_narrow(&z3, &z3, x);
  nsquare_narrow(&z5, &z3, 2);
  mul_narrow(&z5, &z5, &z2);
  nsquare_narrow(&z10, &z5, 5);
  mul_narrow(&z10, &z10, &z5);
  square_narrow(&z11, &z10);
  mul_narrow(&z11, &z11, x);
  nsquare_narrow(&z22, &z11, 11);
  mul_narrow(&z22, &z22, &z11);
  nsquare_narrow(&result_t, &z22, 4);
  mul_narrow(result, &result_t, x);
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

  square_narrow(&z2, x);
  mul_narrow(&z2, &z2, x);
  square_narrow(&z3, &z2);
  mul_narrow(&z3, &z3, x);
  nsquare_narrow(&z5, &z3, 2);
  mul_narrow(&z5, &z5, &z2);
  nsquare_narrow(&z10, &z5, 5);
  mul_narrow(&z10, &z10, &z5);
  square_narrow(&z11, &z10);
  mul_narrow(&z11, &z11, x);
  nsquare_narrow(&z22, &z11, 11);
  mul_narrow(&z22, &z22, &z11);
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
  residue_narrow_t x_t;
  residue_narrow_t phi_8_x_t;
  residue_narrow_t phi_8_x_t_t;

  raise_to_t_minus_1_over_4(&x_t_minus_1_over_4, x);
  nsquare_narrow(&x_t_minus_1, &x_t_minus_1_over_4, 2);
  mul_narrow(&x_t, &x_t_minus_1, x);
  raise_to_phi_t(&phi_8_x_t, &x_t, 8);
  raise_to_t(&phi_8_x_t_t, &phi_8_x_t);
  mul_narrow(result, &phi_8_x_t_t, &x_t_minus_1);
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
