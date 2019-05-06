// Types and functions for manipulating field elements

#ifndef F11_260_H
#define F11_260_H
#include <stdint.h>

#define NLIMBS_REDUCED 10
#define NLIMBS 12
#define T ((1 << 26) - 15)
#define TBITS 26
#define TMASK ((1 << 26) - 1)
#define T_CBITS 4

// Reduced to 10 limbs. For smaller tables.
typedef struct residue_narrow_reduced {
  int32_t limbs[10];
} residue_narrow_reduced_t;

// 11 limbs. Limb 10 is placed in slot 0, and slot 11.
typedef struct residue_narrow {
  int32_t limbs[12];
} residue_narrow_t;

// 11 limbs. Limb 10 is placed in slot 0 and slot 11. Wider for vector
// compatibility.
typedef struct residue_wide {
  int64_t limbs[12];
} residue_wide_t;

// Shrink to 32 bits. Assumes reduction has already occurred, and wide storage
// is being used for vector compatibility.
void narrow(residue_narrow_t *result, const residue_wide_t * __restrict w);

// Reduce to 10 limbs. Suitable for table entries.
void narrow_reduce(
  residue_narrow_reduced_t *result, const residue_narrow_t * __restrict w);

// Reduce to unique representative.
// This is expensive. Only used for final signature or DH Key
void narrow_complete(
  residue_narrow_reduced_t *result, const residue_narrow_t * __restrict w);

// Reduce to mostly unique representative.
// All coefficients are reduced to 0 <= xi <= t
// Unique up to carries (xi == t) => (xi = 0; x[i+1] += 1);
// This is sufficient to determine if x is even or odd.
// Still pretty expensive. Used in point compression.
void narrow_partial_complete(
  residue_narrow_reduced_t *result, const residue_narrow_t * __restrict w);

// Produce a 32-bit entry with 11 limbs
void unnarrow_reduce(
  residue_narrow_t *result, const residue_narrow_reduced_t * __restrict x);

// Produce a 64-bit residue
void widen(
  residue_wide_t *result, const residue_narrow_t * __restrict x);

// Copy a 64-bit residue
void copy_wide(
  residue_wide_t *result, const residue_wide_t * __restrict x);

// Subtract 2 12x64-bit residues.
void sub_wide(
  residue_wide_t *result, const residue_wide_t * __restrict x,
  const residue_wide_t * __restrict y);

// Add 2 12x32-bit residues.
void add_narrow(
  residue_narrow_t *result, const residue_narrow_t * __restrict x,
  const residue_narrow_t * __restrict y);

// Add 2 12x64-bit residues.
void add_wide(
  residue_wide_t *result, const residue_wide_t * __restrict x,
  const residue_wide_t * __restrict y);

// Scale a wide residue by 2.
void double_wide(
  residue_wide_t *result, const residue_wide_t * __restrict x);

// Multiply two wide residues, and produce a wide result. The result is reduced
// to 32 bits, but not narrowed for performance reasons.
void mul_wide(
  residue_wide_t *result, const residue_wide_t * __restrict x,
  const residue_wide_t * __restrict y);
// Multiply a wide residues by a narrow and produce a wide result. The result is
// reduced to 32 bits, but not narrowed for performance reasons.
void mul_wide_narrow(
  residue_wide_t *result, const residue_wide_t * __restrict x,
  const residue_narrow_t * __restrict y);
// Multiply two narrow residues and produce a wide result. The result is reduced
// to 32 bits, but not narrowed for performance reasons.
void mul_narrow(
  residue_wide_t *result, const residue_narrow_t * __restrict x,
  const residue_narrow_t * __restrict y);

// Multiply a wide residue by a constant.
void mul_wide_const(
  residue_wide_t *result, const residue_wide_t * __restrict x, int32_t d);

// Multiply a narrow residue by a constant, producing a wide result
void mul_narrow_const(
  residue_wide_t *result, const residue_narrow_t *x, int32_t d);

// Square a wide residue and produce a wide result. The result is reduced to 32
// bits but not narrowed for performance reasons.
void square_wide(
  residue_wide_t *result, const residue_wide_t * __restrict x);

// Square a narrow residue and produce a wide result. The result is reduced to
// 32 bits but not narrowed for performance reasons.
void square_narrow(
  residue_wide_t *result, const residue_narrow_t * __restrict x);

// Approximately divide each coefficient by t. Carry the results.
void reduce_step_narrow(
  residue_narrow_t *result, const residue_narrow_t *x);

// Approximately divide each coefficient by t. Carry the results.
void reduce_step_wide(
  residue_wide_t *result, const residue_wide_t *x);

// Invert via fermat's theorem
void invert_wide(
  residue_wide_t *result, const residue_wide_t * __restrict x);

// Compute combined inverse and square root
// returns true if x/y was a quadratic residue, and false otherwise.
int sqrt_inv_wide(
  residue_wide_t *result, const residue_wide_t * __restrict x,
  const residue_wide_t * __restrict y);

// Returns true if x == y. Computes in constant time.
int equal_wide(const residue_wide_t * x, const residue_wide_t * y);

void encode(uint8_t *out, const residue_narrow_reduced_t * __restrict x);

void decode(residue_narrow_reduced_t *out, const uint8_t *in);
#endif
