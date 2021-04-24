// Types and functions for manipulating field elements

#ifndef F11_260_H
#define F11_260_H
#include <stdint.h>

#define NLIMBS_REDUCED 10
#define NLIMBS 11
#define T ((1 << 26) - 15)
#define TBITS 26
#define TMASK ((1 << 26) - 1)
#define T_CBITS 4
#define RESIDUE_LENGTH_BYTES 33

// Reduced to 10 limbs. For final results.
typedef struct residue_narrow_reduced {
  __attribute__((__aligned__(8)))
  int32_t limbs[NLIMBS_REDUCED];
} residue_narrow_reduced_t;

// 11 limbs.
typedef struct residue_narrow {
  __attribute__((__aligned__(64)))
  int32_t limbs[NLIMBS];
  int32_t pad[16 - NLIMBS];
} residue_narrow_t;

// 11 limbs.
// compatibility.
typedef struct residue_wide {
  __attribute__((__aligned__(64)))
  int64_t limbs[NLIMBS];
  int64_t pad[16 - NLIMBS];
} residue_wide_t;

residue_wide_t zero_wide;
residue_wide_t one_wide;
residue_narrow_t zero_narrow;
residue_narrow_t one_narrow;

// Shrink to 32 bits. Assumes reduction has already occurred, and wide storage
// is being used for vector compatibility.
void narrow(residue_narrow_t *result, const residue_wide_t * __restrict w);

// Reduce to 10 limbs. Useful for debugging
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

int is_odd(residue_narrow_reduced_t *x);

// Produce a 32-bit entry with 11 limbs
static inline void unnarrow_reduce(
  residue_narrow_t *result, const residue_narrow_reduced_t * __restrict x) {

  result->limbs[10] = 0;
  for (int i = 0; i < NLIMBS_REDUCED; ++i) {
    result->limbs[i] = x->limbs[i];
  }
}

// Produce a 64-bit residue
void widen(
  residue_wide_t *result, const residue_narrow_t * __restrict x);

// Copy a 64-bit residue
void copy_wide(
  residue_wide_t *result, const residue_wide_t * __restrict x);

// Copy a 32-bit residue
void copy_narrow(
  residue_narrow_t *result, const residue_narrow_t * __restrict x);

void copy_narrow_reduced(
  residue_narrow_reduced_t *result,
  const residue_narrow_reduced_t * __restrict x);

// Subtract 2 11x32-bit residues.
void sub_narrow(
  residue_narrow_t *result, const residue_narrow_t * __restrict x,
  const residue_narrow_t * __restrict y);

void negate_wide(residue_wide_t *result, const residue_wide_t *x);

void negate_narrow(residue_narrow_t *result, const residue_narrow_t *x);

// Add 2 11x32-bit residues.
void add_narrow(
  residue_narrow_t *result, const residue_narrow_t * __restrict x,
  const residue_narrow_t * __restrict y);

// Add 2 11x64-bit residues.
void add_wide(
  residue_wide_t *result, const residue_wide_t * __restrict x,
  const residue_wide_t * __restrict y);

// Scale a wide residue by 2.
void double_narrow(
  residue_narrow_t *result, const residue_narrow_t * __restrict x);

// Multiply two narrow residues and produce a wide result. The result is reduced
// to 32 bits.
void mul_narrow(
  residue_narrow_t *result, const residue_narrow_t *x,
  const residue_narrow_t *y);

// Multiply a residue by a constant.
void mul_narrow_const(
  residue_narrow_t *result, const residue_narrow_t * __restrict x, int32_t d);

// Square a narrow residue and produce a narrow result. The result is reduced to
// 32 bits.
void square_narrow(
  residue_narrow_t *result, const residue_narrow_t *x);

// Approximately divide each coefficient by t. Carry the results.
void reduce_step_narrow(
  residue_narrow_t *result, const residue_narrow_t *x);

// Approximately divide each coefficient by t. Carry the results.
void reduce_step_wide(
  residue_wide_t *result, const residue_wide_t *x);

// Invert via fermat's theorem
void invert_narrow(
  residue_narrow_t *result, const residue_narrow_t * __restrict x);

// Compute combined inverse and square root
// returns true if x/y was a quadratic residue, and false otherwise.
int sqrt_inv_narrow(
  residue_narrow_t *result, const residue_narrow_t * __restrict x,
  const residue_narrow_t * __restrict y);

// Returns true if x == y. Computes in constant time.
int equal_narrow(const residue_narrow_t * x, const residue_narrow_t * y);

int equal_narrow_reduced(
  const residue_narrow_reduced_t * x, const residue_narrow_reduced_t * y);

void encode(uint8_t *out, const residue_narrow_reduced_t * __restrict x);
void encode_compressed(
  uint8_t *out, const residue_narrow_reduced_t * __restrict x, int is_odd);

void decode(residue_narrow_reduced_t *out, const uint8_t *in);
#endif
