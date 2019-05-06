#define _DEFAULT_SOURCE
#include <string.h>
#include "comb.h"
#include "curve.h"
#include "constant_time.h"
#include "f11_260.h"

#include <stdio.h>
static void print_wide(const residue_wide_t *x) {
  printf("[");
  for (int i = 0; i < NLIMBS; ++i) {
    printf(" %#lx,", x->limbs[i]);
    // printf("x[%d]: %d\n", i, x[i]);
  }
  printf(" ]");
  printf("\n");
}

static void print_as_narrow_reduced(
  const char *label, const residue_wide_t *x) {
  residue_narrow_t x_n;
  residue_narrow_reduced_t x_nr;
  narrow(&x_n, x);
  narrow_reduce(&x_nr, &x_n);
  printf("%s: [", label);
  for (int i = 0; i < NLIMBS_REDUCED; ++i) {
    printf(" %#x,", x_nr.limbs[i]);
    // printf("x[%d]: %d\n", i, x[i]);
  }
  printf(" ]");
  printf("\n");
}

static void print_projective_pt_wide(const projective_pt_wide_t *x) {
  print_as_narrow_reduced("x", &x->x);
  print_as_narrow_reduced("y", &x->y);
  print_as_narrow_reduced("z", &x->z);
}

static void print_extended_pt_wide(const extended_pt_wide_t *x) {
  print_as_narrow_reduced("x", &x->x);
  print_as_narrow_reduced("y", &x->y);
  print_as_narrow_reduced("t", &x->t);
  print_as_narrow_reduced("z", &x->z);
}

static void print_extended_pt_readd_wide(const extended_pt_readd_wide_t *x) {
  print_as_narrow_reduced("x", &x->x);
  print_as_narrow_reduced("y", &x->y);
  print_as_narrow_reduced("dt", &x->dt);
  print_as_narrow_reduced("z", &x->z);
}

void compute_comb_set(
  sabs_comb_set_t *result, const affine_pt_narrow_reduced_t *base_pt) {
  sabs_comb_set_wide_t result_t;
  teeth_set_t teeth_sets[4];
  extended_pt_wide_t everything_pts[4];

  projective_pt_wide_t temp;
  extended_pt_wide_t temp_ext;
  affine_double_extended(&temp_ext, base_pt);
  printf("single double:\n");
  print_extended_pt_wide(&temp_ext);
  extended_to_readd_wide_neg(&teeth_sets[0].teeth[0], &temp_ext);
  printf("first tooth:\n");
  print_extended_pt_readd_wide(&teeth_sets[0].teeth[0]);
  affine_narrow_reduced_to_extended(&everything_pts[0], base_pt);
  printf("initial everything pt #0:\n");
  print_extended_pt_wide(&everything_pts[0]);
  for (int i = 1; i < COMB_TEETH * COMB_COUNT; ++i) {
    extended_to_projective_wide(&temp, &temp_ext);
    for (int j = 0; j < COMB_SEPARATION - 2; ++j) {
      projective_double(&temp, &temp);
    }
    projective_double_extended(&temp_ext, &temp);
    printf("after 12th double:\n");
    print_extended_pt_wide(&temp_ext);
    if (i % COMB_TEETH == 0) {
      copy_extended_pt_wide(&everything_pts[i/COMB_TEETH], &temp_ext);
      printf("initial everything pt #%d:\n", i/COMB_TEETH);
      print_extended_pt_wide(&everything_pts[i/COMB_TEETH]);
    } else {
      extended_add_extended(
        &everything_pts[i/COMB_TEETH], &everything_pts[i/COMB_TEETH],
        &temp_ext);
      printf("updated everything pt #%d:\n", i/COMB_TEETH);
      print_extended_pt_wide(&everything_pts[i/COMB_TEETH]);
    }
    extended_double_extended(&temp_ext, &temp_ext);
    extended_to_readd_wide_neg(
      &teeth_sets[i/COMB_TEETH].teeth[i % COMB_TEETH], &temp_ext);
    printf("next tooth:\n");
    print_extended_pt_readd_wide(&teeth_sets[i / COMB_TEETH].teeth[ i % COMB_TEETH]);
  }
  // We now have all the precomputation necessary to walk the gray codes and
  // compute the table entries.
  int entry = COMB_TABLE_SIZE - 1;
  for (int i = 0; i < COMB_COUNT; ++i) {
    extended_to_projective_wide(
      &result_t.combs[i].table[entry], &everything_pts[i]);
  }
  printf("initial entry is: %#x\n", entry);
  for (int i = 1; i < COMB_TABLE_SIZE; ++i) {
    int j;
    for (j = 0; j < COMB_TEETH - 1; ++j) {
      int bit = 1 << (j + 1);
      int mask = bit - 1;
      int half_bit = bit >> 1;
      if ((i & mask) == half_bit) {
        entry ^= half_bit;
        printf("entry is: %#x\n", entry);
        break;
      }
    }
    for (int k = 0; k < COMB_COUNT; ++k) {
      extended_readd_wide_extended(
        &everything_pts[k], &everything_pts[k], &teeth_sets[k].teeth[j]);
      negate_extended_pt_readd_wide(
        &teeth_sets[k].teeth[j], &teeth_sets[k].teeth[j]);
      extended_to_projective_wide(
        &result_t.combs[k].table[entry], &everything_pts[k]);
    }
  }
  reduce_comb_set(result, &result_t);
}

// Note not const. Stomps on source
// Leaves source with invalid z-values, but if they are all set to one, it will
// be correct again. Shrinks combs down to affine with a single inversion. Also
// narrows and reduces to 10 limbs.
void reduce_comb_set(sabs_comb_set_t *result, sabs_comb_set_wide_t *source) {
  residue_wide_t z_left;
  residue_wide_t z_right;

  residue_wide_t z_inv;

  copy_wide(&z_left, &source->combs[0].table[0].z);
  for (int i = 1; i < COMB_TABLE_SIZE * COMB_COUNT; ++i) {
    mul_wide(&source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].x,
             &source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].x,
             &z_left);
    mul_wide(&source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].y,
             &source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].y,
             &z_left);
    mul_wide(&z_left, &z_left,
             &source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].z);
  }

  copy_wide(&z_right,
      &source->combs[COMB_COUNT - 1].table[COMB_TABLE_SIZE - 1].z);

  for (int i = COMB_TABLE_SIZE * COMB_COUNT - 2; i >= 0; --i) {
    mul_wide(&source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].x,
             &source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].x,
             &z_right);
    mul_wide(&source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].y,
             &source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].y,
             &z_right);
    mul_wide(&z_right, &z_right,
             &source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].z);
  }

  invert_wide(&z_inv, &z_right);
  for (int i = 0; i < COMB_TABLE_SIZE * COMB_COUNT; ++i) {
    residue_wide_t xy;
    residue_narrow_t temp_narrow;

    mul_wide(&source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].x,
             &source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].x,
             &z_inv);
    narrow(&temp_narrow,
           &source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].x);
    narrow_reduce(
      &result->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].x,
      &temp_narrow);

    mul_wide(&source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].y,
             &source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].y,
             &z_inv);
    narrow(&temp_narrow,
           &source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].y);
    narrow_reduce(
      &result->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].y,
      &temp_narrow);

    mul_wide(&xy,
             &source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].x,
             &source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].y);
    mul_wide_const(&xy, &xy, D);
    narrow(&temp_narrow, &xy);
    narrow_reduce(
      &result->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].dt,
      &temp_narrow);
  }
}

void scalar_comb_multiply(
  projective_pt_wide_t *result, const sabs_comb_set_t * __restrict comb,
  const scalar_t * __restrict n) {

  scalar_t sabs_n;
  convert_to_sabs(&sabs_n, n);

  extended_pt_wide_t temp;
  extended_affine_pt_readd_narrow_t table_pt;

  for (int i = COMB_SEPARATION - 1; i >= 0; --i) {
    if (i != COMB_SEPARATION - 1) {
      extended_double_extended(&temp, &temp);
    }
    for (int j = 0; j < COMB_COUNT; ++j) {
      int entry = 0;

      for (int k = 0; k < COMB_TEETH; ++k) {
        int bit = i + COMB_SEPARATION * (k + j * COMB_TEETH);
        if (bit < SCALAR_LIMB_BITS) {
          entry |= (sabs_n.limbs[bit / SCALAR_LIMB_BITS] >>
              (bit % SCALAR_LIMB_BITS) & 1) << k;
        }
      }

      int32_t invert = (entry >> (COMB_TEETH - 1)) - 1;
      entry ^= invert;

      constant_time_extended_affine_narrow_reduced_lookup(
        &table_pt, entry & COMB_LOOKUP_MASK, COMB_TABLE_SIZE,
        comb->combs[j].table);

      constant_time_cond_extended_affine_negate(&table_pt, invert);

      if (i == (COMB_SEPARATION - 1) && j == 0) {
        affine_readd_to_extended(&temp, &table_pt);
      } else {
        extended_readd_affine_narrow_extended(
          &temp, &temp, &table_pt);
      }
    }
  }

  extended_to_projective_wide(result, &temp);
  explicit_bzero(&sabs_n, sizeof(sabs_n));
  explicit_bzero(&table_pt, sizeof(table_pt));
  explicit_bzero(&temp, sizeof(temp));
}
