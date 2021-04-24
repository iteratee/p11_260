#define _DEFAULT_SOURCE
#include <string.h>
#include "f11_260.h"
#include "scalar.h"
#include "curve.h"
#include "constant_time.h"

__attribute__((__aligned__(64)))
const affine_pt_narrow_t B = {
  .x = {
    .limbs = {
      0x2862b8b, 0x0f08ed2, 0x06e65ee, 0x0c05991, 0x2b12b17,
      0x0049432, 0x33a3707, 0x16e5186, 0x2947e71, 0x0ed9bab,
      0,
    },
  },
  .y = {
    .limbs = {
      0x4,
    },
  },
};

void copy_projective_pt_narrow(
  projective_pt_narrow_t *result, const projective_pt_narrow_t *source) {

  for(int i = 0; i < NLIMBS; ++i) {
    result->x.limbs[i] = source->x.limbs[i];
    result->y.limbs[i] = source->y.limbs[i];
    result->z.limbs[i] = source->z.limbs[i];
  }
}

void copy_extended_pt_narrow(
  extended_pt_narrow_t *result,
  const extended_pt_narrow_t *source) {

  for(int i = 0; i < NLIMBS; ++i) {
    result->x.limbs[i] = source->x.limbs[i];
    result->y.limbs[i] = source->y.limbs[i];
    result->t.limbs[i] = source->t.limbs[i];
    result->z.limbs[i] = source->z.limbs[i];
  }
}

void copy_extended_pt_readd_narrow(
  extended_pt_readd_narrow_t *result,
  const extended_pt_readd_narrow_t *source) {
  for(int i = 0; i < NLIMBS; ++i) {
    result->x.limbs[i] = source->x.limbs[i];
    result->y.limbs[i] = source->y.limbs[i];
    result->dt.limbs[i] = source->dt.limbs[i];
    result->z.limbs[i] = source->z.limbs[i];
  }
}

void copy_extended_affine_pt_readd_narrow(
  extended_affine_pt_readd_narrow_t *result,
  const extended_affine_pt_readd_narrow_t *source) {
  for(int i = 0; i < NLIMBS; ++i) {
    result->x.limbs[i] = source->x.limbs[i];
    result->y.limbs[i] = source->y.limbs[i];
    result->dt.limbs[i] = source->dt.limbs[i];
  }
}

void negate_extended_pt_readd_narrow(
  extended_pt_readd_narrow_t *result,
  const extended_pt_readd_narrow_t *source) {
  for(int i = 0; i < NLIMBS; ++i) {
    result->x.limbs[i] = -source->x.limbs[i];
    result->y.limbs[i] = source->y.limbs[i];
    result->dt.limbs[i] = -source->dt.limbs[i];
    result->z.limbs[i] = source->z.limbs[i];
  }
}

void negate_extended_affine_pt_readd_narrow(
  extended_affine_pt_readd_narrow_t *result,
  const extended_affine_pt_readd_narrow_t *source) {
  for(int i = 0; i < NLIMBS; ++i) {
    result->x.limbs[i] = -source->x.limbs[i];
    result->dt.limbs[i] = -source->dt.limbs[i];
    result->y.limbs[i] = source->y.limbs[i];
  }
}

void affine_narrow_to_extended(
  extended_pt_narrow_t *result,
  const affine_pt_narrow_t * __restrict x) {

  for(int i = 0; i < NLIMBS; ++i) {
    result->x.limbs[i] = x->x.limbs[i];
    result->y.limbs[i] = x->y.limbs[i];
    result->z.limbs[i] = 0;
  }
  result->z.limbs[0] = 1;
  mul_narrow(&result->t, &result->x, &result->y);
}

void extended_to_projective_narrow(
  projective_pt_narrow_t *result, const extended_pt_narrow_t * __restrict x) {
  for(int i = 0; i < NLIMBS; ++i) {
    result->x.limbs[i] = x->x.limbs[i];
    result->y.limbs[i] = x->y.limbs[i];
    result->z.limbs[i] = x->z.limbs[i];
  }
}

void affine_to_readd_narrow(
  extended_pt_readd_narrow_t *result,
  const affine_pt_narrow_t * __restrict x) {

  for(int i = 0; i < NLIMBS; ++i) {
    result->x.limbs[i] = x->x.limbs[i];
    result->y.limbs[i] = x->y.limbs[i];
    result->z.limbs[i] = 0;
  }
  result->z.limbs[0] = 1;

  residue_narrow_t xy;
  mul_narrow(&xy, &x->x, &x->y);
  mul_narrow_const(&result->dt, &xy, D);
}

void extended_to_readd_narrow_neg(
  extended_pt_readd_narrow_t *result,
  const extended_pt_narrow_t * __restrict x) {

  for(int i = 0; i < NLIMBS; ++i) {
    result->x.limbs[i] = -(x->x.limbs[i]);
    result->y.limbs[i] = x->y.limbs[i];
    result->z.limbs[i] = x->z.limbs[i];
  }
  mul_narrow_const(&result->dt, &x->t, -D);
}

void affine_double(
  projective_pt_narrow_t *result,
  const affine_pt_narrow_t * __restrict x) {

  residue_narrow_t x_plus_y;
  residue_narrow_t a, b, e, e_tmp, g, g_minus_2, h;
  square_narrow(&a, &x->x);
  square_narrow(&b, &x->y);

  add_narrow(&x_plus_y, &x->x, &x->y);

  square_narrow(&e, &x_plus_y);
  sub_narrow(&e_tmp, &e, &a);
  sub_narrow(&e, &e_tmp, &b);
  add_narrow(&g, &a, &b);

  for (int i = 0; i < NLIMBS; ++i) {
    g_minus_2.limbs[i] = g.limbs[i];
  }
  g_minus_2.limbs[0] -= 2;

  sub_narrow(&h, &a, &b);
  mul_narrow(&result->x, &e, &g_minus_2);
  mul_narrow(&result->y, &g, &h);
  mul_narrow(&result->z, &g, &g_minus_2);
}

void affine_double_extended(
  extended_pt_narrow_t *result, const affine_pt_narrow_t * __restrict x) {

  residue_narrow_t x_plus_y;
  residue_narrow_t a, b, e, e_tmp, g, g_minus_2, h;
  square_narrow(&a, &x->x);
  square_narrow(&b, &x->y);

  add_narrow(&x_plus_y, &x->x, &x->y);
  square_narrow(&e, &x_plus_y);
  sub_narrow(&e_tmp, &e, &a);
  sub_narrow(&e, &e_tmp, &b);
  add_narrow(&g, &a, &b);

  for (int i = 0; i < NLIMBS; ++i) {
    g_minus_2.limbs[i] = g.limbs[i];
  }
  g_minus_2.limbs[0] -= 2;

  sub_narrow(&h, &a, &b);
  mul_narrow(&result->x, &e, &g_minus_2);
  mul_narrow(&result->y, &g, &h);
  mul_narrow(&result->t, &e, &h);
  mul_narrow(&result->z, &g, &g_minus_2);
}

void projective_double(
  projective_pt_narrow_t *result, const projective_pt_narrow_t *x) {

  residue_narrow_t x_plus_y;
  residue_narrow_t a, b, c, c_temp, e, e_tmp, f, g, h;
  add_narrow(&x_plus_y, &x->x, &x->y);
  square_narrow(&a, &x->x);
  square_narrow(&b, &x->y);
  square_narrow(&c_temp, &x->z);
  double_narrow(&c, &c_temp);

  square_narrow(&e, &x_plus_y);
  sub_narrow(&e_tmp, &e, &a);
  sub_narrow(&e, &e_tmp, &b);
  add_narrow(&g, &a, &b);
  sub_narrow(&f, &g, &c);
  sub_narrow(&h, &a, &b);

  mul_narrow(&result->x, &e, &f);
  mul_narrow(&result->y, &g, &h);
  mul_narrow(&result->z, &f, &g);
}

void projective_double_extended(
  extended_pt_narrow_t *result, const projective_pt_narrow_t * __restrict x) {

  residue_narrow_t x_plus_y;
  residue_narrow_t a, b, c, c_temp, e, e_tmp, f, g, h;
  add_narrow(&x_plus_y, &x->x, &x->y);
  square_narrow(&a, &x->x);
  square_narrow(&b, &x->y);
  square_narrow(&c_temp, &x->z);
  double_narrow(&c, &c_temp);

  square_narrow(&e, &x_plus_y);
  sub_narrow(&e_tmp, &e, &a);
  sub_narrow(&e, &e_tmp, &b);
  add_narrow(&g, &a, &b);
  sub_narrow(&f, &g, &c);
  sub_narrow(&h, &a, &b);

  mul_narrow(&result->x, &e, &f);
  mul_narrow(&result->y, &g, &h);
  mul_narrow(&result->t, &e, &h);
  mul_narrow(&result->z, &f, &g);
}

void extended_double_extended(
  extended_pt_narrow_t *result, const extended_pt_narrow_t *x) {

  residue_narrow_t x_plus_y;
  residue_narrow_t a, b, c, c_temp, e, e_tmp, f, g, h;
  add_narrow(&x_plus_y, &x->x, &x->y);
  square_narrow(&a, &x->x);
  square_narrow(&b, &x->y);
  square_narrow(&c_temp, &x->z);
  double_narrow(&c, &c_temp);

  square_narrow(&e, &x_plus_y);
  sub_narrow(&e_tmp, &e, &a);
  sub_narrow(&e, &e_tmp, &b);
  add_narrow(&g, &a, &b);
  sub_narrow(&f, &g, &c);
  sub_narrow(&h, &a, &b);

  mul_narrow(&result->x, &e, &f);
  mul_narrow(&result->z, &f, &g);
  mul_narrow(&result->y, &g, &h);
  mul_narrow(&result->t, &e, &h);
}

void projective_add(
  projective_pt_narrow_t *result, const projective_pt_narrow_t * __restrict x1,
  const projective_pt_narrow_t * __restrict x2) {

  residue_narrow_t x1_plus_y1, x2_plus_y2;
  residue_narrow_t a, b, c, d, e, e_temp, f, g, t1, t2;

  mul_narrow(&a, &x1->z, &x2->z);
  square_narrow(&b, &a);
  mul_narrow(&c, &x1->x, &x2->x);
  mul_narrow(&d, &x1->y, &x2->y);
  mul_narrow_const(&e_temp, &c, D);
  mul_narrow(&e, &e_temp, &d);

  sub_narrow(&f, &b, &e);
  add_narrow(&g, &b, &e);
  add_narrow(&x1_plus_y1, &x1->x, &x1->y);
  add_narrow(&x2_plus_y2, &x2->x, &x2->y);

  mul_narrow(&t1, &x1_plus_y1, &x2_plus_y2);
  sub_narrow(&t2, &t1, &c);
  sub_narrow(&t1, &t2, &d);
  mul_narrow(&t2, &t1, &f);
  mul_narrow(&result->x, &t2, &a);

  sub_narrow(&t1, &d, &c);
  mul_narrow(&t2, &t1, &g);
  mul_narrow(&result->y, &t2, &a);

  mul_narrow(&result->z, &f, &g);
}

void extended_add(
  projective_pt_narrow_t *result, const extended_pt_narrow_t * __restrict x1,
  const extended_pt_narrow_t * __restrict x2) {

  residue_narrow_t x1_plus_y1, x2_plus_y2;
  residue_narrow_t a, b, c, c_temp, d, e, e_temp, f, g, h;

  mul_narrow(&a, &x1->x, &x2->x);
  mul_narrow(&b, &x1->y, &x2->y);
  mul_narrow_const(&c_temp, &x1->t, D);
  mul_narrow(&c, &c_temp, &x2->t);
  mul_narrow(&d, &x1->z, &x2->z);

  add_narrow(&x1_plus_y1, &x1->x, &x1->y);
  add_narrow(&x2_plus_y2, &x2->x, &x2->y);
  mul_narrow(&e, &x1_plus_y1, &x2_plus_y2);
  sub_narrow(&e_temp, &e, &a);
  sub_narrow(&e, &e_temp, &b);
  sub_narrow(&f, &d, &c);
  add_narrow(&g, &d, &c);
  sub_narrow(&h, &b, &a);

  mul_narrow(&result->x, &e, &f);
  mul_narrow(&result->z, &f, &g);
  mul_narrow(&result->y, &g, &h);
}

void extended_add_extended(
  extended_pt_narrow_t *result, const extended_pt_narrow_t *x1,
  const extended_pt_narrow_t *x2) {

  residue_narrow_t x1_plus_y1, x2_plus_y2;
  residue_narrow_t a, b, c, c_temp, d, e, e_temp, f, g, h;

  mul_narrow(&a, &x1->x, &x2->x);
  mul_narrow(&b, &x1->y, &x2->y);
  mul_narrow_const(&c_temp, &x1->t, D);
  mul_narrow(&c, &c_temp, &x2->t);
  mul_narrow(&d, &x1->z, &x2->z);

  add_narrow(&x1_plus_y1, &x1->x, &x1->y);
  add_narrow(&x2_plus_y2, &x2->x, &x2->y);
  mul_narrow(&e, &x1_plus_y1, &x2_plus_y2);
  sub_narrow(&e_temp, &e, &a);
  sub_narrow(&e, &e_temp, &b);
  sub_narrow(&f, &d, &c);
  add_narrow(&g, &d, &c);
  sub_narrow(&h, &b, &a);

  mul_narrow(&result->x, &e, &f);
  mul_narrow(&result->z, &f, &g);
  mul_narrow(&result->y, &g, &h);
  mul_narrow(&result->t, &e, &h);
}

void extended_readd_narrow_extended(
  extended_pt_narrow_t *result, const extended_pt_narrow_t * __restrict x1,
  const extended_pt_readd_narrow_t * __restrict x2) {

  residue_narrow_t x1_plus_y1;
  residue_narrow_t x2_plus_y2;
  residue_narrow_t a, b, c, d, e, e_temp, f, g, h;

  mul_narrow(&a, &x1->x, &x2->x);
  mul_narrow(&b, &x1->y, &x2->y);
  mul_narrow(&c, &x1->t, &x2->dt);
  mul_narrow(&d, &x1->z, &x2->z);

  add_narrow(&x1_plus_y1, &x1->x, &x1->y);
  add_narrow(&x2_plus_y2, &x2->x, &x2->y);
  mul_narrow(&e, &x1_plus_y1, &x2_plus_y2);
  sub_narrow(&e_temp, &e, &a);
  sub_narrow(&e, &e_temp, &b);
  sub_narrow(&f, &d, &c);
  add_narrow(&g, &d, &c);
  sub_narrow(&h, &b, &a);

  mul_narrow(&result->x, &e, &f);
  mul_narrow(&result->z, &f, &g);
  mul_narrow(&result->y, &g, &h);
  mul_narrow(&result->t, &e, &h);
}

void extended_readd_narrow(
  projective_pt_narrow_t *result, const extended_pt_narrow_t * __restrict x1,
  const extended_pt_readd_narrow_t * __restrict x2) {

  residue_narrow_t x1_plus_y1;
  residue_narrow_t x2_plus_y2;
  residue_narrow_t a, b, c, d, e, e_temp, f, g, h;

  mul_narrow(&a, &x1->x, &x2->x);
  mul_narrow(&b, &x1->y, &x2->y);
  mul_narrow(&c, &x1->t, &x2->dt);
  mul_narrow(&d, &x1->z, &x2->z);

  add_narrow(&x1_plus_y1, &x1->x, &x1->y);
  add_narrow(&x2_plus_y2, &x2->x, &x2->y);
  mul_narrow(&e, &x1_plus_y1, &x2_plus_y2);
  sub_narrow(&e_temp, &e, &a);
  sub_narrow(&e, &e_temp, &b);
  sub_narrow(&f, &d, &c);
  add_narrow(&g, &d, &c);
  sub_narrow(&h, &b, &a);

  mul_narrow(&result->x, &e, &f);
  mul_narrow(&result->z, &f, &g);
  mul_narrow(&result->y, &g, &h);
}

#include <stdio.h>
// static void print_narrow(const residue_narrow_t *x, const char *prefix) {
//   printf("%s\n", prefix);
//   for (int i = 0; i < NLIMBS; ++i) {
//     printf("%#x\n", x->limbs[i]);
//   }
// }

void extended_readd_affine_narrow_extended(
  extended_pt_narrow_t *result, const extended_pt_narrow_t *x1,
  const extended_affine_pt_readd_narrow_t * __restrict x2) {

  residue_narrow_t x1_plus_y1;
  residue_narrow_t x2_plus_y2;
  residue_narrow_t a, b, c, e, e_temp, f, g, h;

  mul_narrow(&a, &x1->x, &x2->x);
  mul_narrow(&b, &x1->y, &x2->y);
  mul_narrow(&c, &x1->t, &x2->dt);

  add_narrow(&x1_plus_y1, &x1->x, &x1->y);
  add_narrow(&x2_plus_y2, &x2->x, &x2->y);
  mul_narrow(&e, &x1_plus_y1, &x2_plus_y2);
  sub_narrow(&e_temp, &e, &a);
  sub_narrow(&e, &e_temp, &b);
  sub_narrow(&f, &x1->z, &c);
  add_narrow(&g, &x1->z, &c);
  sub_narrow(&h, &b, &a);

  mul_narrow(&result->x, &e, &f);
  mul_narrow(&result->z, &f, &g);
  mul_narrow(&result->y, &g, &h);
  mul_narrow(&result->t, &e, &h);
}

void extended_readd_readd_narrow(
  extended_pt_readd_narrow_t *result,
  const extended_pt_narrow_t * __restrict x1,
  const extended_pt_readd_narrow_t * __restrict x2) {

  residue_narrow_t x1_plus_y1;
  residue_narrow_t x2_plus_y2;
  residue_narrow_t a, b, c, d, e, e_temp, f, g, h, t3;

  mul_narrow(&a, &x1->x, &x2->x);
  mul_narrow(&b, &x1->y, &x2->y);
  mul_narrow(&c, &x1->t, &x2->dt);
  mul_narrow(&d, &x1->z, &x2->z);

  add_narrow(&x1_plus_y1, &x1->x, &x1->y);
  add_narrow(&x2_plus_y2, &x2->x, &x2->y);
  mul_narrow(&e, &x1_plus_y1, &x2_plus_y2);
  sub_narrow(&e_temp, &e, &a);
  sub_narrow(&e, &e_temp, &b);
  sub_narrow(&f, &d, &c);
  add_narrow(&g, &d, &c);
  sub_narrow(&h, &b, &a);

  mul_narrow(&result->x, &e, &f);
  mul_narrow(&result->z, &f, &g);
  mul_narrow(&result->y, &g, &h);
  mul_narrow(&t3, &e, &h);

  mul_narrow_const(&result->dt, &t3, D);
}

void readd_to_projective(
  projective_pt_narrow_t *result,
  const extended_pt_readd_narrow_t * __restrict x) {

  copy_narrow(&result->x, &x->x);
  copy_narrow(&result->y, &x->y);
  copy_narrow(&result->z, &x->z);
}

void affine_readd_to_extended(
  extended_pt_narrow_t *result,
  const extended_affine_pt_readd_narrow_t * __restrict x) {

  copy_narrow(&result->x, &x->x);
  copy_narrow(&result->y, &x->y);
  mul_narrow(&result->t, &x->x, &x->y);
  for (int i = 0; i < NLIMBS; ++i) {
    result->z.limbs[i] = 0;
  }
  result->z.limbs[0] = 1;
}

void scalar_multiply(
  projective_pt_narrow_t *result, const affine_pt_narrow_t * __restrict x,
  const scalar_t * __restrict n) {

  scalar_t sabs_n;
  convert_to_sabs(&sabs_n, n);

  const int WINDOW_BITS = 5;
  const uint32_t WINDOW_MASK = (1 << WINDOW_BITS) - 1;
  const uint32_t LOOKUP_MASK = WINDOW_MASK >> 1;
  const int TABLE_SIZE = 16;
  extended_pt_readd_narrow_t table[TABLE_SIZE];

  extended_pt_narrow_t x2;
  affine_double_extended(&x2, x);
  affine_to_readd_narrow(&table[0], x);
  for (int i = 1; i < TABLE_SIZE; ++i) {
    extended_readd_readd_narrow(&table[i], &x2, &table[i-1]);
  }

  int i;
  int first = 1;
  // Set i to the highest i such that
  // a) i < SCALAR_BITS
  // b) i % WINDOW_BITS = 0

  projective_pt_narrow_t temp;
  extended_pt_narrow_t temp_ext;
  extended_pt_readd_narrow_t window_pt;

  i = SCALAR_BITS - ((SCALAR_BITS - 1) % WINDOW_BITS) - 1;
  for (; i >= 0; i -= WINDOW_BITS) {
    uint32_t bits = sabs_n.limbs[i/SCALAR_LIMB_BITS] >> (i % SCALAR_LIMB_BITS);
    if (i % SCALAR_LIMB_BITS > (SCALAR_LIMB_BITS - WINDOW_BITS) &&
        i / SCALAR_LIMB_BITS < SCALAR_LIMBS - 1) {

      bits |= sabs_n.limbs[i/SCALAR_LIMB_BITS + 1] <<
        (SCALAR_LIMB_BITS - i % SCALAR_LIMB_BITS);
    }

    bits &= WINDOW_MASK;
    int32_t invert = (bits >> (WINDOW_BITS - 1)) - 1;
    bits ^= invert;

    constant_time_extended_narrow_lookup(
      &window_pt, bits & LOOKUP_MASK, TABLE_SIZE, table);
    constant_time_cond_extended_negate(&window_pt, invert);

    if (first) {
      readd_to_projective(&temp, &window_pt);
      first = 0;
    } else {
      for (int i = 0; i < WINDOW_BITS - 1; ++i) {
        projective_double(&temp, &temp);
      }
      projective_double_extended(&temp_ext, &temp);
      extended_readd_narrow(&temp, &temp_ext, &window_pt);
    }
  }

  copy_projective_pt_narrow(result, &temp);
  explicit_bzero(&sabs_n, sizeof(sabs_n));
  explicit_bzero(&window_pt, sizeof(window_pt));
  explicit_bzero(table, sizeof(table));
  explicit_bzero(&temp, sizeof(temp));
  explicit_bzero(&temp_ext, sizeof(temp_ext));
}

void scalar_multiply_unsafe(
  projective_pt_narrow_t *result, const affine_pt_narrow_t * __restrict x,
  const scalar_t * __restrict n) {

  scalar_t sabs_n;
  convert_to_sabs(&sabs_n, n);

  const int WINDOW_BITS = 5;
  const uint32_t WINDOW_MASK = (1 << WINDOW_BITS) - 1;
  const uint32_t LOOKUP_MASK = WINDOW_MASK >> 1;
  const int TABLE_SIZE = 16;
  extended_pt_readd_narrow_t table[TABLE_SIZE];

  extended_pt_narrow_t x2;
  affine_double_extended(&x2, x);
  affine_to_readd_narrow(&table[0], x);
  for (int i = 1; i < TABLE_SIZE; ++i) {
    extended_readd_readd_narrow(&table[i], &x2, &table[i-1]);
  }

  int i;
  int first = 1;
  // Set i to the highest i such that
  // a) i < SCALAR_BITS
  // b) i % WINDOW_BITS = 0

  projective_pt_narrow_t temp;
  extended_pt_narrow_t temp_ext;
  extended_pt_readd_narrow_t window_pt;

  i = SCALAR_BITS - ((SCALAR_BITS - 1) % WINDOW_BITS) - 1;
  for (; i >= 0; i -= WINDOW_BITS) {
    uint32_t bits = sabs_n.limbs[i/SCALAR_LIMB_BITS] >> (i % SCALAR_LIMB_BITS);
    if (i % SCALAR_LIMB_BITS > (SCALAR_LIMB_BITS - WINDOW_BITS) &&
        i / SCALAR_LIMB_BITS < SCALAR_LIMBS - 1) {

      bits |= sabs_n.limbs[i/SCALAR_LIMB_BITS + 1] <<
        (SCALAR_LIMB_BITS - i % SCALAR_LIMB_BITS);
    }

    bits &= WINDOW_MASK;
    int32_t invert = (bits >> (WINDOW_BITS - 1)) - 1;
    bits ^= invert;

    copy_extended_pt_readd_narrow(&window_pt, &table[bits & LOOKUP_MASK]);
    if (invert) {
      negate_extended_pt_readd_narrow(&window_pt, &window_pt);
    }

    if (first) {
      readd_to_projective(&temp, &window_pt);
      first = 0;
    } else {
      for (int i = 0; i < WINDOW_BITS - 1; ++i) {
        projective_double(&temp, &temp);
      }
      projective_double_extended(&temp_ext, &temp);
      extended_readd_narrow(&temp, &temp_ext, &window_pt);
    }
  }

  copy_projective_pt_narrow(result, &temp);
}

int point_decompress(
  affine_pt_narrow_t *result,
  residue_narrow_reduced_t *y, int low_bit) {

  residue_narrow_t y_n;

  residue_narrow_t u;
  residue_narrow_t v;

  residue_narrow_t y2;
  residue_narrow_reduced_t temp;

  unnarrow_reduce(&y_n, y);
  square_narrow(&y2, &y_n);
  copy_narrow(&result->y, &y_n);

  sub_narrow(&u, &one_narrow, &y2);
  mul_narrow_const(&y2, &y2, D);
  sub_narrow(&v, &one_narrow, &y2);

  if (sqrt_inv_narrow(&result->x, &u, &v)) {
    narrow_partial_complete(&temp, &result->x);

    int x_is_odd = is_odd(&temp);
    if ((x_is_odd && !low_bit) || (low_bit && !x_is_odd)) {
      negate_narrow(&result->x, &result->x);
    }

    return 1;
  }

  return 0;
}
