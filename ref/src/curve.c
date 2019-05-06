#define _DEFAULT_SOURCE
#include <string.h>
#include "f11_260.h"
#include "scalar.h"
#include "curve.h"
#include "constant_time.h"

const affine_pt_narrow_reduced_t B = {
  .x = {
    .limbs = {
      0x2862b8b, 0x0f08ed2, 0x06e65ee, 0x0c05991, 0x2b12b17,
      0x0049432, 0x33a3707, 0x16e5186, 0x2947e71, 0x0ed9bab,
    },
  },
  .y = {
    .limbs = {
      0x4, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
    },
  },
};

void copy_projective_pt_wide(
  projective_pt_wide_t *result, const projective_pt_wide_t *source) {

  for(int i = 0; i < NLIMBS; ++i) {
    result->x.limbs[i] = source->x.limbs[i];
    result->y.limbs[i] = source->y.limbs[i];
    result->z.limbs[i] = source->z.limbs[i];
  }
}

void copy_extended_pt_wide(
  extended_pt_wide_t *result,
  const extended_pt_wide_t *source) {

  for(int i = 0; i < NLIMBS; ++i) {
    result->x.limbs[i] = source->x.limbs[i];
    result->y.limbs[i] = source->y.limbs[i];
    result->t.limbs[i] = source->t.limbs[i];
    result->z.limbs[i] = source->z.limbs[i];
  }
}

void copy_extended_pt_readd_wide(
  extended_pt_readd_wide_t *result,
  const extended_pt_readd_wide_t *source) {

  for(int i = 0; i < NLIMBS; ++i) {
    result->x.limbs[i] = source->x.limbs[i];
    result->y.limbs[i] = source->y.limbs[i];
    result->dt.limbs[i] = source->dt.limbs[i];
    result->z.limbs[i] = source->z.limbs[i];
  }
}

void negate_extended_pt_readd_wide(
  extended_pt_readd_wide_t *result,
  const extended_pt_readd_wide_t *source) {
  for(int i = 0; i < NLIMBS; ++i) {
    result->x.limbs[i] = -source->x.limbs[i];
    result->y.limbs[i] = source->y.limbs[i];
    result->dt.limbs[i] = -source->dt.limbs[i];
    result->z.limbs[i] = source->z.limbs[i];
  }
}

void affine_narrow_reduced_to_extended(
  extended_pt_wide_t *result,
  const affine_pt_narrow_reduced_t * __restrict x) {

  for(int i = 0; i < NLIMBS_REDUCED; ++i) {
    result->x.limbs[i+1] = x->x.limbs[i];
    result->y.limbs[i+1] = x->y.limbs[i];
    result->z.limbs[i+1] = 0;
  }
  result->x.limbs[0] = result->x.limbs[NLIMBS - 1] = 0;
  result->y.limbs[0] = result->y.limbs[NLIMBS - 1] = 0;
  result->z.limbs[0] = result->z.limbs[NLIMBS - 1] = 0;
  result->z.limbs[1] = 1;
  mul_wide(&result->t, &result->x, &result->y);
}

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

static void print_narrow(const residue_narrow_t *x) {
  printf("[");
  for (int i = 0; i < NLIMBS; ++i) {
    printf(" %#x,", x->limbs[i]);
    // printf("x[%d]: %d\n", i, x[i]);
  }
  printf(" ]");
  printf("\n");
}

void extended_to_projective_wide(
  projective_pt_wide_t *result, const extended_pt_wide_t * __restrict x) {
  for(int i = 0; i < NLIMBS; ++i) {
    result->x.limbs[i] = x->x.limbs[i];
    result->y.limbs[i] = x->y.limbs[i];
    result->z.limbs[i] = x->z.limbs[i];
  }
}

void affine_to_readd_narrow_reduced(
  extended_pt_readd_narrow_reduced_t *result,
  const affine_pt_narrow_reduced_t * __restrict x) {

  for(int i = 0; i < NLIMBS_REDUCED; ++i) {
    result->x.limbs[i] = x->x.limbs[i];
    result->y.limbs[i] = x->y.limbs[i];
    result->z.limbs[i] = 0;
  }
  result->z.limbs[0] = 1;

  residue_wide_t xy;
  residue_wide_t dt_wide;
  residue_narrow_t dt_narrow;
  residue_narrow_t x_n, y_n;
  unnarrow_reduce(&x_n, &x->x);
  unnarrow_reduce(&y_n, &x->y);
  mul_narrow(&xy, &x_n, &y_n);
  mul_wide_const(&dt_wide, &xy, D);
  narrow(&dt_narrow, &dt_wide);
  narrow_reduce(&result->dt, &dt_narrow);
}

void affine_to_readd_wide(
  extended_pt_readd_wide_t *result,
  const affine_pt_narrow_reduced_t * __restrict x) {

  result->x.limbs[0] = 0;
  result->y.limbs[0] = 0;
  result->z.limbs[0] = 0;

  result->x.limbs[NLIMBS - 1] = 0;
  result->y.limbs[NLIMBS - 1] = 0;
  result->z.limbs[NLIMBS - 1] = 0;
  for(int i = 0; i < NLIMBS_REDUCED; ++i) {
    result->x.limbs[i+1] = x->x.limbs[i];
    result->y.limbs[i+1] = x->y.limbs[i];
    result->z.limbs[i+1] = 0;
  }
  result->z.limbs[1] = 1;

  residue_wide_t xy;
  residue_narrow_t x_n, y_n;
  unnarrow_reduce(&x_n, &x->x);
  unnarrow_reduce(&y_n, &x->y);
  mul_narrow(&xy, &x_n, &y_n);
  mul_wide_const(&result->dt, &xy, D);
}

void extended_to_readd_wide_neg(
  extended_pt_readd_wide_t *result,
  const extended_pt_wide_t * __restrict x) {

  for(int i = 0; i < NLIMBS; ++i) {
    result->x.limbs[i] = -(x->x.limbs[i]);
    result->y.limbs[i] = x->y.limbs[i];
    result->z.limbs[i] = x->z.limbs[i];
  }
  mul_wide_const(&result->dt, &x->t, -D);
}

void affine_double(
  projective_pt_wide_t *result,
  const affine_pt_narrow_reduced_t * __restrict x) {

  residue_narrow_t x_pad, y_pad;
  unnarrow_reduce(&x_pad, &x->x);
  unnarrow_reduce(&y_pad, &x->y);

  residue_narrow_t x_plus_y;
  residue_wide_t a, b, e, e_tmp, g, g_minus_2, h;
  square_narrow(&a, &x_pad);
  square_narrow(&b, &y_pad);

  add_narrow(&x_plus_y, &x_pad, &y_pad);
  square_narrow(&e, &x_plus_y);
  sub_wide(&e_tmp, &e, &a);
  sub_wide(&e, &e_tmp, &b);
  add_wide(&g, &a, &b);

  for (int i = 0; i < NLIMBS; ++i) {
    g_minus_2.limbs[i] = g.limbs[i];
  }
  g_minus_2.limbs[1] -= 2;

  sub_wide(&h, &a, &b);
  mul_wide(&result->x, &e, &g_minus_2);
  mul_wide(&result->y, &g, &h);
  mul_wide(&result->z, &g, &g_minus_2);
}

void affine_double_extended(
  extended_pt_wide_t *result, const affine_pt_narrow_reduced_t * __restrict x) {

  residue_narrow_t x_pad, y_pad;
  unnarrow_reduce(&x_pad, &x->x);
  unnarrow_reduce(&y_pad, &x->y);

  residue_narrow_t x_plus_y;
  residue_wide_t a, b, e, e_tmp, g, g_minus_2, h;
  square_narrow(&a, &x_pad);
  square_narrow(&b, &y_pad);

  add_narrow(&x_plus_y, &x_pad, &y_pad);
  square_narrow(&e, &x_plus_y);
  sub_wide(&e_tmp, &e, &a);
  sub_wide(&e, &e_tmp, &b);
  add_wide(&g, &a, &b);

  for (int i = 0; i < NLIMBS; ++i) {
    g_minus_2.limbs[i] = g.limbs[i];
  }
  g_minus_2.limbs[1] -= 2;

  sub_wide(&h, &a, &b);
  mul_wide(&result->x, &e, &g_minus_2);
  mul_wide(&result->y, &g, &h);
  mul_wide(&result->t, &e, &h);
  mul_wide(&result->z, &g, &g_minus_2);
}

void projective_double(
  projective_pt_wide_t *result, const projective_pt_wide_t *x) {

  residue_wide_t x_plus_y;
  residue_wide_t a, b, c, c_temp, e, e_tmp, f, g, h;
  square_wide(&a, &x->x);
  square_wide(&b, &x->y);
  square_wide(&c_temp, &x->z);
  double_wide(&c, &c_temp);

  add_wide(&x_plus_y, &x->x, &x->y);
  square_wide(&e, &x_plus_y);
  sub_wide(&e_tmp, &e, &a);
  sub_wide(&e, &e_tmp, &b);
  add_wide(&g, &a, &b);
  sub_wide(&f, &g, &c);
  sub_wide(&h, &a, &b);

  mul_wide(&result->x, &e, &f);
  mul_wide(&result->y, &g, &h);
  mul_wide(&result->z, &f, &g);
}

void projective_double_extended(
  extended_pt_wide_t *result, const projective_pt_wide_t * __restrict x) {

  residue_wide_t x_plus_y;
  residue_wide_t a, b, c, c_temp, e, e_tmp, f, g, h;
  square_wide(&a, &x->x);
  square_wide(&b, &x->y);
  square_wide(&c_temp, &x->z);
  double_wide(&c, &c_temp);

  add_wide(&x_plus_y, &x->x, &x->y);
  square_wide(&e, &x_plus_y);
  sub_wide(&e_tmp, &e, &a);
  sub_wide(&e, &e_tmp, &b);
  add_wide(&g, &a, &b);
  sub_wide(&f, &g, &c);
  sub_wide(&h, &a, &b);

  mul_wide(&result->x, &e, &f);
  mul_wide(&result->y, &g, &h);
  mul_wide(&result->t, &e, &h);
  mul_wide(&result->z, &f, &g);
}

void extended_double_extended(
  extended_pt_wide_t *result, const extended_pt_wide_t *x) {

  residue_wide_t x_plus_y;
  residue_wide_t a, b, c, c_temp, e, e_tmp, f, g, h;
  square_wide(&a, &x->x);
  square_wide(&b, &x->y);
  square_wide(&c_temp, &x->z);
  double_wide(&c, &c_temp);

  add_wide(&x_plus_y, &x->x, &x->y);
  square_wide(&e, &x_plus_y);
  sub_wide(&e_tmp, &e, &a);
  sub_wide(&e, &e_tmp, &b);
  add_wide(&g, &a, &b);
  sub_wide(&f, &g, &c);
  sub_wide(&h, &a, &b);

  mul_wide(&result->x, &e, &f);
  mul_wide(&result->y, &g, &h);
  mul_wide(&result->t, &e, &h);
  mul_wide(&result->z, &f, &g);
}

void extended_add(
  projective_pt_wide_t *result, const extended_pt_wide_t * __restrict x1,
  const extended_pt_wide_t * __restrict x2) {

  residue_wide_t x1_plus_y1, x2_plus_y2;
  residue_wide_t a, b, c, c_temp, d, e, e_temp, f, g, h;

  mul_wide(&a, &x1->x, &x2->x);
  mul_wide(&b, &x1->y, &x2->y);
  mul_wide_const(&c_temp, &x1->t, D);
  mul_wide(&c, &c_temp, &x2->t);
  mul_wide(&d, &x1->z, &x2->z);

  add_wide(&x1_plus_y1, &x1->x, &x1->y);
  add_wide(&x2_plus_y2, &x2->x, &x2->y);
  mul_wide(&e, &x1_plus_y1, &x2_plus_y2);
  sub_wide(&e_temp, &e, &a);
  sub_wide(&e, &e_temp, &b);
  sub_wide(&f, &d, &c);
  add_wide(&g, &d, &c);
  sub_wide(&h, &b, &a);

  mul_wide(&result->x, &e, &f);
  mul_wide(&result->y, &g, &h);
  mul_wide(&result->z, &f, &g);
}

void extended_add_extended(
  extended_pt_wide_t *result, const extended_pt_wide_t *x1,
  const extended_pt_wide_t *x2) {

  residue_wide_t x1_plus_y1, x2_plus_y2;
  residue_wide_t a, b, c, c_temp, d, e, e_temp, f, g, h;

  mul_wide(&a, &x1->x, &x2->x);
  mul_wide(&b, &x1->y, &x2->y);
  mul_wide_const(&c_temp, &x1->t, D);
  mul_wide(&c, &c_temp, &x2->t);
  mul_wide(&d, &x1->z, &x2->z);

  add_wide(&x1_plus_y1, &x1->x, &x1->y);
  add_wide(&x2_plus_y2, &x2->x, &x2->y);
  mul_wide(&e, &x1_plus_y1, &x2_plus_y2);
  sub_wide(&e_temp, &e, &a);
  sub_wide(&e, &e_temp, &b);
  sub_wide(&f, &d, &c);
  add_wide(&g, &d, &c);
  sub_wide(&h, &b, &a);

  mul_wide(&result->x, &e, &f);
  mul_wide(&result->y, &g, &h);
  mul_wide(&result->t, &e, &h);
  mul_wide(&result->z, &f, &g);
}

void extended_readd_wide_extended(
  extended_pt_wide_t *result,
  const extended_pt_wide_t *x1,
  const extended_pt_readd_wide_t * __restrict x2) {

  residue_wide_t x1_plus_y1, x2_plus_y2;
  residue_wide_t a, b, c, d, e, e_temp, f, g, h;

  mul_wide(&a, &x1->x, &x2->x);
  mul_wide(&b, &x1->y, &x2->y);
  mul_wide(&c, &x1->t, &x2->dt);
  mul_wide(&d, &x1->z, &x2->z);

  add_wide(&x1_plus_y1, &x1->x, &x1->y);
  add_wide(&x2_plus_y2, &x2->x, &x2->y);
  mul_wide(&e, &x1_plus_y1, &x2_plus_y2);
  sub_wide(&e_temp, &e, &a);
  sub_wide(&e, &e_temp, &b);
  sub_wide(&f, &d, &c);
  add_wide(&g, &d, &c);
  sub_wide(&h, &b, &a);

  mul_wide(&result->x, &e, &f);
  mul_wide(&result->y, &g, &h);
  mul_wide(&result->z, &f, &g);
}

void extended_readd_narrow_extended(
  extended_pt_wide_t *result, const extended_pt_wide_t * __restrict x1,
  const extended_pt_readd_narrow_t * __restrict x2) {

  residue_wide_t x1_plus_y1;
  residue_narrow_t x2_plus_y2;
  residue_wide_t a, b, c, d, e, e_temp, f, g, h;

  mul_wide_narrow(&a, &x1->x, &x2->x);
  mul_wide_narrow(&b, &x1->y, &x2->y);
  mul_wide_narrow(&c, &x1->t, &x2->dt);
  mul_wide_narrow(&d, &x1->z, &x2->z);

  add_wide(&x1_plus_y1, &x1->x, &x1->y);
  add_narrow(&x2_plus_y2, &x2->x, &x2->y);
  mul_wide_narrow(&e, &x1_plus_y1, &x2_plus_y2);
  sub_wide(&e_temp, &e, &a);
  sub_wide(&e, &e_temp, &b);
  sub_wide(&f, &d, &c);
  add_wide(&g, &d, &c);
  sub_wide(&h, &b, &a);

  mul_wide(&result->x, &e, &f);
  mul_wide(&result->y, &g, &h);
  mul_wide(&result->t, &e, &h);
  mul_wide(&result->z, &f, &g);
}

void extended_readd_narrow(
  projective_pt_wide_t *result, const extended_pt_wide_t * __restrict x1,
  const extended_pt_readd_narrow_t * __restrict x2) {

  residue_wide_t x1_plus_y1;
  residue_narrow_t x2_plus_y2;
  residue_wide_t a, b, c, d, e, e_temp, f, g, h;

  mul_wide_narrow(&a, &x1->x, &x2->x);
  mul_wide_narrow(&b, &x1->y, &x2->y);
  mul_wide_narrow(&c, &x1->t, &x2->dt);
  mul_wide_narrow(&d, &x1->z, &x2->z);

  add_wide(&x1_plus_y1, &x1->x, &x1->y);
  add_narrow(&x2_plus_y2, &x2->x, &x2->y);
  mul_wide_narrow(&e, &x1_plus_y1, &x2_plus_y2);
  sub_wide(&e_temp, &e, &a);
  sub_wide(&e, &e_temp, &b);
  sub_wide(&f, &d, &c);
  add_wide(&g, &d, &c);
  sub_wide(&h, &b, &a);

  mul_wide(&result->x, &e, &f);
  mul_wide(&result->y, &g, &h);
  mul_wide(&result->z, &f, &g);
}

void extended_readd_affine_narrow_extended(
  extended_pt_wide_t *result, const extended_pt_wide_t *x1,
  const extended_affine_pt_readd_narrow_t * __restrict x2) {

  residue_wide_t x1_plus_y1;
  residue_narrow_t x2_plus_y2;
  residue_wide_t a, b, c, e, e_temp, f, g, h;

  mul_wide_narrow(&a, &x1->x, &x2->x);
  mul_wide_narrow(&b, &x1->y, &x2->y);
  mul_wide_narrow(&c, &x1->t, &x2->dt);

  add_wide(&x1_plus_y1, &x1->x, &x1->y);
  add_narrow(&x2_plus_y2, &x2->x, &x2->y);
  mul_wide_narrow(&e, &x1_plus_y1, &x2_plus_y2);
  sub_wide(&e_temp, &e, &a);
  sub_wide(&e, &e_temp, &b);
  sub_wide(&f, &x1->z, &c);
  add_wide(&g, &x1->z, &c);
  sub_wide(&h, &b, &a);

  mul_wide(&result->x, &e, &f);
  mul_wide(&result->y, &g, &h);
  mul_wide(&result->t, &e, &h);
  mul_wide(&result->z, &f, &g);
}

void extended_readd_readd_narrow_reduce(
  extended_pt_readd_narrow_reduced_t *result,
  const extended_pt_wide_t * __restrict x1,
  const extended_pt_readd_narrow_reduced_t * __restrict x2) {

  residue_narrow_t x2_narrow;
  residue_narrow_t y2_narrow;
  residue_narrow_t dt2_narrow;
  residue_narrow_t z2_narrow;
  residue_wide_t x1_plus_y1;
  residue_narrow_t x2_plus_y2;
  residue_wide_t a, b, c, d, e, e_temp, f, g, h, x3, y3, t3, dt3, z3;
  residue_narrow_t x3_narrow, y3_narrow, dt3_narrow, z3_narrow;

  unnarrow_reduce(&x2_narrow, &x2->x);
  unnarrow_reduce(&y2_narrow, &x2->y);
  unnarrow_reduce(&dt2_narrow, &x2->dt);
  unnarrow_reduce(&z2_narrow, &x2->z);

  mul_wide_narrow(&a, &x1->x, &x2_narrow);
  mul_wide_narrow(&b, &x1->y, &y2_narrow);
  mul_wide_narrow(&c, &x1->t, &dt2_narrow);
  mul_wide_narrow(&d, &x1->z, &z2_narrow);

  add_wide(&x1_plus_y1, &x1->x, &x1->y);
  add_narrow(&x2_plus_y2, &x2_narrow, &y2_narrow);
  mul_wide_narrow(&e, &x1_plus_y1, &x2_plus_y2);
  sub_wide(&e_temp, &e, &a);
  sub_wide(&e, &e_temp, &b);
  sub_wide(&f, &d, &c);
  add_wide(&g, &d, &c);
  sub_wide(&h, &b, &a);

  mul_wide(&x3, &e, &f);
  mul_wide(&y3, &g, &h);
  mul_wide(&t3, &e, &h);
  mul_wide_const(&dt3, &t3, D);
  mul_wide(&z3, &f, &g);

  narrow(&x3_narrow, &x3);
  narrow(&y3_narrow, &y3);
  narrow(&dt3_narrow, &dt3);
  narrow(&z3_narrow, &z3);

  narrow_reduce(&result->x, &x3_narrow);
  narrow_reduce(&result->y, &y3_narrow);
  narrow_reduce(&result->dt, &dt3_narrow);
  narrow_reduce(&result->z, &z3_narrow);
}

void readd_to_projective(
  projective_pt_wide_t *result,
  const extended_pt_readd_narrow_t * __restrict x) {

  widen(&result->x, &x->x);
  widen(&result->y, &x->y);
  widen(&result->z, &x->z);
}

void affine_readd_to_extended(
  extended_pt_wide_t *result,
  const extended_affine_pt_readd_narrow_t * __restrict x) {

  widen(&result->x, &x->x);
  widen(&result->y, &x->y);
  mul_narrow(&result->t, &x->x, &x->y);
  for (int i = 0; i < NLIMBS; ++i) {
    result->z.limbs[i] = 0;
  }
  result->z.limbs[1] = 1;
}


#include <stdio.h>
static void print_scalar(const scalar_t *x) {
  printf("[");
  for (int i = 0; i < SCALAR_LIMBS; ++i) {
    printf(" %#x,", x->limbs[i]);
    // printf("x[%d]: %d\n", i, x[i]);
  }
  printf(" ]");
  printf("\n");
}
static void print_projective_pt_wide(const projective_pt_wide_t *x) {
  printf("x: [");
  for (int i = 0; i < NLIMBS; ++i) {
    printf(" %#lx,", x->x.limbs[i]);
    // printf("x[%d]: %d\n", i, x[i]);
  }
  printf(" ]");
  printf("\n");
  printf("y: [");
  for (int i = 0; i < NLIMBS; ++i) {
    printf(" %#lx,", x->y.limbs[i]);
    // printf("x[%d]: %d\n", i, x[i]);
  }
  printf(" ]");
  printf("\n");
  printf("z: [");
  for (int i = 0; i < NLIMBS; ++i) {
    printf(" %#lx,", x->z.limbs[i]);
    // printf("x[%d]: %d\n", i, x[i]);
  }
  printf(" ]");
  printf("\n");
}

static void print_extended_pt_wide(const extended_pt_wide_t *x) {
  printf("x: [");
  for (int i = 0; i < NLIMBS; ++i) {
    printf(" %#lx,", x->x.limbs[i]);
    // printf("x[%d]: %d\n", i, x[i]);
  }
  printf(" ]");
  printf("\n");
  printf("y: [");
  for (int i = 0; i < NLIMBS; ++i) {
    printf(" %#lx,", x->y.limbs[i]);
    // printf("x[%d]: %d\n", i, x[i]);
  }
  printf(" ]");
  printf("\n");
  printf("t: [");
  for (int i = 0; i < NLIMBS; ++i) {
    printf(" %#lx,", x->t.limbs[i]);
    // printf("x[%d]: %d\n", i, x[i]);
  }
  printf(" ]");
  printf("\n");
  printf("z: [");
  for (int i = 0; i < NLIMBS; ++i) {
    printf(" %#lx,", x->z.limbs[i]);
    // printf("x[%d]: %d\n", i, x[i]);
  }
  printf(" ]");
  printf("\n");
}

static void print_readd_reduced(const extended_pt_readd_narrow_reduced_t *x) {
  printf("x: [");
  for (int i = 0; i < NLIMBS_REDUCED; ++i) {
    printf(" %#x,", x->x.limbs[i]);
    // printf("x[%d]: %d\n", i, x[i]);
  }
  printf(" ]");
  printf("\n");
  printf("y: [");
  for (int i = 0; i < NLIMBS_REDUCED; ++i) {
    printf(" %#x,", x->y.limbs[i]);
    // printf("x[%d]: %d\n", i, x[i]);
  }
  printf(" ]");
  printf("\n");
  printf("dt: [");
  for (int i = 0; i < NLIMBS_REDUCED; ++i) {
    printf(" %#x,", x->dt.limbs[i]);
    // printf("x[%d]: %d\n", i, x[i]);
  }
  printf(" ]");
  printf("\n");
  printf("z: [");
  for (int i = 0; i < NLIMBS_REDUCED; ++i) {
    printf(" %#x,", x->z.limbs[i]);
    // printf("x[%d]: %d\n", i, x[i]);
  }
  printf(" ]");
  printf("\n");
}

static void print_readd(const extended_pt_readd_narrow_t *x) {
  printf("x: [");
  for (int i = 0; i < NLIMBS; ++i) {
    printf(" %#x,", x->x.limbs[i]);
    // printf("x[%d]: %d\n", i, x[i]);
  }
  printf(" ]");
  printf("\n");
  printf("y: [");
  for (int i = 0; i < NLIMBS; ++i) {
    printf(" %#x,", x->y.limbs[i]);
    // printf("x[%d]: %d\n", i, x[i]);
  }
  printf(" ]");
  printf("\n");
  printf("dt: [");
  for (int i = 0; i < NLIMBS; ++i) {
    printf(" %#x,", x->dt.limbs[i]);
    // printf("x[%d]: %d\n", i, x[i]);
  }
  printf(" ]");
  printf("\n");
  printf("z: [");
  for (int i = 0; i < NLIMBS; ++i) {
    printf(" %#x,", x->z.limbs[i]);
    // printf("x[%d]: %d\n", i, x[i]);
  }
  printf(" ]");
  printf("\n");
}

void scalar_multiply(
  projective_pt_wide_t *result, const affine_pt_narrow_reduced_t * __restrict x,
  const scalar_t * __restrict n) {

  scalar_t sabs_n;
  convert_to_sabs(&sabs_n, n);

  const int WINDOW_BITS = 5;
  const uint32_t WINDOW_MASK = (1 << WINDOW_BITS) - 1;
  const uint32_t LOOKUP_MASK = WINDOW_MASK >> 1;
  const int TABLE_SIZE = 16;
  extended_pt_readd_narrow_reduced_t table[TABLE_SIZE];

  extended_pt_wide_t x2;
  affine_double_extended(&x2, x);
  affine_to_readd_narrow_reduced(&table[0], x);
  for (int i = 1; i < TABLE_SIZE; ++i) {
    extended_readd_readd_narrow_reduce(&table[i], &x2, &table[i-1]);
  }

  int i;
  int first = 1;
  // Set i to the highest i such that
  // a) i < SCALAR_BITS
  // b) i % WINDOW_BITS = 0

  projective_pt_wide_t temp;
  extended_pt_wide_t temp_ext;
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

    constant_time_extended_narrow_reduced_lookup(
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

  copy_projective_pt_wide(result, &temp);
  explicit_bzero(&sabs_n, sizeof(sabs_n));
  explicit_bzero(&window_pt, sizeof(window_pt));
  explicit_bzero(table, sizeof(table));
  explicit_bzero(&temp, sizeof(temp));
  explicit_bzero(&temp_ext, sizeof(temp_ext));
}
