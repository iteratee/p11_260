#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "comb.h"
#include "curve.h"
#include "f11_260.h"
#include "gen.h"
#include "scalar.h"
#include "sign.h"

int main(int _argc, char **argv) {
  residue_narrow_t x = {
    .limbs = {
      0x14e8b6e, 0x3553e74, 0x0464e4c, 0x61de408,
      0x006a30e, 0x6e9b25b, 0x3e6f39e, 0x19ec754,
      0x5c71cc3, 0x2bc1c0e, 0x554338e, 0x14e8b6e,
    },
  };

  residue_wide_t two = {
    .limbs = {
      0x0, 0x2, 0x0, 0x0,
      0x0, 0x0, 0x0, 0x0,
      0x0, 0x0, 0x0, 0x0,
    },
  };

  residue_wide_t x_plus_two;

  residue_narrow_reduced_t x_narrow_reduced = {
    .limbs = {
      0x206b305, 0x2f7c2ce, 0x0cf58a7, 0x2b81791, 0x19b26fa,
      0x2986830, 0x0503be5, 0x0789163, 0x16d90a0, 0x005a82e,
    },
  };

  residue_wide_t x_wide;

  residue_narrow_t y = {
    .limbs = {
      0x56ed38e, 0x5f5b0e1, 0x4668277, 0x0f7d85a,
      0x4515e42, 0x00cb559, 0x3f8a910, 0x6655708,
      0x3085b4d, 0x581ceff, 0x3324c03, 0x56ed38e,
    },
  };

  residue_narrow_reduced_t y_narrow_reduced = {
    .limbs = {
      0x086dd54, 0x2f7aedb, 0x38904ae, 0x2e28aa4, 0x29de1ad,
      0x289d572, 0x0f6837a, 0x19987b1, 0x012fb71, 0x1c37867,
    },
  };

  residue_wide_t y_wide;

  residue_wide_t mul_expected = {
    .limbs = {
      0x06e9e1d, 0x1c508c4, 0x3eeb85d, 0x04bc914,
      0x0a57e1c, 0x1f13f9a, 0x2d8aa7d, 0x232cce3,
      0x31e92c4, 0x04fb073, 0x2582507, 0x06e9e1d,
    },
  };

  residue_wide_t square_expected = {
    .limbs = {
      0x3088d3c, 0x2073353, 0x18e5de4, 0x320a4ab,
      0x3ee123a, 0x2d88419, 0x3d1ae13, 0x02b3dcf,
      0x2997027, 0x3d550a2, 0x220a052, 0x3088d3c,
    },
  };

  residue_narrow_t negative_one_redundant = {
    .limbs = {
      0x000000e, 0x3ffffff, 0x3ffffff, 0x3ffffff,
      0x3ffffff, 0x3ffffff, 0x3ffffff, 0x3ffffff,
      0x3ffffff, 0x3ffffff, 0x3ffffff, 0x000000e,
    },
  };

  residue_narrow_t negative_t2_plus_one = {
    .limbs = {
      0x000000e,
      0x3ffffff, 0x000000e, 0x3ffffff,
      0x3ffffff, 0x3ffffff, 0x3ffffff, 0x3ffffff,
      0x3ffffff, 0x3ffffff, 0x3ffffff, 0x000000e,
    },
  };

  residue_narrow_reduced_t negative_t2_plus_one_partial = {
    .limbs = {
      0x3fffff1, 0x0000000, 0x3fffff1, 0x3fffff1,
      0x3fffff1, 0x3fffff1, 0x3fffff1, 0x3fffff1,
      0x3fffff1, 0x3fffff1,
    },
  };

  residue_narrow_reduced_t negative_t2_plus_one_complete = {
    .limbs = {
      0x0000000, 0x0000001, 0x3fffff1, 0x3fffff1,
      0x3fffff1, 0x3fffff1, 0x3fffff1, 0x3fffff1,
      0x3fffff1, 0x3fffff1,
    },
  };

  residue_wide_t sqrt_x_plus_2_over_y = {
    .limbs = {
      0x040bbb0, 0x3fa8549, 0x0706e5c, 0x3b33dc9,
      0x3401712, 0x3a58fb3, 0x076ec4f, 0x3347ad0,
      0x16ca1b0, 0x26ed559, 0x06033f0, 0x040bbb0,
    },
  };

  residue_wide_t x_inverse = {
    .limbs = {
      0x09fd09b, 0x17a9f53, 0x22e2983, 0x0f09456,
      0x11fb41e, 0x1e47b3f, 0x37dd25f, 0x3bc6938,
      0x2b654cd, 0x233a0b2, 0x3f8c25b, 0x09fd09b,
    },
  };

  #if 1
  residue_wide_t result;
  residue_narrow_t result_narrow;
  residue_narrow_reduced_t result_narrow_reduced;

  mul_narrow(&result, &x, &y);
  for (int i = 0; i < NLIMBS; ++i) {
    assert(mul_expected.limbs[i] == result.limbs[i]);
  }

  widen(&x_wide, &x);
  for (int i = 0; i < NLIMBS; ++i) {
    assert(x.limbs[i] == x_wide.limbs[i]);
  }

  mul_wide_narrow(&result, &x_wide, &y);
  for (int i = 0; i < NLIMBS; ++i) {
    assert(mul_expected.limbs[i] == result.limbs[i]);
  }

  widen(&y_wide, &y);
  mul_wide(&result, &x_wide, &y_wide);
  for (int i = 0; i < NLIMBS; ++i) {
    assert(mul_expected.limbs[i] == result.limbs[i]);
  }

  square_narrow(&result, &x);
  for (int i = 0; i < NLIMBS; ++i) {
    assert(square_expected.limbs[i] == result.limbs[i]);
  }

  square_wide(&result, &x_wide);
  for (int i = 0; i < NLIMBS; ++i) {
    assert(square_expected.limbs[i] == result.limbs[i]);
  }

  // The reduction function doesn't reduce this redundant version of negative
  // one any more.
  reduce_step_narrow(&result_narrow, &negative_one_redundant);
  for (int i = 0; i < NLIMBS; ++i) {
    assert(negative_one_redundant.limbs[i] == result_narrow.limbs[i]);
  }

  reduce_step_narrow(&result_narrow, &negative_t2_plus_one);
  for (int i = 0; i < NLIMBS; ++i) {
    assert(negative_t2_plus_one.limbs[i] == result_narrow.limbs[i]);
  }

  narrow_partial_complete(&result_narrow_reduced, &negative_t2_plus_one);
  for (int i = 0; i < NLIMBS_REDUCED; ++i) {
    assert(negative_t2_plus_one_partial.limbs[i] ==
        result_narrow_reduced.limbs[i]);
  }

  narrow_complete(&result_narrow_reduced, &negative_t2_plus_one);
  for (int i = 0; i < NLIMBS_REDUCED; ++i) {
    assert(negative_t2_plus_one_complete.limbs[i] ==
        result_narrow_reduced.limbs[i]);
  }

  scalar_t scalar_result;
  scalar_t scalar_x = {
    .limbs = {
      0xa46168f9, 0x4cbf07a5, 0x62cf2928, 0xfd04242b, 0x3b12d23f,
      0x355e9e63, 0xc22e849e, 0x6331c34a, 0x1,
    },
  };
  scalar_t scalar_y = {
    .limbs = {
      0x148b9452, 0xaca9b6bb, 0xe0eeb33d, 0x7e64c899, 0xd61c602a,
      0x96dcbb6b, 0x6a037c88, 0x39fbbaf0, 0x0,
    },
  };
  scalar_t scalar_x_plus_y = {
    .limbs = {
      0xb8ecfd4b, 0xf968be60, 0x43bddc65, 0x7b68ecc5, 0x112f326a,
      0xcc3b59cf, 0x2c320126, 0x9d2d7e3b, 0x1,
    },
  };
  scalar_t scalar_x_plus_x = {
    .limbs = {
      0x48c2d1f2, 0x997e0f4b, 0xc59e5250, 0xfa084856, 0x7625a47f,
      0x6abd3cc6, 0x845d093c, 0xc6638695, 0x2,
    },
  };
  scalar_t scalar_x_plus_x_plus_x_plus_y = {
    .limbs = {
      0xd90232fc, 0xac09d5c3, 0xd4a42a06, 0x1a7823b2, 0x2a5e47bb,
      0x24a61ea1, 0xa6cd4ac4, 0x639199d0, 0x0,
    },
  };
  scalar_t scalar_x_minus_y = {
    .limbs = {
      0x8fd5d4a7, 0xa01550ea, 0x81e075ea, 0x7e9f5b91, 0x64f67215,
      0x9e81e2f7, 0x582b0815, 0x2936085a, 0x1,
    },
  };
  scalar_t scalar_y_minus_x = {
    .limbs = {
      0x98d7c79a, 0x46c7a6fd, 0xb2d78ec5, 0xdc59b5d7, 0xf8001d19,
      0x73d094fc, 0xb196b789, 0xd6c962a5, 0x2,
    },
  };
  scalar_t scalar_x_times_y = {
    .limbs = {
      0x30b3d35a, 0x9ca90acf, 0x6926efdd, 0x80620b0a, 0x52e190e7,
      0x8011b9b8, 0x8c7d8f43, 0x90491703, 0x3,
    },
  };
  scalar_t scalar_x_sabs = {
    .limbs = {
      0x80d57bfa, 0x58a59402, 0x47f78b34, 0x488fef43, 0xe39c4ac1,
      0xf60a5f48, 0x4d93c310, 0xb19a0ba5, 0x0,
    },
  };
  scalar_t scalar_y_sabs = {
    .limbs = {
      0x4d415fc7, 0xfc096781, 0x21635296, 0x36bcca2f, 0x5f9c594e,
      0xaff2a9c7, 0x265f1ed5, 0x1cfebcf8, 0x2,
    },
  };
  scalar_hash_t scalar_hash_val = {
    .limbs = {
      0xcbbc3de7, 0xa212405d, 0x5c85f47c, 0x79aa991c,
      0xfe310944, 0x54075530, 0xd5ef6878, 0x72e57186,
      0x36dcac18, 0xb72461e2, 0x5405caca, 0x4e9e0bff,
      0x8d67a990, 0xf62f262c, 0x6df205dd, 0x24d78573,
    },
  };
  scalar_t reduced_hash_val = {
    .limbs = {
      0xef1d4f9d, 0xd832a3a5, 0xdf1682be, 0x8d257e79, 0x41b1f2ca,
      0x5be9564c, 0x320d4cb6, 0x108f8d04, 0x3,
    },
  };

  uint8_t buffer[33];
  uint8_t encode_x[33] = {
    0x05, 0xb3, 0x06, 0x3a, 0x0b, 0xdf,
    0x7b, 0x8a, 0xf5, 0x4c, 0xe4, 0x05, 0xae,
    0xfa, 0x26, 0x9b, 0xc1, 0xa0, 0x61,
    0x5a, 0xbe, 0x03, 0xc5, 0x58, 0x24, 0x1e,
    0xa0, 0x90, 0x6d, 0xb9, 0xa0, 0x16, 0x00,
  };
  uint8_t encode_y[33] = {
    0x54, 0xdd, 0x86, 0x6c, 0xbb, 0xde,
    0xeb, 0x4a, 0x90, 0x38, 0xa9, 0xa2, 0xb8,
    0xad, 0xe1, 0x9d, 0xca, 0x55, 0x27,
    0xaa, 0x37, 0x68, 0x4f, 0xec, 0x61, 0x66,
    0x71, 0xfb, 0x12, 0x9c, 0xe1, 0x0d, 0x07,
  };

  add_mod_l(&scalar_result, &scalar_x, &scalar_y);
  for (int i = 0; i < SCALAR_LIMBS; ++i) {
    assert(scalar_x_plus_y.limbs[i] == scalar_result.limbs[i]);
  }
  add_mod_l(&scalar_result, &scalar_x, &scalar_x);
  for (int i = 0; i < SCALAR_LIMBS; ++i) {
    assert(scalar_x_plus_x.limbs[i] == scalar_result.limbs[i]);
  }
  add_mod_l(&scalar_result, &scalar_x_plus_x, &scalar_x_plus_y);
  for (int i = 0; i < SCALAR_LIMBS; ++i) {
    assert(scalar_x_plus_x_plus_x_plus_y.limbs[i] == scalar_result.limbs[i]);
  }
  sub_mod_l(&scalar_result, &scalar_x, &scalar_y);
  for (int i = 0; i < SCALAR_LIMBS; ++i) {
    assert(scalar_x_minus_y.limbs[i] == scalar_result.limbs[i]);
  }

  sub_mod_l(&scalar_result, &scalar_y, &scalar_x);
  for (int i = 0; i < SCALAR_LIMBS; ++i) {
    assert(scalar_y_minus_x.limbs[i] == scalar_result.limbs[i]);
  }

  mult_mod_l(&scalar_result, &scalar_x, &scalar_y);
  for (int i = 0; i < SCALAR_LIMBS; ++i) {
    assert(scalar_x_times_y.limbs[i] == scalar_result.limbs[i]);
  }

  convert_to_sabs(&scalar_result, &scalar_x);
  for (int i = 0; i < SCALAR_LIMBS; ++i) {
    assert(scalar_x_sabs.limbs[i] == scalar_result.limbs[i]);
  }

  convert_to_sabs(&scalar_result, &scalar_y);
  for (int i = 0; i < SCALAR_LIMBS; ++i) {
    assert(scalar_y_sabs.limbs[i] == scalar_result.limbs[i]);
  }

  encode(buffer, &x_narrow_reduced);
  for (int i = 0; i < 33; ++i) {
    assert(encode_x[i] == buffer[i]);
  }

  encode(buffer, &y_narrow_reduced);
  for (int i = 0; i < 33; ++i) {
    assert(encode_y[i] == buffer[i]);
  }

  decode(&result_narrow_reduced, encode_x);
  for (int i = 0; i < NLIMBS_REDUCED; ++i) {
    assert(x_narrow_reduced.limbs[i] == result_narrow_reduced.limbs[i]);
  }

  decode(&result_narrow_reduced, encode_y);
  for (int i = 0; i < NLIMBS_REDUCED; ++i) {
    assert(y_narrow_reduced.limbs[i] == result_narrow_reduced.limbs[i]);
  }

  //x/y is not a quadratic residue, but (x+2)/y is.
  assert(!sqrt_inv_wide(&result, &x_wide, &y_wide));
  add_wide(&x_plus_two, &x_wide, &two);
  assert(sqrt_inv_wide(&result, &x_plus_two, &y_wide));
  for (int i = 0; i < NLIMBS; ++i) {
    assert(sqrt_x_plus_2_over_y.limbs[i] == result.limbs[i]);
  }

  invert_wide(&result, &x_wide);
  for (int i = 0; i < NLIMBS; ++i) {
    assert(x_inverse.limbs[i] == result.limbs[i]);
  }

  reduce_hash_mod_l(&scalar_result, &scalar_hash_val);
  for (int i = 0; i < SCALAR_LIMBS; ++i) {
    assert(reduced_hash_val.limbs[i] == scalar_result.limbs[i]);
  }

  scalar_t mult_scalar = {
    .limbs = {
      0x55f0b9a3, 0x82b106c5, 0xcb2e2b7d, 0x30735cbc,
      0xa512a8ba, 0x4c5cd391, 0xe9d0c788, 0x92bb2562, 0x3,
    },
  };
  projective_pt_wide_t expected_scalar_mult = {
    .x = {
      .limbs = {
        0x0350abe, 0x1267d8d, 0x39a3cd3, 0x09e1275, 0x2d21378, 0x24771d9,
        0x3558a1d, 0x3bdca9b, 0x0dd862d, 0x0bb230a, 0x1668292, 0x0350abe,
      },
    },
    .y = {
      .limbs = {
        0x0b090d6, 0x04d69fd, 0x03e739d, 0x36ce258, 0x0b6464b, 0x19dab22,
        0x249c1a8, 0x1d28c7d, 0x1591dbc, 0x085ebab, 0x0e8274f, 0x0b090d6,
      },
    },
    .z = {
      .limbs = {0, 0x1},
    },
  };
  projective_pt_wide_t result_pt;

  for (int i = 0; i<1; ++i) {
    scalar_multiply(&result_pt, &B, &mult_scalar);
  }
  {
    residue_wide_t tmp;
    mul_wide(&tmp, &expected_scalar_mult.x, &result_pt.z);
    assert(equal_wide(&tmp, &result_pt.x));
    mul_wide(&tmp, &expected_scalar_mult.y, &result_pt.z);
    assert(equal_wide(&tmp, &result_pt.y));
  }

  affine_pt_narrow_t expected_everything0 = {
    .x = {
      .limbs = {
        0, 0x20eef1a, 0x3c30e66, 0x0d710f0, 0x248a6fa, 0x30c967f,
        0x3ce302c, 0x0ccd1f2, 0x197e993, 0x2ebaef3, 0x0f2f019, 0,
      },
    },
    .y = {
      .limbs = {
        0, 0x3017cc0, 0x02a5110, 0x06d37e5, 0x283a64a, 0x01484b5,
        0x196f37b, 0x13de2d2, 0x0da32d1, 0x392e0fc, 0x221d742, 0,
      },
    },
  };

  affine_pt_narrow_t expected_everything1 = {
    .x = {
      .limbs = {
        0, 0x0e35d45, 0x038f90c, 0x0283483, 0x01ee50a, 0x1e364f9,
        0x362414c, 0x156b1ed, 0x006fff6, 0x271f9ed, 0x0ffa45d, 0,
      },
    },
    .y = {
      .limbs = {
        0, 0x156ae67, 0x27941ab, 0x19a3000, 0x3572ab5, 0x2b90ce3,
        0x136156c, 0x0727496, 0x0edae82, 0x0fa5dfd, 0x16f293c, 0,
      },
    },
  };

  affine_pt_narrow_t expected_everything2 = {
    .x = {
      .limbs = {
        0, 0x37fcb1b, 0x16004b9, 0x1d18743, 0x0bce648, 0x0d78db6,
        0x35b1d65, 0x23bb620, 0x2fbc323, 0x1a9a586, 0x3b22577, 0,
      },
    },
    .y = {
      .limbs = {
        0, 0x082fb15, 0x03487d6, 0x3d1c2c9, 0x2c9e7ad, 0x187be10,
        0x2e9b6ba, 0x15b8f89, 0x243ae4c, 0x328bb11, 0x00b12a9, 0,
      },
    },
  };


  affine_pt_narrow_t expected_everything3 = {
    .x = {
      .limbs = {
        0, 0x3e79b25, 0x2ca71b7, 0x2b2ea3c, 0x0de7ac4, 0x3026d10,
        0x2bce79e, 0x1153866, 0x03e5a80, 0x22b9a37, 0x03e9c59, 0,
      },
    },
    .y = {
      .limbs = {
        0, 0x20100d6, 0x2330974, 0x3402585, 0x172cfd6, 0x275a21c,
        0x213e87c, 0x29989f2, 0x155e437, 0x096a378, 0x3a674eb, 0,
      },
    },
  };

  affine_pt_narrow_t expected_gray_code_end0 = {
    .x = {
      .limbs = {
        0, 0x14dd884, 0x12c9e33, 0x2d42122, 0x26f0b14, 0x1b9ea17,
        0x3779e94, 0x2562a88, 0x0be34f0, 0x192ead9, 0x089ec45, 0,
      },
    },
    .y = {
      .limbs = {
        0, 0x1de5221, 0x172f820, 0x28c1b33, 0x08003c6, 0x0e65926,
        0x188cd49, 0x3bb39fd, 0x1b9d8d7, 0x03d5020, 0x045742b, 0,
      },
    },
  };

  affine_pt_narrow_t expected_gray_code_end1 = {
    .x = {
      .limbs = {
        0, 0x1d1cf29, 0x2e289d7, 0x1a83709, 0x2252d11, 0x3d6411c,
        0x3fd73ad, 0x2737d9c, 0x2ca9eba, 0x058f290, 0x3879a7c, 0,
      },
    },
    .y = {
      .limbs = {
        0, 0x357399d, 0x0276752, 0x0d5199f, 0x1bbd3a0, 0x39044f1,
        0x0c5e83a, 0x1a99cdd, 0x0dcb61f, 0x35b7272, 0x1184cff, 0,
      },
    },
  };

  affine_pt_narrow_t expected_gray_code_end2 = {
    .x = {
      .limbs = {
        0, 0x1ea3c19, 0x081dc9e, 0x1a0b337, 0x1d7f3f4, 0x295a0aa,
        0x1ebff45, 0x0956bf0, 0x17aae80, 0x05d8632, 0x3082c9a, 0,
      },
    },
    .y = {
      .limbs = {
        0, 0x22ad91f, 0x1ffcc65, 0x37b4f5c, 0x29c51ab, 0x3f9bd02,
        0x296aaf9, 0x2a58b82, 0x2c54e16, 0x2a7672c, 0x21486e2, 0,
      },
    },
  };

  affine_pt_narrow_t expected_gray_code_end3 = {
    .x = {
      .limbs = {
        0, 0x06b9c9d, 0x3d00674, 0x10a73fc, 0x30fda83, 0x139185c,
        0x043e082, 0x3c67915, 0x208192a, 0x025e451, 0x258a566, 0,
      },
    },
    .y = {
      .limbs = {
        0, 0x3d2a04f, 0x1314c36, 0x131c7a3, 0x1882ef3, 0x1a0a5e8,
        0x1919356, 0x0a5616a, 0x1eea31d, 0x2c216b3, 0x18ba4aa, 0,
      },
    },
  };

  sabs_comb_set_t computed_base_comb;
  compute_comb_set(&computed_base_comb, &B);
  for (int i = 0; i < NLIMBS_REDUCED; ++i) {
    //assert(computed_base_comb.combs[0].table[COMB_TABLE_SIZE - 1].x.limbs[i] ==
      //expected_everything0.x.limbs[i]);
    //assert(computed_base_comb.combs[0].table[COMB_TABLE_SIZE - 1].y.limbs[i] ==
      //expected_everything0.y.limbs[i]);
    //assert(computed_base_comb.combs[1].table[COMB_TABLE_SIZE - 1].x.limbs[i] ==
      //expected_everything1.x.limbs[i]);
    //assert(computed_base_comb.combs[1].table[COMB_TABLE_SIZE - 1].y.limbs[i] ==
      //expected_everything1.y.limbs[i]);
    //assert(computed_base_comb.combs[2].table[COMB_TABLE_SIZE - 1].x.limbs[i] ==
      //expected_everything2.x.limbs[i]);
    //assert(computed_base_comb.combs[2].table[COMB_TABLE_SIZE - 1].y.limbs[i] ==
      //expected_everything2.y.limbs[i]);
    //assert(computed_base_comb.combs[3].table[COMB_TABLE_SIZE - 1].x.limbs[i] ==
      //expected_everything3.x.limbs[i]);
    //assert(computed_base_comb.combs[3].table[COMB_TABLE_SIZE - 1].y.limbs[i] ==
      //expected_everything3.y.limbs[i]);
  }

  for (int i = 0; i < NLIMBS_REDUCED; ++i) {
    //assert(computed_base_comb.combs[0].table[7].x.limbs[i] ==
      //expected_gray_code_end0.x.limbs[i]);
    //assert(computed_base_comb.combs[0].table[7].y.limbs[i] ==
      //expected_gray_code_end0.y.limbs[i]);
    //assert(computed_base_comb.combs[1].table[7].x.limbs[i] ==
      //expected_gray_code_end1.x.limbs[i]);
    //assert(computed_base_comb.combs[1].table[7].y.limbs[i] ==
      //expected_gray_code_end1.y.limbs[i]);
    //assert(computed_base_comb.combs[2].table[7].x.limbs[i] ==
      //expected_gray_code_end2.x.limbs[i]);
    //assert(computed_base_comb.combs[2].table[7].y.limbs[i] ==
      //expected_gray_code_end2.y.limbs[i]);
    //assert(computed_base_comb.combs[3].table[7].x.limbs[i] ==
      //expected_gray_code_end3.x.limbs[i]);
    //assert(computed_base_comb.combs[3].table[7].y.limbs[i] ==
      //expected_gray_code_end3.y.limbs[i]);
  }
  #endif

  #if 0
  for (int i = 0; i<1; ++i) {
    scalar_comb_multiply(&result_pt, &base_comb, &mult_scalar);
  }
  {
    residue_wide_t tmp;
    mul_wide(&tmp, &expected_scalar_mult.x, &result_pt.z);
    assert(equal_wide(&tmp, &result_pt.x));
    mul_wide(&tmp, &expected_scalar_mult.y, &result_pt.z);
    assert(equal_wide(&tmp, &result_pt.y));
  }
  #endif
  #if 0
  for (int i = 0; i<100000; ++i) {
    scalar_t priv_key;
    affine_pt_narrow_reduced_t pub_key;
    gen_key(&priv_key, &pub_key);
  }
  #endif
  for (int i = 0; i < 100000; ++i) {
    uint8_t encoded_sk[66];
    uint8_t encoded_sig[65];
    const uint8_t *msg = (uint8_t *) "Hello World!";
    const size_t msglen = 13;
    scalar_t priv_key;
    scalar_t priv_key_decoded;
    affine_pt_narrow_t pub_key;
    affine_pt_narrow_t pub_key_decoded;
    gen_key(&priv_key, &pub_key);
    memcpy(encoded_sk, &priv_key, SCALAR_BYTES);
    encode_pub_key(encoded_sk + SCALAR_BYTES, &pub_key);
    priv_key_decoded.limbs[SCALAR_LIMBS - 1] = 0;
    memcpy(&priv_key_decoded, encoded_sk, SCALAR_BYTES);
    for (int j = 0; j < SCALAR_LIMBS; ++j) {
      assert(priv_key.limbs[j] == priv_key_decoded.limbs[j]);
    }
    signature_t result;
    sign(&result, &priv_key_decoded, encoded_sk + SCALAR_BYTES, msg, msglen);
    encode_sig(encoded_sig, &result);
    signature_t result_decoded;
    decode_sig(&result_decoded, encoded_sig);
    for (int j = 0; j < SCALAR_LIMBS; ++j) {
      assert(result.s.limbs[j] == result_decoded.s.limbs[j]);
    }
    for (int j = 0; j < NLIMBS_REDUCED; ++j) {
      assert(result.y.limbs[j] == result_decoded.y.limbs[j]);
    }
    assert(decode_pub_key(&pub_key_decoded, encoded_sk + SCALAR_BYTES));

    uint8_t y_buf[RESIDUE_LENGTH_BYTES];
    encode(y_buf, &result_decoded.y);
    if(!verify(&result, y_buf, encoded_sk + SCALAR_BYTES, &pub_key_decoded, msg,
               msglen)) {
      printf("verification failed\n");
      exit(1);
    }
  }
}
