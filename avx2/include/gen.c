#include <stdlib.h>
#include <string.h>
#include "comb.h"
#include "curve.h"
#include "gen.h"
#include "scalar.h"

void gen_key(scalar_t * __restrict priv_key,
             affine_pt_narrow_t * __restrict pub_key) {
  scalar_hash_t large_key;
  char *large_key_ptr = (char *) &large_key;
  arc4random_buf(large_key_ptr, sizeof(large_key));

  // It's just as random to use montgomery reduction as to correct for the
  // montgomery factor.
  mont_reduce_hash_mod_l(priv_key, &large_key);

  projective_pt_wide_t result_pt;
  scalar_comb_multiply(&result_pt, &base_comb, priv_key);

  residue_wide_t z_inv;

  invert_wide(&z_inv, &result_pt.z);
  mul_wide(&result_pt.x, &result_pt.x, &z_inv);
  mul_wide(&result_pt.y, &result_pt.y, &z_inv);

  residue_narrow_t temp_narrow;
  narrow(&pub_key->x, &result_pt.x);

  narrow(&pub_key->y, &result_pt.y);

  explicit_bzero(&large_key, sizeof(large_key));
  explicit_bzero(&result_pt, sizeof(result_pt));
  explicit_bzero(&z_inv, sizeof(z_inv));
  explicit_bzero(&temp_narrow, sizeof(temp_narrow));
}

void encode_pub_key(uint8_t *result, const affine_pt_narrow_t *pub_key) {
  residue_narrow_reduced_t y_reduced;
  residue_narrow_reduced_t x_reduced;
  narrow_complete(&y_reduced, &pub_key->y);
  narrow_partial_complete(&x_reduced, &pub_key->x);

  y_reduced.limbs[NLIMBS_REDUCED - 1] |= is_odd(&x_reduced) << TBITS;
  encode(result, &y_reduced);
}

int decode_pub_key(affine_pt_narrow_t *result, const uint8_t *encoded_key) {
  residue_narrow_reduced_t y_decoded;
  decode(&y_decoded, encoded_key);
  int is_odd = y_decoded.limbs[NLIMBS_REDUCED - 1] >> TBITS;
  y_decoded.limbs[NLIMBS_REDUCED - 1] &= TMASK;
  return point_decompress(result, &y_decoded, is_odd);
}
