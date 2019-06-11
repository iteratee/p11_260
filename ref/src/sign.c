#define _DEFAULT_SOURCE
#include <blake2.h>
#include <stdlib.h>
#include <string.h>

#include "comb.h"
#include "curve.h"
#include "f11_260.h"
#include "scalar.h"

#include "sign.h"
void sign(signature_t *result, scalar_t *priv_key,
  const uint8_t *pub_key, const uint8_t *msg, size_t msg_len) {
  blake2b_state hash_ctxt;

  char session_key_wash[16];

  scalar_hash_t scalar_large;
  scalar_t session_key;

  arc4random_buf(session_key_wash, sizeof(session_key_wash));
  blake2b_init_key(&hash_ctxt, 64, session_key_wash, sizeof(session_key_wash));
  blake2b_update(&hash_ctxt, (uint8_t *) priv_key, SCALAR_BYTES);
  blake2b_update(&hash_ctxt, (uint8_t *) msg, msg_len);
  blake2b_final(&hash_ctxt, (uint8_t *) &scalar_large, sizeof(scalar_hash_t));

  reduce_hash_mod_l(&session_key, &scalar_large);

  projective_pt_wide_t result_pt;
  scalar_comb_multiply(&result_pt, &base_comb, &session_key);
  residue_wide_t z_inv;

  invert_wide(&z_inv, &result_pt.z);
  mul_wide(&result_pt.x, &result_pt.x, &z_inv);
  mul_wide(&result_pt.y, &result_pt.y, &z_inv);

  residue_narrow_t temp_narrow;
  narrow(&temp_narrow, &result_pt.y);
  narrow_complete(&result->y, &temp_narrow);

  residue_narrow_reduced_t temp_narrow_reduced;
  narrow(&temp_narrow, &result_pt.x);
  narrow_partial_complete(&temp_narrow_reduced, &temp_narrow);
  result->y.limbs[NLIMBS_REDUCED - 1] |=
      is_odd(&temp_narrow_reduced) << (TBITS);

  uint8_t y_buf[RESIDUE_LENGTH_BYTES];
  encode(y_buf, &result->y);

  blake2b_init(&hash_ctxt, 64);
  blake2b_update(&hash_ctxt, y_buf, RESIDUE_LENGTH_BYTES);
  blake2b_update(&hash_ctxt, pub_key, RESIDUE_LENGTH_BYTES);
  blake2b_update(&hash_ctxt, msg, msg_len);
  blake2b_final(&hash_ctxt, (uint8_t *) &scalar_large, sizeof(scalar_hash_t));

  scalar_t hash_scalar;
  reduce_hash_mod_l(&hash_scalar, &scalar_large);
  mult_mod_l(&hash_scalar, &hash_scalar, priv_key);
  sub_mod_l(&result->s, &session_key, &hash_scalar);

  explicit_bzero(&session_key, sizeof(session_key));
  explicit_bzero(&hash_scalar, sizeof(hash_scalar));
  explicit_bzero(&session_key_wash, sizeof(session_key_wash));
}

int verify(
  const signature_t *sig, const uint8_t *r_bytes, const uint8_t *pub_key_bytes,
  const affine_pt_narrow_reduced_t *pub_key_pt, const uint8_t *msg,
  size_t msg_len) {

  projective_pt_wide_t sB;
  projective_pt_wide_t hA;
  projective_pt_wide_t result_pt;
  residue_narrow_reduced_t result_y;

  scalar_hash_t scalar_large;
  blake2b_state hash_ctxt;
  blake2b_init(&hash_ctxt, 64);
  blake2b_update(&hash_ctxt, r_bytes, RESIDUE_LENGTH_BYTES);
  blake2b_update(&hash_ctxt, pub_key_bytes, RESIDUE_LENGTH_BYTES);
  blake2b_update(&hash_ctxt, msg, msg_len);
  blake2b_final(&hash_ctxt, (uint8_t *) &scalar_large, sizeof(scalar_hash_t));

  scalar_t hash_scalar;
  reduce_hash_mod_l(&hash_scalar, &scalar_large);

  // Can use non-const version for both of these.
  scalar_comb_multiply(&sB, &base_comb, &sig->s);
  scalar_multiply(&hA, pub_key_pt, &hash_scalar);
  projective_add(&result_pt, &sB, &hA);

  // Everything below except the comparison should eventually be in helper
  // functions: Point affinization, and point compression bit-for-bit.
  // Same applies for the signing.
  residue_wide_t z_inv;

  invert_wide(&z_inv, &result_pt.z);
  mul_wide(&result_pt.x, &result_pt.x, &z_inv);
  mul_wide(&result_pt.y, &result_pt.y, &z_inv);

  residue_narrow_t temp_narrow;
  narrow(&temp_narrow, &result_pt.y);
  narrow_complete(&result_y, &temp_narrow);

  residue_narrow_reduced_t temp_narrow_reduced;
  narrow(&temp_narrow, &result_pt.x);
  narrow_partial_complete(&temp_narrow_reduced, &temp_narrow);
  result_y.limbs[NLIMBS_REDUCED - 1] |=
      is_odd(&temp_narrow_reduced) << (TBITS);

  return equal_narrow_reduced(&sig->y, &result_y);
}
