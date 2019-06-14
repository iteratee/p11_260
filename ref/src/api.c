#define _DEFAULT_SOURCE
#include <string.h>
#include "crypto_sign.h"
#include "curve.h"
#include "gen.h"
#include "scalar.h"
#include "sign.h"

int crypto_sign_keypair(unsigned char *pk, unsigned char *sk) {
  affine_pt_narrow_t pub_key_pt;
  scalar_t priv_key;
  gen_key(&priv_key, &pub_key_pt);
  encode_pub_key(pk, &pub_key_pt);
  memcpy(sk, &priv_key, SCALAR_BYTES);
  memcpy(sk + SCALAR_BYTES, pk, RESIDUE_LENGTH_BYTES);
  explicit_bzero(&priv_key, sizeof(priv_key));
  return 0;
}

int crypto_sign(
    unsigned char *sm,unsigned long long *smlen,
    const unsigned char *m,unsigned long long mlen,
    const unsigned char *sk) {
  signature_t sig_struct;
  scalar_t priv_key;
  priv_key.limbs[SCALAR_LIMBS - 1] = 0;
  memcpy(&priv_key, sk, SCALAR_BYTES);
  sign(&sig_struct, &priv_key, sk + SCALAR_BYTES);

  *smlen = mlen + SIG_LENGTH;
  encode_sig(sm, &sig_struct);
  memcpy(sm + SIG_LENGTH, m, mlen);
  return 0;
}

int crypto_sign_open(
    unsigned char *m,unsigned long long *mlen,
    const unsigned char *sm,unsigned long long smlen,
    const unsigned char *pk) {
  signature_t sig_struct;
  decode_sig(&sig_struct, sm);
  affine_pt_narrow_t pub_key_pt;
  if (!decode_pub_key(&pub_key_pt, pk)) {
    return -1;
  }
  if (!verify(&sig_struct, sm, pk, &pub_key_pt, sm + SIG_LENGTH,
              smlen - SIG_LENGTH)) {
    return -2;
  }
  *mlen = smlen - SIG_LENGTH;
  memcpy(m, sm, smlen - SIG_LENGTH);
  return 0;
}
