#ifndef SIGN_H
#define SIGN_H
#include "curve.h"
#include "scalar.h"

#define SIG_LENGTH 65

typedef struct signature {
  residue_narrow_reduced_t y;
  scalar_t s;
} signature_t;

void sign(signature_t *result, scalar_t *priv_key,
  const uint8_t *pub_key, const uint8_t *msg, size_t msg_len);

int verify(
  const signature_t *sig, const uint8_t *r_bytes, const uint8_t *pub_key_bytes,
  const affine_pt_narrow_t *pub_key_pt, const uint8_t *msg,
  size_t msg_len);

void encode_sig(uint8_t *result, const signature_t *sig);
void decode_sig(signature_t *result, const uint8_t *encoded_sig);
#endif
