#ifndef GEN_H
#define GEN_H

#include "scalar.h"
#include "curve.h"

void gen_key(scalar_t * __restrict priv_key,
             affine_pt_narrow_t * __restrict pub_key);
void encode_pub_key(uint8_t *result, const affine_pt_narrow_t *pub_key);
int decode_pub_key(affine_pt_narrow_t *result, const uint8_t *encoded_key);
#endif
