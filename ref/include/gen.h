#ifndef GEN_H
#define GEN_H

#include "scalar.h"
#include "curve.h"

void gen_key(scalar_t * __restrict priv_key,
             affine_pt_narrow_reduced_t * __restrict pub_key);
#endif
