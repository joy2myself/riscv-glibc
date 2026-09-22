/* Single-precision SLEEF-component composite for log2p1f.

   SLEEF has xlog2f and xlog1pf, but no native whole-routine log2p1f.  The
   comparator therefore rounds dp1 = x + 1 in binary32 and invokes the
   target-owned xlog2f translation.  The rounded-input contract, including
   its tiny-input and -0 consequences, is intentionally explicit.

   Copyright (C) 2026 Free Software Foundation, Inc.  */

#include <vecmath_config.h>
#include <v_math.h>

#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

#include "v_f_log2f_original_skeleton.c"

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT dp1 = VFADD_VF (x, 1.0f, vlen);
  VFLOAT result = RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_log2f_kernel_m, LMUL) (dp1, vlen);

  return result;
}

#undef LMUL
#undef SIMDLEN
