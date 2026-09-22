/* Single-precision SLEEF-component composite for log10p1f.

   SLEEF has xlog10f and xlog1pf, but no native whole-routine log10p1f.  The
   comparison contract used here is deliberately the simple source-shaped
   composition: round dp1 = x + 1 in binary32, then run the target-owned
   xlog10f translation on dp1.  This makes the rounded-1+x and -0 behavior
   visible in comparison results instead of silently claiming a direct port.

   Copyright (C) 2026 Free Software Foundation, Inc.  */

#include <vecmath_config.h>
#include <v_math.h>

#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

#define RVV_SLEEF_ORIG_LOG10F_HELPERS_ONLY
#define RVV_SLEEF_ORIG_LOG10F_SOURCE_KERNEL_ONLY
#include "v_f_log10f_original_skeleton.c"
#undef RVV_SLEEF_ORIG_LOG10F_SOURCE_KERNEL_ONLY
#undef RVV_SLEEF_ORIG_LOG10F_HELPERS_ONLY

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  /* Composite boundary: this is a binary32 add, as it would be at the
     caller boundary before invoking the source xlog10f component.  */
  VFLOAT dp1 = VFADD_VF (x, 1.0f, vlen);
  VFLOAT result = RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_log10f_source_kernel_m, LMUL) (dp1, vlen);

  return result;
}

#undef LMUL
#undef SIMDLEN
