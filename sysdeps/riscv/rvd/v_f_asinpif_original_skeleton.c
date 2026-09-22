/* Single-precision SLEEF-component composite for asinpif.

   SLEEF supplies xasinf_u1, but not a native asinpif routine.  This
   comparator runs the target-owned xasinf_u1 translation and applies a
   target-side split 1/pi multiplication.  It deliberately does not import
   the clean implementation's EH quarantine or endpoint repair.

   Copyright (C) 2026 Free Software Foundation, Inc.  */

#include <vecmath_config.h>
#include <v_math.h>

#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

#include "v_f_asinf_original_skeleton.c"

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT asin_value = RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_asinf_kernel_m, LMUL) (x, vlen);

  VFLOAT product_hi, product_lo;
  RVV_SLEEF_ORIG_DFMUL_VF_VF
    (asin_value, VFMV_VF (RVV_ASINPIF_ONE_OVER_PI_HI_F, vlen),
     product_hi, product_lo, vlen);
  product_lo = VFADD_VV
    (product_lo,
     VFMUL_VF (asin_value, RVV_ASINPIF_ONE_OVER_PI_LO_F, vlen),
     vlen);
  VFLOAT result = VFADD_VV (product_hi, product_lo, vlen);
  /* The source component preserves the sign of its input.  The positive
     normalized-1/pi composite must preserve it too; in particular, the
     split-product residual can otherwise turn asin(-0) into +0 when the
     pair is collapsed.  */
  return __riscv_vfsgnj (result, x, vlen);
}

#undef LMUL
#undef SIMDLEN
