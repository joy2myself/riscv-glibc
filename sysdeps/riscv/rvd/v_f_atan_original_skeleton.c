/* Single-precision source-original SLEEF atan-family skeleton.

   The finite path follows SLEEF xatanf_u1.  The *_original wrappers are
   comparison-only symbols and deliberately do not inherit clean's early EH
   quarantine, fallback, or endpoint repair layer.  */

#include <vecmath_config.h>
#include <v_math.h>

#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT x_orig = x;
  VFLOAT r_hi, r_lo;
  /* xatanf passes (|d|, 0) as y and (1, 0) as x. */
  RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_atan2kf_u1_m, LMUL)
    (__riscv_vfabs (x_orig, vlen), VFMV_VF (0.0f, vlen),
     VFMV_VF (1.0f, vlen), VFMV_VF (0.0f, vlen),
     &r_hi, &r_lo, vlen);

#if defined(COMPILE_FOR_ATANPI_ORIGINAL)
  /* atanpif is a target-side normalized-pi composition.  Preserve the
     source double-float angle until after the split 1/pi multiply; collapsing
     first loses two ULPs for a small ordinary finite subset. */
  VFLOAT scaled_hi, scaled_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
     (r_hi, r_lo,
     VFMV_VF (RVV_ATANF_INV_PI_HI_F, vlen),
     VFMV_VF (RVV_ATANF_INV_PI_LO_F, vlen),
     scaled_hi, scaled_lo, vlen);
  VFLOAT r = VFADD_VV (scaled_hi, scaled_lo, vlen);
#else
  VFLOAT r = VFADD_VV (r_hi, r_lo, vlen);
#endif
  VUINT abs_bits = __riscv_vand (F_AS_U (x_orig), 0x7fffffff, vlen);
  VBOOL x_inf = __riscv_vmseq (abs_bits, 0x7f800000, vlen);
#if defined(COMPILE_FOR_ATANPI_ORIGINAL)
  r = __riscv_vmerge (r, VFMV_VF (RVV_ATANF_NORM_PI2_F, vlen), x_inf, vlen);
#else
  r = __riscv_vmerge (r, VFMV_VF (RVV_ATANF_PI2_F, vlen),
                      x_inf, vlen);
#endif
  r = __riscv_vfsgnjx (r, x_orig, vlen);

  return r;
}

#undef LMUL
#undef SIMDLEN
