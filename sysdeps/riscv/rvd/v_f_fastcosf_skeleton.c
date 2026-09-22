/* Single-precision vector fastcosf function.
   fastcosf(x) computes the pinned SLEEF xfastcosf_u3500 approximation.

   The finite kernel and private xcosf fallback retain the source algorithm.
   The public entry adds the clean VECLIBM-style classify/prepare/sanitize/
   merge envelope so NaN/Inf and conversion-dangerous lanes do not enter the
   source-shaped finite core.  */

#include <vecmath_config.h>
#include <v_math.h>

#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

/* Reuse the target-owned source xcosf translation without importing the
   source-comparator public fastcosf symbol into this translation unit.  The
   helper is static in this TU and therefore does not couple public dispatch
   or the benchmark symbols to *_original.  The embedded skeleton selects
   helper-only and clean-EH behavior from V_MATH_FUNC in v_math.h.  */
#include "v_f_fastcosf_original_skeleton.c"

/* The source vrint_vf_vf/vint conversion points are explicit RNE.  Keep the
   conversion rounding independent of the caller's FRM; the surrounding
   finite arithmetic remains inside the public FRM envelope.  */
static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (rvv_libmvec_fastcosf_source_core, LMUL)
  (VFLOAT d, size_t vlen)
{
  VFLOAT s, u;
  VINT q;
  float m1_pi_scalar = RVV_FASTCOSF_CLEAN_M1_PI_F;
  float pi_scalar = RVV_FASTCOSF_CLEAN_PI_F;
  float pi_neg_scalar = -pi_scalar;
  float pi_half_scalar = pi_scalar * 0.5f;

  /* Pinned xfastcosf_u3500 reduction.  */
  s = VFMADD_VV (d,
                 VFMV_VF (m1_pi_scalar, vlen),
                 VFMV_VF (-0.5f, vlen), vlen);
  u = __riscv_vfcvt_f
    (RVV_VFCVT_X_F_RNE (s, vlen), vlen);
  q = RVV_VFCVT_X_F_RNE (s, vlen);

  /* Keep source's separate subtraction followed by the fused reduction.  */
  VFLOAT d_reduced = VFSUB_VF (d, pi_half_scalar, vlen);
  d = VFMADD_VV
    (u, VFMV_VF (pi_neg_scalar, vlen), d_reduced, vlen);

  s = VFMUL_VV (d, d, vlen);

  /* Pinned xfastcosf_u3500 reduced-precision polynomial.  */
  u = VFMV_VF (RVV_FASTCOSF_CLEAN_POLY0_F, vlen);
  u = VFMADD_VV (u, s, VFMV_VF (RVV_FASTCOSF_CLEAN_POLY1_F, vlen), vlen);
  u = VFMADD_VV (u, s, VFMV_VF (RVV_FASTCOSF_CLEAN_POLY2_F, vlen), vlen);
  u = VFMADD_VV (VFMUL_VV (s, d, vlen), u, d, vlen);

  /* Source flips the result when q is even.  */
  VBOOL even = __riscv_vmseq (__riscv_vand (q, 1, vlen), 0, vlen);
  VUINT sign_mask = __riscv_vmerge
    (VMVU_VX (0, vlen), VMVU_VX (0x80000000u, vlen), even, vlen);
  return U_AS_F (__riscv_vxor (F_AS_U (u), sign_mask, vlen));
}

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT t = x;

  SET_ROUNDTONEAREST;

  /* Classify before any source conversion.  Ordered comparisons and the
     source-shaped conversion core only see finite, conversion-safe lanes.  */
  VUINT vclass = __riscv_vfclass (t, vlen);
  VBOOL NaN_Inf;
  IDENTIFY_F (vclass, class_NaN | class_Inf, NaN_Inf, vlen);

  VBOOL finite = __riscv_vmnot (NaN_Inf, vlen);
  VBOOL fast_range = __riscv_vmflt
    (__riscv_vfabs (t, vlen), RVV_FASTCOSF_FAST_RANGE_F, vlen);
  VBOOL fast_args = __riscv_vmand (finite, fast_range, vlen);
  VBOOL fallback_args = __riscv_vmand
    (finite, __riscv_vmnot (fast_range, vlen), vlen);

  VBOOL not_fast = __riscv_vmnot (fast_args, vlen);
  VBOOL not_fallback = __riscv_vmnot (fallback_args, vlen);
  VBOOL not_special = __riscv_vmnot (NaN_Inf, vlen);

  /* Prepare independent providers with harmless zero inputs on lanes that
     they will not own.  This keeps conversion/fallback work quarantined and
     leaves all selection to the final merge.  */
  VFLOAT fast_input = VFMERGE_VFM (t, fp_posZero_f, not_fast, vlen);
  VFLOAT fallback_input = VFMERGE_VFM (t, fp_posZero_f, not_fallback, vlen);
  VFLOAT special_input = VFMERGE_VFM (t, fp_posZero_f, not_special, vlen);

  VFLOAT fast_result = VFMV_VF (0.0f, vlen);
  VFLOAT fallback_result = VFMV_VF (0.0f, vlen);
  VFLOAT special_result = VFMV_VF (0.0f, vlen);

  UINT nb_fast = __riscv_vcpop (fast_args, vlen);
  UINT nb_fallback = __riscv_vcpop (fallback_args, vlen);
  UINT nb_special = __riscv_vcpop (NaN_Inf, vlen);

  if (nb_fast > 0)
    fast_result = RVV_LIBMVEC_IMPL_NAME
      (rvv_libmvec_fastcosf_source_core, LMUL) (fast_input, vlen);

  if (nb_fallback > 0)
    fallback_result = RVV_LIBMVEC_IMPL_NAME
      (rvv_sleef_orig_fastcosf_xcosf, LMUL)
      (fallback_input, vlen);

  /* Reproduce the clean cosf special-result provider locally.  Ordinary
     finite fallback lanes still use the private source xcosf algorithm above;
     no TU-local clean cosf implementation is called across object files.  */
  if (nb_special > 0)
    {
      VBOOL inf_args;
      IDENTIFY_F (vclass, class_Inf, inf_args, vlen);
      VFLOAT special_arg = __riscv_vfmerge
        (special_input, fp_sNaN_f, inf_args, vlen);
      special_result = __riscv_vfadd
        (NaN_Inf, special_arg, special_arg, vlen);
    }

  VFLOAT result = VMERGE_VVM_F
    (fallback_result, fast_result, fast_args, vlen);
  result = VMERGE_VVM_F (result, special_result, NaN_Inf, vlen);

  RESTORE_FRM;
  return result;
}

#undef LMUL
#undef SIMDLEN
