/* Single-precision source-original SLEEF log10f skeleton.

   This is a target-owned translation of SLEEF xlog10f.  It keeps the source
   reduction, constants, double-float operation order, and source fixup order
   separate from the clean implementation's EH and precision adaptations.  */

#include <vecmath_config.h>
#include <v_math.h>

#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

static inline VINT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_log10f_vilogb2k_m, LMUL)
  (VFLOAT d, size_t vlen)
{
  VUINT q = __riscv_vsrl (F_AS_U (d), MAN_LEN_F, vlen);
  q = __riscv_vand (q, 0xff, vlen);
  return __riscv_vsub (U_AS_I (q), EXP_BIAS_F, vlen);
}

static inline VFLOAT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_log10f_vldexp3_m, LMUL)
  (VFLOAT d, VINT q, size_t vlen)
{
  /* SLEEF vldexp3 adds the exponent to an already-normalized float.  */
  return I_AS_F (__riscv_vadd (
      F_AS_I (d), __riscv_vsll (q, MAN_LEN_F, vlen), vlen));
}

/* The source and clean entry points share one source-shaped finite graph,
   but use different scalar constant representations.  `_original` performs
   SLEEF's double-to-float casts under the caller's FRM; clean receives
   already-materialized float constants so its VECLIBM EH path does not gain
   source-cast NX side effects.  */

/* Reusable target-owned translation of SLEEF xlog10f.  */
static VFLOAT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_log10f_kernel_with_constants_m, LMUL)
  (VFLOAT d, float poly0, float poly1, float poly2, size_t vlen)
{
  /* SLEEF xlog10f normalizes positive subnormals before extracting e/m.  */
  VBOOL subnormal = __riscv_vmflt (d, RVV_LOG10F_SOURCE_FLT_MIN_F, vlen);
  VFLOAT d_scaled = VFMUL_VF (d, RVV_LOG10F_SOURCE_SUBNORMAL_SCALE_F, vlen);
  d = VMERGE_VVM_F (d, d_scaled, subnormal, vlen);

  /* e = ilogb((4/3)*d), m = d * 2^(-e).  */
  VFLOAT d_adj = VFMUL_VF (d, RVV_LOG10F_SOURCE_INV_075, vlen);
  VINT e = RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_log10f_vilogb2k_m, LMUL) (d_adj, vlen);
  VFLOAT m = RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_log10f_vldexp3_m, LMUL)
    (d, __riscv_vneg (e, vlen), vlen);
  e = VMERGE_VVM_I
    (e, __riscv_vsub (e, RVV_LOG10F_SOURCE_SUBNORMAL_EXP_ADJUST, vlen),
     subnormal, vlen);

  /* x = (m - 1) / (m + 1), retaining the source double-float pair.  */
  VFLOAT num_hi, num_lo;
  VFLOAT den_hi, den_lo;
  RVV_SLEEF_ORIG_DFADD2_VF_VF (VFMV_VF (-1.0f, vlen), m,
                               num_hi, num_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF_VF (VFMV_VF (1.0f, vlen), m,
                               den_hi, den_lo, vlen);

  VFLOAT x_hi, x_lo;
  RVV_SLEEF_ORIG_DFDIV_VF2_VF2 (num_hi, num_lo, den_hi, den_lo,
                                x_hi, x_lo, vlen);
  VFLOAT x2 = VFMUL_VV (x_hi, x_hi, vlen);

  VFLOAT t = VFMV_VF (poly0, vlen);
  t = VFMADD_VV (t, x2,
                 VFMV_VF (poly1, vlen),
                 vlen);
  t = VFMADD_VV (t, x2,
                 VFMV_VF (poly2, vlen),
                 vlen);

  /* s = e*log10(2) + x*(2*log10(e)) + x^3*t.  */
  VFLOAT ef = __riscv_vfcvt_f (e, vlen);
  VFLOAT s_hi, s_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF (
    VFMV_VF (RVV_SLEEF_ORIG_F64_TO_F32 (RVV_LOG10F_SOURCE_LOG10_2_HI), vlen),
    VFMV_VF (RVV_SLEEF_ORIG_F64_TO_F32 (RVV_LOG10F_SOURCE_LOG10_2_LO), vlen),
    ef,
    s_hi, s_lo, vlen);

  VFLOAT x_scale_hi, x_scale_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2 (
    x_hi, x_lo,
    VFMV_VF
      (RVV_SLEEF_ORIG_F64_TO_F32 (RVV_LOG10F_SOURCE_TWO_LOG10_E_HI), vlen),
    VFMV_VF
      (RVV_SLEEF_ORIG_F64_TO_F32 (RVV_LOG10F_SOURCE_TWO_LOG10_E_LO), vlen),
    x_scale_hi, x_scale_lo, vlen);
  RVV_SLEEF_ORIG_DFADD_VF2_VF2 (s_hi, s_lo, x_scale_hi, x_scale_lo,
                                s_hi, s_lo, vlen);

  VFLOAT poly_tail = VFMUL_VV (VFMUL_VV (x2, x_hi, vlen), t, vlen);
  RVV_SLEEF_ORIG_DFADD_VF2_VF2_VF (s_hi, s_lo, poly_tail,
                                   s_hi, s_lo, vlen);

  VFLOAT result = VFADD_VV (s_hi, s_lo, vlen);

  /* SLEEF xlog10f fixup order: +Inf, negative/NaN, then zero.  */
  VUINT vclass = __riscv_vfclass (d, vlen);
  VBOOL pos_inf, nan, zero;
  IDENTIFY_F (vclass, class_Inf, pos_inf, vlen);
  VBOOL neg = __riscv_vmflt (d, fp_posZero_f, vlen);
  IDENTIFY_F (vclass, class_NaN, nan, vlen);
  IDENTIFY_F (vclass, class_Zero, zero, vlen);
  VBOOL invalid = __riscv_vmor (neg, nan, vlen);

  result = VMERGE_VVM_F (result, VFMV_VF (fp_posInf_f, vlen),
                         pos_inf, vlen);
  result = VMERGE_VVM_F (result, VFMV_VF (fp_qNaN_f, vlen),
                         invalid, vlen);
  result = VMERGE_VVM_F (result, VFMV_VF (fp_negInf_f, vlen),
                         zero, vlen);

  return result;
}

/* Source-facing wrapper.  The casts are performed before entering the finite
   graph; all observe the caller's FRM, and the sticky exception state is
   equivalent to the source cast sites.  */
#ifndef RVV_SLEEF_ORIG_LOG10F_CLEAN_KERNEL_ONLY
static VFLOAT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_log10f_source_kernel_m, LMUL)
  (VFLOAT d, size_t vlen)
{
  /* Retain the pre-kernel scalar conversions and their exception effects.
     Their discarded values previously occupied unused structure fields.  */
  (void) RVV_SLEEF_ORIG_F64_TO_F32 (RVV_LOG10F_SOURCE_INV_075);
  float poly0 = RVV_SLEEF_ORIG_F64_TO_F32 (RVV_LOG10F_SOURCE_POLY0_D);
  float poly1 = RVV_SLEEF_ORIG_F64_TO_F32 (RVV_LOG10F_SOURCE_POLY1_D);
  float poly2 = RVV_SLEEF_ORIG_F64_TO_F32 (RVV_LOG10F_SOURCE_POLY2_D);
  (void) RVV_SLEEF_ORIG_F64_TO_F32 (RVV_LOG10F_SOURCE_LOG10_2_HI);
  (void) RVV_SLEEF_ORIG_F64_TO_F32 (RVV_LOG10F_SOURCE_LOG10_2_LO);
  (void) RVV_SLEEF_ORIG_F64_TO_F32 (RVV_LOG10F_SOURCE_TWO_LOG10_E_HI);
  (void) RVV_SLEEF_ORIG_F64_TO_F32 (RVV_LOG10F_SOURCE_TWO_LOG10_E_LO);
  return RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_log10f_kernel_with_constants_m, LMUL)
    (d, poly0, poly1, poly2, vlen);
}
#endif

/* Clean-facing wrapper.  It uses the clean float constants but the same
   source-shaped finite graph, so no source double-to-float conversion is
   introduced into the clean VECLIBM EH path.  */
#ifndef RVV_SLEEF_ORIG_LOG10F_SOURCE_KERNEL_ONLY
static VFLOAT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_log10f_clean_kernel_m, LMUL)
  (VFLOAT d, size_t vlen)
{
  return RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_log10f_kernel_with_constants_m, LMUL)
    (d, RVV_LOG10F_CLEAN_POLY0_F, RVV_LOG10F_CLEAN_POLY1_F,
     RVV_LOG10F_CLEAN_POLY2_F, vlen);
}
#endif

#ifndef RVV_SLEEF_ORIG_LOG10F_HELPERS_ONLY
static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT d, size_t vlen)
{
  return RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_log10f_source_kernel_m, LMUL) (d, vlen);
}
#endif

#ifndef RVV_SLEEF_ORIG_LOG10F_HELPERS_ONLY
#undef LMUL
#undef SIMDLEN
#endif
