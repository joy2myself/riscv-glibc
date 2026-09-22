/* Single-precision source-original SLEEF exp10f skeleton.

   This is a target-owned translation of SLEEF xexp10f (u10).  The source
   reduction, six-term polynomial, double-float reconstruction, and split
   exponent scaling are kept in source order.  Clean-only EH quarantine and
   endpoint repair are intentionally not imported here.  */

#include <vecmath_config.h>
#include <v_math.h>

#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

static inline VFLOAT
RVV_SLEEF_ORIG_EXP10F_VPOW2I (VINT e, size_t vlen)
{
  VINT biased = __riscv_vadd (e, EXP_BIAS_F, vlen);
  return I_AS_F (__riscv_vsll (biased, MAN_LEN_F, vlen));
}

static inline VFLOAT
RVV_SLEEF_ORIG_EXP10F_VLDEXP2 (VFLOAT d, VINT e, size_t vlen)
{
  VINT e_half = __riscv_vsra (e, 1, vlen);
  VINT e_rest = __riscv_vsub (e, e_half, vlen);
  VFLOAT p_half = RVV_SLEEF_ORIG_EXP10F_VPOW2I (e_half, vlen);
  VFLOAT p_rest = RVV_SLEEF_ORIG_EXP10F_VPOW2I (e_rest, vlen);
  return VFMUL_VV (VFMUL_VV (d, p_half, vlen), p_rest, vlen);
}

/* Keep the translated xexp10f kernel reusable by source-component
   composites such as exp10m1f.  */
static VFLOAT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_exp10f_kernel_m, LMUL)
  (VFLOAT d, size_t vlen)
{
  /* Pinned SLEEF spells this as vmul(d, vcast(LOG10_2)), not as the
     scalar-vector form.  Keep the broadcast explicit so the target-owned
     source port retains the source operation shape under dynamic FRM.  */
  VFLOAT scaled = VFMUL_VV
    (d, VFMV_VF (RVV_SLEEF_ORIG_F64_TO_F32 (RVV_EXP10F_SOURCE_LOG10_2),
                 vlen),
     vlen);
  VFLOAT u = __riscv_vfcvt_f
    (RVV_SLEEF_ORIG_EXP10F_VFCVT_X_RNE (scaled, vlen), vlen);
  VINT q = RVV_SLEEF_ORIG_EXP10F_VFCVT_X_RNE (u, vlen);

  VFLOAT s = VFMADD_VV (
    u, VFMV_VF (-RVV_EXP10F_SOURCE_L10U_F, vlen), d, vlen);
  s = VFMADD_VV (
    u, VFMV_VF (-RVV_EXP10F_SOURCE_L10L_F, vlen), s, vlen);

  u = VFMV_VF
    (RVV_SLEEF_ORIG_F64_TO_F32 (RVV_EXP10F_SOURCE_POLY0_D), vlen);
  u = VFMADD_VV
    (u, s, VFMV_VF
      (RVV_SLEEF_ORIG_F64_TO_F32 (RVV_EXP10F_SOURCE_POLY1_D), vlen), vlen);
  u = VFMADD_VV
    (u, s, VFMV_VF
      (RVV_SLEEF_ORIG_F64_TO_F32 (RVV_EXP10F_SOURCE_POLY2_D), vlen), vlen);
  u = VFMADD_VV
    (u, s, VFMV_VF
      (RVV_SLEEF_ORIG_F64_TO_F32 (RVV_EXP10F_SOURCE_POLY3_D), vlen), vlen);
  u = VFMADD_VV
    (u, s, VFMV_VF
      (RVV_SLEEF_ORIG_F64_TO_F32 (RVV_EXP10F_SOURCE_POLY4_D), vlen), vlen);
  u = VFMADD_VV
    (u, s, VFMV_VF
      (RVV_SLEEF_ORIG_F64_TO_F32 (RVV_EXP10F_SOURCE_POLY5_D), vlen), vlen);

  VFLOAT x_hi, x_lo;
  RVV_SLEEF_ORIG_DFADD_VF2_VF2_VF (
    VFMV_VF (RVV_SLEEF_ORIG_F64_TO_F32 (RVV_EXP10F_SOURCE_LN10_HI), vlen),
    VFMV_VF (RVV_SLEEF_ORIG_F64_TO_F32 (RVV_EXP10F_SOURCE_LN10_LO), vlen),
    VFMUL_VV (u, s, vlen), x_hi, x_lo, vlen);

  VFLOAT t_hi, t_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF (x_hi, x_lo, s, t_hi, t_lo, vlen);

  VFLOAT sum_hi, sum_lo;
  RVV_SLEEF_ORIG_DFADD_VF_VF2 (VFMV_VF (1.0f, vlen), t_hi, t_lo,
                               sum_hi, sum_lo, vlen);
  RVV_SLEEF_ORIG_DFNORMALIZE_VF2 (sum_hi, sum_lo, sum_hi, sum_lo, vlen);
  u = sum_hi;
  u = RVV_SLEEF_ORIG_EXP10F_VLDEXP2 (u, q, vlen);

  VBOOL overflow = __riscv_vmfgt (d, RVV_EXP10F_SOURCE_MAX_INPUT_F, vlen);
  VBOOL underflow = __riscv_vmflt (d, RVV_EXP10F_SOURCE_MIN_INPUT_F, vlen);
  u = VFMERGE_VFM (u, fp_posInf_f, overflow, vlen);
  u = VFMERGE_VFM (u, fp_posZero_f, underflow, vlen);

  return u;
}

#ifndef RVV_SLEEF_ORIG_EXP10F_HELPERS_ONLY
static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT d, size_t vlen)
{
  return RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_exp10f_kernel_m, LMUL) (d, vlen);
}
#endif

#ifndef RVV_SLEEF_ORIG_EXP10F_HELPERS_ONLY
#undef LMUL
#undef SIMDLEN
#endif
