/* Single-precision source-original SLEEF exp2f skeleton.

   The finite kernel is a target-owned translation of SLEEF xexp2f.  It keeps
   the source reduction, polynomial, and vldexp2 operation order.  It is a
   comparison implementation: it intentionally does not import clean's
   special-lane quarantine or endpoint repair.  */

#include <vecmath_config.h>
#include <v_math.h>

#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

/* SLEEF's vrint_vf_vf and vrint_vi2_vf request RNE at the conversion
   points only; ordinary floating-point arithmetic retains caller FRM.  */
static inline VFLOAT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_exp2f_vpow2i_m, LMUL)
  (VINT e, size_t vlen)
{
  VINT biased = __riscv_vadd (e, EXP_BIAS_F, vlen);
  return I_AS_F (__riscv_vsll (biased, MAN_LEN_F, vlen));
}

static inline VFLOAT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_exp2f_vldexp2_m, LMUL)
  (VFLOAT d, VINT e, size_t vlen)
{
  /* Keep SLEEF vldexp2's split exponent and two multiply operations.  */
  VINT e_half = __riscv_vsra (e, 1, vlen);
  VINT e_rest = __riscv_vsub (e, e_half, vlen);
  VFLOAT p_half = RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_exp2f_vpow2i_m, LMUL) (e_half, vlen);
  VFLOAT p_rest = RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_exp2f_vpow2i_m, LMUL) (e_rest, vlen);
  return VFMUL_VV (VFMUL_VV (d, p_half, vlen), p_rest, vlen);
}

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT d, size_t vlen)
{
  VFLOAT u = __riscv_vfcvt_f (RVV_VFCVT_X_F_RNE (d, vlen), vlen);
  VINT q = RVV_VFCVT_X_F_RNE (u, vlen);

  VFLOAT s = VFSUB_VV (d, u, vlen);

  /* Keep the source's unsuffixed decimal literals as double data and convert
     each one at the source operation point under the caller's FRM.  */
  u = VFMV_VF (RVV_SLEEF_ORIG_F64_TO_F32 (RVV_EXP2F_SOURCE_POLY0_D), vlen);
  u = VFMADD_VV (u, s, VFMV_VF
    (RVV_SLEEF_ORIG_F64_TO_F32 (RVV_EXP2F_SOURCE_POLY1_D), vlen), vlen);
  u = VFMADD_VV (u, s, VFMV_VF
    (RVV_SLEEF_ORIG_F64_TO_F32 (RVV_EXP2F_SOURCE_POLY2_D), vlen), vlen);
  u = VFMADD_VV (u, s, VFMV_VF
    (RVV_SLEEF_ORIG_F64_TO_F32 (RVV_EXP2F_SOURCE_POLY3_D), vlen), vlen);
  u = VFMADD_VV (u, s, VFMV_VF
    (RVV_SLEEF_ORIG_F64_TO_F32 (RVV_EXP2F_SOURCE_POLY4_D), vlen), vlen);
  u = VFMADD_VV (u, s, VFMV_VF
    (RVV_SLEEF_ORIG_F64_TO_F32 (RVV_EXP2F_SOURCE_POLY5_D), vlen), vlen);
  u = VFMADD_VV (u, s, VFMV_VF (1, vlen), vlen);

  u = RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_exp2f_vldexp2_m, LMUL) (u, q, vlen);

  VBOOL overflow = __riscv_vmfge (d, RVV_EXP2F_SOURCE_MAX_INPUT_F, vlen);
  VBOOL underflow = __riscv_vmflt (d, RVV_EXP2F_SOURCE_MIN_INPUT_F, vlen);
  u = VFMERGE_VFM (u, fp_posInf_f, overflow, vlen);
  u = VFMERGE_VFM (u, fp_posZero_f, underflow, vlen);

  return u;
}

#undef LMUL
#undef SIMDLEN
