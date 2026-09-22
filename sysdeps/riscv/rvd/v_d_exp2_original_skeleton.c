/* Double-precision source-original SLEEF exp2 skeleton.

   This is a target-owned translation of SLEEF xexp2.  It keeps the source
   RNE reduction, Estrin polynomial, two-step vldexp2 reconstruction, and
   source range-fixup order separate from the clean implementation's EH and
   FAST_LDEXP adaptations.  */

#include <vecmath_config.h>
#include <v_math.h>

#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (rvv_sleef_orig_exp2_vpow2i, LMUL)
  (VINT q, size_t vlen)
{
  /* SLEEF vpow2i constructs an exact power-of-two exponent field.  */
  VINT biased = __riscv_vadd (q, EXP_BIAS, vlen);
  return I_AS_F (__riscv_vsll (biased, MAN_LEN, vlen));
}

static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (rvv_sleef_orig_exp2_vldexp2, LMUL)
  (VFLOAT d, VINT q, size_t vlen)
{
  /* Keep source vldexp2's q>>1 / q-(q>>1) split and multiply order.  */
  VINT q_half = __riscv_vsra (q, 1, vlen);
  VINT q_rest = __riscv_vsub (q, q_half, vlen);
  VFLOAT p_half = RVV_LIBMVEC_IMPL_NAME
    (rvv_sleef_orig_exp2_vpow2i, LMUL) (q_half, vlen);
  VFLOAT p_rest = RVV_LIBMVEC_IMPL_NAME
    (rvv_sleef_orig_exp2_vpow2i, LMUL) (q_rest, vlen);
  return VFMUL_VV (VFMUL_VV (d, p_half, vlen), p_rest, vlen);
}

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT d, size_t vlen)
{
  /* Source xexp2 performs a floating RNE first, then an integer RNE.  */
  VINT rounded_once = RVV_VFCVT_X_F_RNE (d, vlen);
  VFLOAT u = __riscv_vfcvt_f (rounded_once, vlen);
  VINT q = RVV_VFCVT_X_F_RNE (u, vlen);
  VFLOAT s = VFSUB_VV (d, u, vlen);

  VFLOAT s2 = VFMUL_VV (s, s, vlen);
  VFLOAT s4 = VFMUL_VV (s2, s2, vlen);
  VFLOAT s8 = VFMUL_VV (s4, s4, vlen);

  /* POLY10(s,s2,s4,s8,c9,...,c0), expanded with source Estrin grouping.  */
  VFLOAT p_hi = VFMADD_VV (
    s, VFMV_VF (RVV_EXP2_POLY_HI1_D, vlen),
    VFMV_VF (RVV_EXP2_POLY_HI0_D, vlen), vlen);

  VFLOAT p_mid_hi_a = VFMADD_VV (
    s, VFMV_VF (RVV_EXP2_POLY_MID_HI_A1_D, vlen),
    VFMV_VF (RVV_EXP2_POLY_MID_HI_A0_D, vlen), vlen);
  VFLOAT p_mid_hi_b = VFMADD_VV (
    s, VFMV_VF (RVV_EXP2_POLY_MID_HI_B1_D, vlen),
    VFMV_VF (RVV_EXP2_POLY_MID_HI_B0_D, vlen), vlen);
  VFLOAT p_mid_hi = VFMADD_VV (s2, p_mid_hi_a, p_mid_hi_b, vlen);

  VFLOAT p_mid_lo_a = VFMADD_VV (
    s, VFMV_VF (RVV_EXP2_POLY_MID_LO_A1_D, vlen),
    VFMV_VF (RVV_EXP2_POLY_MID_LO_A0_D, vlen), vlen);
  VFLOAT p_mid_lo_b = VFMADD_VV (
    s, VFMV_VF (RVV_EXP2_POLY_MID_LO_B1_D, vlen),
    VFMV_VF (RVV_EXP2_POLY_MID_LO_B0_D, vlen), vlen);
  VFLOAT p_mid_lo = VFMADD_VV (s2, p_mid_lo_a, p_mid_lo_b, vlen);

  VFLOAT poly8 = VFMADD_VV (s4, p_mid_hi, p_mid_lo, vlen);
  VFLOAT poly10 = VFMADD_VV (s8, p_hi, poly8, vlen);
  poly10 = VFMADD_VV (poly10, s,
                       VFMV_VF (RVV_LN2_HI_D, vlen), vlen);
  poly10 = VFMADD_VV (poly10, s, VFMV_VF (fp_posOne, vlen), vlen);

  VFLOAT result = RVV_LIBMVEC_IMPL_NAME
    (rvv_sleef_orig_exp2_vldexp2, LMUL) (poly10, q, vlen);

  /* SLEEF source range fixups are deliberately applied last.  */
  VBOOL overflow = __riscv_vmfge (d, RVV_EXP2_OVERFLOW_BOUND_D, vlen);
  VBOOL underflow = __riscv_vmflt (d, RVV_EXP2_UNDERFLOW_BOUND_D, vlen);
  result = VFMERGE_VFM (result, fp_posInf, overflow, vlen);
  result = VFMERGE_VFM (result, fp_posZero, underflow, vlen);

  return result;
}

#undef LMUL
#undef SIMDLEN
