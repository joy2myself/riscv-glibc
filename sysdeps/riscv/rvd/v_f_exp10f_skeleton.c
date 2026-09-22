/* Single-precision vector exp10 function.
   exp10f(x) computes 10^x

   Based on SLEEF xexp10f (u10 = 1.0 ULP precision).

   Copyright (C) 2026 Free Software Foundation, Inc.
   This file is part of the GNU C Library. */

#include <vecmath_config.h>
#include <v_math.h>

#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

/* Keep SLEEF's split-exponent scaling in the finite kernel.  The clean
   wrapper supplies the EH quarantine below; this helper is otherwise the
   source vldexp2_vf_vf_vi2 operation, including its two ordered multiplies.  */
static inline VFLOAT
RVV_LIBMVEC_EXP10F_VPOW2I (VINT e, size_t vlen)
{
  VINT biased = __riscv_vadd (e, EXP_BIAS_F, vlen);
  return I_AS_F (__riscv_vsll (biased, MAN_LEN_F, vlen));
}

static inline VFLOAT
RVV_LIBMVEC_EXP10F_VLDEXP2 (VFLOAT d, VINT e, size_t vlen)
{
  VINT e_half = __riscv_vsra (e, 1, vlen);
  VINT e_rest = __riscv_vsub (e, e_half, vlen);
  VFLOAT p_half = RVV_LIBMVEC_EXP10F_VPOW2I (e_half, vlen);
  VFLOAT p_rest = RVV_LIBMVEC_EXP10F_VPOW2I (e_rest, vlen);
  return VFMUL_VV (VFMUL_VV (d, p_half, vlen), p_rest, vlen);
}

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT vx = x;
  VFLOAT vy, vy_special;
  VBOOL special_args, range_args, overflow, underflow;

  SET_ROUNDTONEAREST;

  VFLOAT vx_orig = vx;

  /* Handle NaN/Inf and out-of-range inputs (align with v_f_logf_skeleton).  */
  VUINT vclass = __riscv_vfclass (vx, vlen);
  IDENTIFY (vclass, class_NaN | class_Inf, special_args, vlen);
  VFLOAT vx_cmp = __riscv_vfmerge (vx_orig, fp_posZero, special_args, vlen);

  overflow = __riscv_vmfgt (vx_cmp, RVV_EXP10F_CLEAN_MAX_INPUT_F, vlen);
  underflow = __riscv_vmflt (vx_cmp, RVV_EXP10F_CLEAN_MIN_INPUT_F, vlen);
  range_args = __riscv_vmor (overflow, underflow, vlen);

  vy_special = VFMV_VF (fp_posZero, vlen);

  UINT nb_special_args = __riscv_vcpop (special_args, vlen);
  if (nb_special_args > 0)
    {
      /* Substitute -Inf with +0 and propagate NaNs via x+x.  */
      VBOOL id_mask;
      IDENTIFY (vclass, class_negInf, id_mask, vlen);
      vx = __riscv_vfmerge (vx, fp_posZero, id_mask, vlen);
      vy_special = __riscv_vfadd (special_args, vx, vx, vlen);
      vx = __riscv_vfmerge (vx, fp_posZero, special_args, vlen);
    }

  vy_special = __riscv_vfmerge (vy_special, fp_posInf, overflow, vlen);

  VBOOL special_total = __riscv_vmor (special_args, range_args, vlen);
  vx = __riscv_vfmerge (vx, fp_posZero, range_args, vlen);

  /* Range reduction: 10^x = 2^q * 10^s, with s in [-log10(2)/2, log10(2)/2]. */
  VFLOAT n_flt = __riscv_vfmul (vx, RVV_EXP10F_CLEAN_LOG10_2_F, vlen);
  VINT rounded_n = __riscv_vfcvt_x (n_flt, vlen);
  VFLOAT u = __riscv_vfcvt_f (rounded_n, vlen);
  VINT q = __riscv_vfcvt_x (u, vlen);

  VFLOAT s = __riscv_vfmacc (vx, -RVV_EXP10F_CLEAN_L10U_F, u, vlen);
  s = __riscv_vfmacc (s, -RVV_EXP10F_CLEAN_L10L_F, u, vlen);

  /* Polynomial approximation for exp10(s): 1 + s*(ln10 + s*P(s)).  */
  VFLOAT poly = VFMV_VF (RVV_EXP10F_CLEAN_POLY0_F, vlen);
  poly = __riscv_vfmadd (poly, s,
                         VFMV_VF (RVV_EXP10F_CLEAN_POLY1_F, vlen), vlen);
  poly = __riscv_vfmadd (poly, s,
                         VFMV_VF (RVV_EXP10F_CLEAN_POLY2_F, vlen), vlen);
  poly = __riscv_vfmadd (poly, s,
                         VFMV_VF (RVV_EXP10F_CLEAN_POLY3_F, vlen), vlen);
  poly = __riscv_vfmadd (poly, s,
                         VFMV_VF (RVV_EXP10F_CLEAN_POLY4_F, vlen), vlen);
  poly = __riscv_vfmadd (poly, s,
                         VFMV_VF (RVV_EXP10F_CLEAN_POLY5_F, vlen), vlen);

  VFLOAT us = __riscv_vfmul (poly, s, vlen);

  VFLOAT x_hi, x_lo;
  DFADD_VF2_VF2_VF (VFMV_VF (RVV_EXP10F_CLEAN_LN10_HI_F, vlen),
                    VFMV_VF (RVV_EXP10F_CLEAN_LN10_LO_F, vlen), us,
                    x_hi, x_lo, vlen);

  VFLOAT t_hi, t_lo;
  DFMUL_VF2_VF2_VF (x_hi, x_lo, s, t_hi, t_lo, vlen);

  VFLOAT sum_hi, sum_lo;
  DFADD_VF2_VF_VF2 (VFMV_VF (1.0f, vlen), t_hi, t_lo, sum_hi, sum_lo, vlen);
  DFNORMALIZE_VF2_VF2 (sum_hi, sum_lo, sum_hi, sum_lo, vlen);

  vy = sum_hi;

  /* Scale by 2^q.  */
  vy = RVV_LIBMVEC_EXP10F_VLDEXP2 (vy, q, vlen);

  /* Incorporate special results.  */
  vy = __riscv_vmerge (vy, vy_special, special_total, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
