/* Single-precision vector cbrt function.
   cbrtf(x) = x^(1/3)

   Based on SLEEF xcbrtf_u1 implementation.

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

/* Shared float-float helpers live in v_math.h.  */

static inline VINT
RVV_LIBMVEC_IMPL_NAME (vilogbk, LMUL)
  (VFLOAT d, size_t vlen)
{
  VBOOL o = __riscv_vmflt
    (d, VFMV_VF (RVV_SLEEF_CBRTF_VILOGBK_THRESH, vlen), vlen);
  VFLOAT d_scaled = __riscv_vfmul
    (VFMV_VF (RVV_SLEEF_CBRTF_VILOGBK_SCALE, vlen), d, vlen);
  d = __riscv_vmerge (d, d_scaled, o, vlen);

  VUINT q = __riscv_vsrl (F_AS_U (d), MAN_LEN_F, vlen);
  q = __riscv_vand (q, RVV_SLEEF_CBRTF_VILOGBK_MASK, vlen);

  VINT bias = __riscv_vmerge
    (VMVI_VX (RVV_SLEEF_CBRTF_VILOGBK_SMALL_BIAS, vlen),
     VMVI_VX (RVV_SLEEF_CBRTF_VILOGBK_LARGE_BIAS, vlen), o, vlen);
  return __riscv_vsub (U_AS_I (q), bias, vlen);
}

static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (vpow2i, LMUL) (VINT q, size_t vlen)
{
  VINT q_bias = __riscv_vadd (q, VMVI_VX (EXP_BIAS_F, vlen), vlen);
  q_bias = __riscv_vsll (q_bias, MAN_LEN_F, vlen);
  return I_AS_F (q_bias);
}

static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (vldexp2, LMUL)
  (VFLOAT d, VINT e, size_t vlen)
{
  VINT e1 = __riscv_vsra (e, 1, vlen);
  VINT e2 = __riscv_vsub (e, e1, vlen);
  VFLOAT s1 = RVV_LIBMVEC_IMPL_NAME (vpow2i, LMUL) (e1, vlen);
  VFLOAT s2 = RVV_LIBMVEC_IMPL_NAME (vpow2i, LMUL) (e2, vlen);
  return __riscv_vfmul (__riscv_vfmul (d, s1, vlen), s2, vlen);
}

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT d, abs_d, z, z_special;
  VBOOL special_args;
  VINT e, qu, re;
  VFLOAT q2_hi, q2_lo;
  SET_ROUNDTONEAREST;

  d = x;
  z_special = x;
  VUINT vclass = __riscv_vfclass (d, vlen);
  IDENTIFY (vclass, class_NaN | class_Inf | class_Zero, special_args, vlen);
  if (__riscv_vcpop (special_args, vlen) > 0)
    {
      z_special = __riscv_vfadd (special_args, d, d, vlen);
      d = __riscv_vfmerge (d, fp_posOne, special_args, vlen);
    }
  if (__riscv_vcpop (special_args, vlen) == vlen)
    {
      RESTORE_FRM;
      return z_special;
    }
  abs_d = __riscv_vfabs (d, vlen);

  e = __riscv_vadd (
    RVV_LIBMVEC_IMPL_NAME (vilogbk, LMUL)
      (abs_d, vlen), RVV_SLEEF_CBRTF_EXPONENT_OFFSET, vlen);
  d = RVV_LIBMVEC_IMPL_NAME (vldexp2, LMUL)
    (d, __riscv_vneg (e, vlen), vlen);

  /* Compute qu = trunc ((e + 6144) / 3) and re = (e + 6144) mod 3 with an
     exact bounded integer division.  */
  VINT t_i = __riscv_vadd (e, RVV_SLEEF_CBRTF_QUOTIENT_OFFSET, vlen);
  qu = __riscv_vsra
    (__riscv_vmul (t_i, RVV_SLEEF_CBRTF_QUOTIENT_MUL, vlen),
     RVV_SLEEF_CBRTF_QUOTIENT_SHIFT, vlen);
  re = __riscv_vsub
    (t_i, __riscv_vmul (qu, RVV_SLEEF_CBRTF_QUOTIENT_DIVISOR, vlen), vlen);

  q2_hi = VFMV_VF (fp_posOne_f, vlen);
  q2_lo = VFMV_VF (fp_posZero_f, vlen);
  VBOOL re_eq1 = __riscv_vmseq (re, RVV_SLEEF_CBRTF_REMAINDER_ONE, vlen);
  VBOOL re_eq2 = __riscv_vmseq (re, RVV_SLEEF_CBRTF_REMAINDER_TWO, vlen);
  q2_hi = __riscv_vmerge
    (q2_hi, VFMV_VF (RVV_SLEEF_CBRTF_Q2_1_HI, vlen), re_eq1, vlen);
  q2_lo = __riscv_vmerge
    (q2_lo, VFMV_VF (RVV_SLEEF_CBRTF_Q2_1_LO, vlen), re_eq1, vlen);
  q2_hi = __riscv_vmerge
    (q2_hi, VFMV_VF (RVV_SLEEF_CBRTF_Q2_2_HI, vlen), re_eq2, vlen);
  q2_lo = __riscv_vmerge
    (q2_lo, VFMV_VF (RVV_SLEEF_CBRTF_Q2_2_LO, vlen), re_eq2, vlen);

  q2_hi = __riscv_vfsgnjx (q2_hi, d, vlen);
  q2_lo = __riscv_vfsgnjx (q2_lo, d, vlen);
  d = __riscv_vfabs (d, vlen);

  VFLOAT xpoly = VFMV_VF (RVV_SLEEF_CBRTF_POLY0_F, vlen);
  xpoly = __riscv_vfmadd (xpoly, d,
                          VFMV_VF (RVV_SLEEF_CBRTF_POLY1_F, vlen), vlen);
  xpoly = __riscv_vfmadd (xpoly, d,
                          VFMV_VF (RVV_SLEEF_CBRTF_POLY2_F, vlen), vlen);
  xpoly = __riscv_vfmadd (xpoly, d,
                          VFMV_VF (RVV_SLEEF_CBRTF_POLY3_F, vlen), vlen);
  xpoly = __riscv_vfmadd (xpoly, d,
                          VFMV_VF (RVV_SLEEF_CBRTF_POLY4_F, vlen), vlen);
  xpoly = __riscv_vfmadd (xpoly, d,
                          VFMV_VF (RVV_SLEEF_CBRTF_POLY5_F, vlen), vlen);

  VFLOAT y = __riscv_vfmul (xpoly, xpoly, vlen);
  y = __riscv_vfmul (y, y, vlen);
  /* Keep the source's rounded residual multiply and subtract separate.  The
     fused x + residual/3 form changes the finite Newton rounding and is not
     an EH adaptation.  */
  VFLOAT residual = __riscv_vfnmsub (d, y, xpoly, vlen);
  VFLOAT correction = __riscv_vfmul
    (residual, RVV_SLEEF_CBRTF_NEG_ONE_THIRD, vlen);
  xpoly = __riscv_vfsub (xpoly, correction, vlen);

  z = xpoly;

  VFLOAT u_hi, u_lo;
  RVV_SLEEF_ORIG_DFMUL_VF_VF
    (xpoly, xpoly, u_hi, u_lo, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (u_hi, u_lo, u_hi, u_lo, u_hi, u_lo, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (u_hi, u_lo, d, u_hi, u_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF
    (u_hi, u_lo, __riscv_vfneg (xpoly, vlen), u_hi, u_lo, vlen);
  y = __riscv_vfadd (u_hi, u_lo, vlen);

  y = __riscv_vfmul
    (VFMV_VF (RVV_SLEEF_CBRTF_NEG_TWO_THIRDS, vlen), y, vlen);
  y = __riscv_vfmul (y, z, vlen);

  VFLOAT v_hi, v_lo;
  RVV_SLEEF_ORIG_DFMUL_VF_VF
    (z, z, v_hi, v_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF
    (v_hi, v_lo, y, v_hi, v_lo, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (v_hi, v_lo, d, v_hi, v_lo, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (v_hi, v_lo, q2_hi, q2_lo, v_hi, v_lo, vlen);

  VFLOAT z_val = __riscv_vfadd (v_hi, v_lo, vlen);
  VINT exp_adj = __riscv_vsub
    (qu, VMVI_VX (RVV_SLEEF_CBRTF_EXP_ADJUST, vlen), vlen);
  z = RVV_LIBMVEC_IMPL_NAME (vldexp2, LMUL) (z_val, exp_adj, vlen);

  z = __riscv_vmerge (z, z_special, special_args, vlen);

  RESTORE_FRM;
  return z;
}

#undef LMUL
#undef SIMDLEN
