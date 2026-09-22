/* Single-precision vector exp10m1f function.
   exp10m1f(x) computes 10^x - 1

   Implemented via expm1f(x * ln(10)) using a SLEEF-derived expm1f kernel.

   Copyright (C) 2026 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <https://www.gnu.org/licenses/>.  */

#include <vecmath_config.h>
#include <v_math.h>

#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

/* Shared float-float helpers live in v_math.h.  */

/* dfsqu: (a_hi, a_lo)^2 -> (r_hi, r_lo).  */

/* expk2f: exp((d_hi, d_lo)) -> (r_hi, r_lo)
   SLEEF expk2f helper for extended precision exponential.  */
static inline void
RVV_LIBMVEC_IMPL_NAME (exp10m1f_expk2f_vf2, LMUL)
  (VFLOAT d_hi, VFLOAT d_lo, VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  /* Argument reduction: exp(d) = 2^q * exp(s) where |s| <= ln(2)/2.  */
  VFLOAT u = VFMUL_VF
    (VFADD_VV (d_hi, d_lo, vlen), RVV_EXP10M1F_SOURCE_R_LN2_F, vlen);
  VINT q = __riscv_vfcvt_x (u, vlen);
  VFLOAT q_flt = __riscv_vfcvt_f (q, vlen);

  /* s = d - q * ln(2) using double-float arithmetic.  */
  VFLOAT s_hi, s_lo;
  VFLOAT tmp = VFMUL_VF (q_flt, -RVV_EXP10M1F_SOURCE_L2U_F, vlen);
  DFADD2_VF2_VF (d_hi, d_lo, tmp, s_hi, s_lo, vlen);
  tmp = VFMUL_VF (q_flt, -RVV_EXP10M1F_SOURCE_L2L_F, vlen);
  DFADD2_VF2_VF (s_hi, s_lo, tmp, s_hi, s_lo, vlen);

  /* Polynomial approximation for exp(s) - 1.  */
  u = VFMV_VF (RVV_EXP10M1F_EXPM1_POLY0_F, vlen);
  u = VFMADD_VV
    (u, s_hi, VFMV_VF (RVV_EXP10M1F_EXPM1_POLY1_F, vlen), vlen);
  u = VFMADD_VV
    (u, s_hi, VFMV_VF (RVV_EXP10M1F_EXPM1_POLY2_F, vlen), vlen);
  u = VFMADD_VV
    (u, s_hi, VFMV_VF (RVV_EXP10M1F_EXPM1_POLY3_F, vlen), vlen);

  /* t = s * u + 0.166666659414234244790680580464.  */
  VFLOAT t_hi, t_lo;
  DFMUL_VF2_VF2_VF (s_hi, s_lo, u, t_hi, t_lo, vlen);
  DFADD2_VF2_VF_S (t_hi, t_lo, RVV_EXP10M1F_EXPM1_TAYLOR_C_F,
                   t_hi, t_lo, vlen);

  /* t = s * t + 0.5.  */
  VFLOAT tmp_hi, tmp_lo;
  DFMUL_VF2_VF2_VF2 (s_hi, s_lo, t_hi, t_lo, tmp_hi, tmp_lo, vlen);
  DFADD2_VF2_VF_S (tmp_hi, tmp_lo, 0.5f, t_hi, t_lo, vlen);

  /* t = s + s^2 * t.  */
  VFLOAT s2_hi, s2_lo;
  DFSQU_VF2_VF2 (s_hi, s_lo, s2_hi, s2_lo, vlen);
  DFMUL_VF2_VF2_VF2 (s2_hi, s2_lo, t_hi, t_lo, tmp_hi, tmp_lo, vlen);
  DFADD2_VF2_VF2_VF2 (s_hi, s_lo, tmp_hi, tmp_lo, t_hi, t_lo, vlen);

  /* t = 1 + t.  */
  DFADD2_VF2_VF_S (t_hi, t_lo, 1.0f, t_hi, t_lo, vlen);

  /* Scale by 2^q.  */
  FAST_LDEXP (t_hi, q, vlen);
  FAST_LDEXP (t_lo, q, vlen);

  /* Handle underflow: if d.x < -104, return 0.  */
  VBOOL underflow = __riscv_vmflt
    (d_hi, RVV_EXP10M1F_SOURCE_UNDERFLOW_F, vlen);
  t_hi = __riscv_vfmerge (t_hi, fp_posZero, underflow, vlen);
  t_lo = __riscv_vfmerge (t_lo, fp_posZero, underflow, vlen);

  *r_hi = t_hi;
  *r_lo = t_lo;
}
static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT vx = x;
  VFLOAT vy;
  VBOOL special_args, overflow_x, underflow_x;

  SET_ROUNDTONEAREST;

  VFLOAT x_orig = vx;
  VFLOAT abs_x = __riscv_vfabs (x_orig, vlen);

  /* Handle NaN/Inf and range shortcuts (align with v_f_logf_skeleton style).  */
  VUINT vclass = __riscv_vfclass (x_orig, vlen);
  IDENTIFY (vclass, class_NaN | class_Inf, special_args, vlen);
  overflow_x = __riscv_vmfgt (x_orig, RVV_EXP10M1F_CLEAN_MAX_INPUT_F, vlen);
  underflow_x = __riscv_vmflt (x_orig, RVV_EXP10M1F_CLEAN_MIN_INPUT_F, vlen);
  VBOOL special_total = VMOR (special_args, VMOR (overflow_x, underflow_x, vlen), vlen);

  UINT nb_special_total = __riscv_vcpop (special_total, vlen);
  if (nb_special_total > 0)
    {
      /* Remove special lanes from core computation.  */
      vx = VFMERGE_VFM (vx, fp_posZero_f, special_total, vlen);
    }

  /* Scale x -> t = x * ln(10) in double-float to preserve accuracy when
     |x| is large (avoids large ULP error from rounding t to float early).  */
  VFLOAT ln10_hi = VFMV_VF (RVV_EXP10M1F_CLEAN_LN10_HI_F, vlen);
  VFLOAT ln10_lo = VFMV_VF (RVV_EXP10M1F_CLEAN_LN10_LO_F, vlen);
  VFLOAT t_hi, t_lo;
  DFMUL_VF2_VF2_VF (ln10_hi, ln10_lo, vx, t_hi, t_lo, vlen);

  /* Float-rounded t for classification and tiny/overflow/underflow checks.  */
  VFLOAT t_orig = VFADD_VV (t_hi, t_lo, vlen);

  /* expm1(t) = expk2f(t) - 1 using double-float arithmetic.  */
  VFLOAT exp_hi, exp_lo;
  RVV_LIBMVEC_IMPL_NAME (exp10m1f_expk2f_vf2, LMUL)
    (t_hi, t_lo, &exp_hi, &exp_lo, vlen);
  VFLOAT y_hi, y_lo;
  DFADD2_VF2_VF_S (exp_hi, exp_lo, -1.0f, y_hi, y_lo, vlen);
  vy = __riscv_vfadd (y_hi, y_lo, vlen);

  /* Use t for tiny inputs to avoid cancellation.  */
  VFLOAT abs_t = __riscv_vfabs (t_orig, vlen);
  VBOOL tiny = __riscv_vmflt
    (abs_t, RVV_EXP10M1F_SOURCE_EXPM1_TINY_F, vlen);
  vy = __riscv_vmerge (vy, t_orig, tiny, vlen);

  /* For subnormals and signed zero, exp10m1f(x) rounds to x * ln(10).  */
  VBOOL tiny_input = __riscv_vmflt
    (abs_x, RVV_EXP10M1F_CLEAN_MIN_NORMAL_F, vlen);
  if (__riscv_vcpop (tiny_input, vlen) > 0)
    {
      VFLOAT tiny_y = __riscv_vfmul
        (vx, RVV_EXP10M1F_CLEAN_LN10_HI_F, vlen);
      tiny_y = __riscv_vfmadd
        (vx, RVV_EXP10M1F_CLEAN_LN10_LO_F, tiny_y, vlen);
      vy = __riscv_vmerge (vy, tiny_y, tiny_input, vlen);
    }

  /* Handle overflow and -0.0; the negative tail is selected in x-domain.  */
  VBOOL overflow = __riscv_vmfgt
    (t_orig, RVV_EXP10M1F_SOURCE_EXPM1_MAX_F, vlen);
  vy = __riscv_vfmerge (vy, fp_posInf_f, overflow, vlen);
  VBOOL negzero = __riscv_vmseq (F_AS_U_F (x_orig), VMVU_VX_F (0x80000000, vlen),
                                 vlen);
  vy = __riscv_vfmerge (vy, fp_negZero_f, negzero, vlen);

  /* Incorporate special results after the finite core.  Keep these merges in
     the destination vector instead of carrying a separate special-result
     vector across the core.  Besides shortening the live vector state, this
     avoids a GCC RVV LMUL spill/reload failure where the long-lived vector
     was reloaded as zero for m2/m4/m8.  The EH classification and sanitized
     core path remain unchanged.  */
  if (nb_special_total > 0)
    {
      VBOOL id_mask;

      /* +Inf or overflow -> +Inf.  */
      IDENTIFY (vclass, class_posInf, id_mask, vlen);
      vy = VFMERGE_VFM (vy, fp_posInf_f, id_mask, vlen);
      vy = VFMERGE_VFM (vy, fp_posInf_f, overflow_x, vlen);

      /* -Inf or underflow -> -1.  */
      IDENTIFY (vclass, class_negInf, id_mask, vlen);
      vy = VFMERGE_VFM (vy, fp_negOne_f, id_mask, vlen);
      vy = VFMERGE_VFM (vy, fp_negOne_f, underflow_x, vlen);

      /* NaNs: propagate via x + x.  */
      IDENTIFY (vclass, class_NaN, id_mask, vlen);
      vy = VFADD_VV_M (id_mask, vy, x_orig, vlen);
      vy = VFADD_VV_M (id_mask, vy, x_orig, vlen);
    }

  /* Once 10^x drops below half an ULP at 1.0f, exp10m1f(x) rounds to -1.  */
  VBOOL round_to_neg_one = __riscv_vmflt
    (x_orig, RVV_EXP10M1F_CLEAN_MIN_INPUT_F, vlen);
  VBOOL finite_x = __riscv_vmnot (special_args, vlen);
  round_to_neg_one = __riscv_vmand (round_to_neg_one, finite_x, vlen);
  vy = __riscv_vfmerge (vy, fp_negOne_f, round_to_neg_one, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
