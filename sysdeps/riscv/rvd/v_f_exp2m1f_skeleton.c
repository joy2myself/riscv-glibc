/* Single-precision vector exp2m1f function.
   exp2m1f(x) computes 2^x - 1

   Implemented via expm1(x * ln(2)) using a SLEEF-derived expm1f kernel.

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

/* dfsqu: (a_hi, a_lo)^2 -> (r_hi, r_lo).  */

static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (exp2m1f_source_vpow2i, LMUL) (VINT q, size_t vlen)
{
  VINT biased = __riscv_vadd (q, EXP_BIAS_F, vlen);
  return I_AS_F (__riscv_vsll (biased, MAN_LEN_F, vlen));
}

static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (exp2m1f_source_vldexp2, LMUL)
  (VFLOAT d, VINT e, size_t vlen)
{
  VINT e_half = __riscv_vsra (e, 1, vlen);
  VINT e_rest = __riscv_vsub (e, e_half, vlen);
  VFLOAT p_half = RVV_LIBMVEC_IMPL_NAME (exp2m1f_source_vpow2i, LMUL)
    (e_half, vlen);
  VFLOAT p_rest = RVV_LIBMVEC_IMPL_NAME (exp2m1f_source_vpow2i, LMUL)
    (e_rest, vlen);
  return VFMUL_VV (VFMUL_VV (d, p_half, vlen), p_rest, vlen);
}

static inline void
RVV_LIBMVEC_IMPL_NAME (exp2m1f_source_dfmul_pair_pair, LMUL)
  (VFLOAT x_hi, VFLOAT x_lo, VFLOAT y_hi, VFLOAT y_lo,
   VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT s = VFMUL_VV (x_hi, y_hi, vlen);
  VFLOAT residual = VFMSUB_VV (x_hi, y_hi, s, vlen);
  VFLOAT cross_hi = VFMADD_VV (x_lo, y_hi, residual, vlen);
  *r_hi = s;
  *r_lo = VFMADD_VV (x_hi, y_lo, cross_hi, vlen);
}

/* expk2f: exp((d_hi, d_lo)) -> (r_hi, r_lo)
   SLEEF expk2f helper for extended precision exponential.  */
static inline void
RVV_LIBMVEC_IMPL_NAME (exp2m1f_expk2f_vf2, LMUL)
  (VFLOAT d_hi, VFLOAT d_lo, VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  /* Argument reduction: exp(d) = 2^q * exp(s) where |s| <= ln(2)/2.  */
  VFLOAT u = VFMUL_VF (VFADD_VV (d_hi, d_lo, vlen),
                       RVV_SLEEF_EXPM1F_R_LN2,
                       vlen);
  VINT q = __riscv_vfcvt_x (u, vlen);
  VFLOAT q_flt = __riscv_vfcvt_f (q, vlen);

  /* s = d - q * ln(2) using double-float arithmetic.  */
  VFLOAT s_hi, s_lo;
  VFLOAT tmp = VFMUL_VF (q_flt, -RVV_SLEEF_EXPM1F_L2U, vlen);
  DFADD2_VF2_VF (d_hi, d_lo, tmp, s_hi, s_lo, vlen);
  tmp = VFMUL_VF (q_flt, -RVV_SLEEF_EXPM1F_L2L, vlen);
  DFADD2_VF2_VF (s_hi, s_lo, tmp, s_hi, s_lo, vlen);

  /* Polynomial approximation for exp(s) - 1.  */
  u = VFMV_VF (RVV_SLEEF_EXPM1F_POLY0, vlen);
  u = VFMADD_VV (u, s_hi, VFMV_VF (RVV_SLEEF_EXPM1F_POLY1, vlen), vlen);
  u = VFMADD_VV (u, s_hi, VFMV_VF (RVV_SLEEF_EXPM1F_POLY2, vlen), vlen);
  u = VFMADD_VV (u, s_hi, VFMV_VF (RVV_SLEEF_EXPM1F_POLY3, vlen), vlen);

  /* t = s * u + 0.166666659414234244790680580464.  */
  VFLOAT t_hi, t_lo;
  DFMUL_VF2_VF2_VF (s_hi, s_lo, u, t_hi, t_lo, vlen);
  DFADD2_VF2_VF_S (t_hi, t_lo, RVV_SLEEF_EXPM1F_TAYLOR_C,
                   t_hi, t_lo, vlen);

  /* t = s * t + 0.5.  */
  VFLOAT tmp_hi, tmp_lo;
  RVV_LIBMVEC_IMPL_NAME (exp2m1f_source_dfmul_pair_pair, LMUL)
    (s_hi, s_lo, t_hi, t_lo, &tmp_hi, &tmp_lo, vlen);
  DFADD2_VF2_VF_S (tmp_hi, tmp_lo, 0.5f, t_hi, t_lo, vlen);

  /* t = s + s^2 * t.  */
  VFLOAT s2_hi, s2_lo;
  DFSQU_VF2_VF2 (s_hi, s_lo, s2_hi, s2_lo, vlen);
  RVV_LIBMVEC_IMPL_NAME (exp2m1f_source_dfmul_pair_pair, LMUL)
    (s2_hi, s2_lo, t_hi, t_lo, &tmp_hi, &tmp_lo, vlen);
  DFADD2_VF2_VF2_VF2 (s_hi, s_lo, tmp_hi, tmp_lo, t_hi, t_lo, vlen);

  /* t = 1 + t.  */
  DFADD_VF2_VF_VF2 (VFMV_VF (1.0f, vlen), t_hi, t_lo,
                    t_hi, t_lo, vlen);

  /* Scale by 2^q.  */
  t_hi = RVV_LIBMVEC_IMPL_NAME (exp2m1f_source_vldexp2, LMUL)
    (t_hi, q, vlen);
  t_lo = RVV_LIBMVEC_IMPL_NAME (exp2m1f_source_vldexp2, LMUL)
    (t_lo, q, vlen);

  /* Handle underflow: if d.x < -104, return 0.  */
  VBOOL underflow = __riscv_vmflt (d_hi, RVV_SLEEF_EXPM1F_UNDERFLOW, vlen);
  t_hi = __riscv_vfmerge (t_hi, fp_posZero, underflow, vlen);
  t_lo = __riscv_vfmerge (t_lo, fp_posZero, underflow, vlen);

  *r_hi = t_hi;
  *r_lo = t_lo;
}
static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT vx = x;
  VFLOAT vy, vy_special;
  VBOOL special_args, overflow_x, underflow_x;

  SET_ROUNDTONEAREST;

  VFLOAT x_orig = vx;

  /* Handle NaN/Inf and out-of-range inputs (align with v_f_logf_skeleton).  */
  VUINT vclass = __riscv_vfclass (x_orig, vlen);
  IDENTIFY (vclass, class_NaN | class_Inf, special_args, vlen);
  VFLOAT x_cmp = __riscv_vfmerge (x_orig, fp_posZero_f, special_args, vlen);
  overflow_x = __riscv_vmfge (x_cmp, RVV_LIBMVEC_EXP2M1F_MAX, vlen);
  underflow_x = __riscv_vmfle (x_cmp, RVV_LIBMVEC_EXP2M1F_MIN, vlen);
  VBOOL special_total
      = VMOR (special_args, VMOR (overflow_x, underflow_x, vlen), vlen);

  /* Identify the small exact-integer lattice without ever converting an
     out-of-range or special lane.  The enclosing SET_ROUNDTONEAREST makes
     the round-trip conversion a safe integer test; fractional lanes may
     raise the inexact bit here, which is already part of their scalar
     contract, while exact-integer lanes do not.  */
  VFLOAT abs_x_cmp = __riscv_vfabs (x_cmp, vlen);
  VBOOL small_x = __riscv_vmfle
    (abs_x_cmp, RVV_LIBMVEC_EXP2M1F_EXACT_INT_LIMIT, vlen);
  VFLOAT x_for_int = __riscv_vfmerge
    (x_cmp, fp_posZero_f, __riscv_vmnot (small_x, vlen), vlen);
  VINT x_int = __riscv_vfcvt_x (x_for_int, vlen);
  VFLOAT x_back = __riscv_vfcvt_f (x_int, vlen);
  VBOOL exact_integer = __riscv_vmfeq (x_for_int, x_back, vlen);
  exact_integer = __riscv_vmand
    (exact_integer, small_x, vlen);
  exact_integer = __riscv_vmand
    (exact_integer, __riscv_vmnot (special_args, vlen), vlen);

  vy_special = VFMV_VF_F (0.0f, vlen);

  UINT nb_special_total = __riscv_vcpop (special_total, vlen);
  if (nb_special_total > 0)
    {
      VBOOL id_mask;

      /* +Inf or overflow -> +Inf.  */
      IDENTIFY (vclass, class_posInf, id_mask, vlen);
      vy_special = VFMERGE_VFM (vy_special, fp_posInf_f, id_mask, vlen);
      vy_special = VFMERGE_VFM (vy_special, fp_posInf_f, overflow_x, vlen);

      /* -Inf or underflow -> -1.  */
      IDENTIFY (vclass, class_negInf, id_mask, vlen);
      vy_special = VFMERGE_VFM (vy_special, fp_negOne_f, id_mask, vlen);
      vy_special = VFMERGE_VFM (vy_special, fp_negOne_f, underflow_x, vlen);

      /* NaNs: propagate via x + x.  */
      IDENTIFY (vclass, class_NaN, id_mask, vlen);
      vy_special = VFADD_VV_M (id_mask, vy_special, x_orig, vlen);
      vy_special = VFADD_VV_M (id_mask, vy_special, x_orig, vlen);

      /* Remove special lanes from core computation.  */
      vx = VFMERGE_VFM (vx, fp_posZero_f, special_total, vlen);
    }

  /* Exact integer lanes are handled by a quiet finite provider below.  Zero
     them before the transformed expm1 core so its internal reduction cannot
     leak inexact for an otherwise exact result.  */
  vx = VFMERGE_VFM (vx, fp_posZero_f, exact_integer, vlen);

  /* Scale x -> t = x * ln(2) in double-float.  */
  VFLOAT ln2_hi = VFMV_VF (RVV_SLEEF_EXPM1F_L2U, vlen);
  VFLOAT ln2_lo = VFMV_VF (RVV_SLEEF_EXPM1F_L2L, vlen);
  VFLOAT t_hi, t_lo;
  DFMUL_VF2_VF2_VF (ln2_hi, ln2_lo, vx, t_hi, t_lo, vlen);
  VFLOAT t_orig = VFADD_VV (t_hi, t_lo, vlen);

  /* exp2m1(x) = expm1(t).  */
  VFLOAT exp_hi, exp_lo;
  RVV_LIBMVEC_IMPL_NAME (exp2m1f_expk2f_vf2, LMUL)
    (t_hi, t_lo, &exp_hi, &exp_lo, vlen);
  VFLOAT y_hi, y_lo;
  DFADD2_VF2_VF_S (exp_hi, exp_lo, -1.0f, y_hi, y_lo, vlen);
  vy = __riscv_vfadd (y_hi, y_lo, vlen);

  /* Handle overflow and -0.0; the negative tail is selected in x-domain.  */
  VBOOL overflow = __riscv_vmfgt (t_orig, RVV_SLEEF_EXPM1F_MAX, vlen);
  vy = __riscv_vfmerge (vy, fp_posInf_f, overflow, vlen);

  /* For -24 <= n <= 24, 2^n - 1 is exactly representable.  Constructing
     the power of two from its exponent field and subtracting one is exact;
     unlike the transformed core, it does not expose a bookkeeping inexact
     flag.  The raw negative-zero merge below remains authoritative for -0.  */
  VINT exact_exp = __riscv_vadd (x_int, EXP_BIAS_F, vlen);
  VFLOAT exact_pow = I_AS_F
    (__riscv_vsll (exact_exp, MAN_LEN_F, vlen));
  VFLOAT exact_result = VFSUB_VF (exact_pow, 1.0f, vlen);
  vy = VMERGE_VVM_F (vy, exact_result, exact_integer, vlen);

  VBOOL negzero = __riscv_vmseq (F_AS_U_F (x_orig),
                                  VMVU_VX_F (0x80000000, vlen), vlen);
  vy = __riscv_vfmerge (vy, fp_negZero_f, negzero, vlen);

  /* Keep the finite range exceptions visible even though their unsafe input
     lanes were sanitized before the finite core.  The probe results are
     discarded by the final special merge below; only their fflags effects
     are architecturally relevant.  */
  if (nb_special_total > 0)
    {
      VFLOAT overflow_probe = VFMV_VF (RVV_SLEEF_EXPM1F_MAX_FINITE, vlen);
      overflow_probe
          = __riscv_vfmul (overflow_x, overflow_probe, 2.0f, vlen);
      vy = VMERGE_VVM_F (vy, overflow_probe, overflow_x, vlen);

      VFLOAT underflow_probe = VFMV_VF (fp_negOne_f, vlen);
      underflow_probe
          = __riscv_vfadd (underflow_x, underflow_probe,
                           RVV_SLEEF_EXPM1F_MIN_SUBNORMAL, vlen);
      vy = VMERGE_VVM_F (vy, underflow_probe, underflow_x, vlen);
    }

  /* Incorporate special results.  */
  vy = VMERGE_VVM_F (vy, vy_special, special_total, vlen);

  /* Once 2^x drops below half an ULP at 1.0f, exp2m1f(x) rounds to -1.  */
  VBOOL round_to_neg_one = underflow_x;
  vy = __riscv_vfmerge (vy, fp_negOne_f, round_to_neg_one, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
