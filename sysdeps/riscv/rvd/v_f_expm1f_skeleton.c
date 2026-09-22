/* Single-precision vector expm1 function.
   expm1f(x) computes e^x - 1

   Based on SLEEF xexpm1f (u10 precision).

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

/* ============================================================================
   Double-float helper macros for RVV (vfloat2 operations)
   ============================================================================ */

/* dfadd2: (a_hi, a_lo) + scalar b -> (r_hi, r_lo) */

/* Shared float-float helpers live in v_math.h.  */

/* dfsqu: (a_hi, a_lo)^2 -> (r_hi, r_lo) */

/* Keep the source SLEEF contracts that are not represented by the generic
   clean helpers local to this routine.  In particular, the source pair
   multiply has a deliberate nested-FMA association, and vldexp2 performs two
   floating-point multiplies rather than exponent-bit injection.  */
static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (expm1f_source_vpow2i, LMUL) (VINT q, size_t vlen)
{
  VINT biased = __riscv_vadd (q, EXP_BIAS_F, vlen);
  return I_AS_F (__riscv_vsll (biased, MAN_LEN_F, vlen));
}

static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (expm1f_source_vldexp2, LMUL)
  (VFLOAT d, VINT e, size_t vlen)
{
  VINT e_half = __riscv_vsra (e, 1, vlen);
  VINT e_rest = __riscv_vsub (e, e_half, vlen);
  VFLOAT p_half = RVV_LIBMVEC_IMPL_NAME (expm1f_source_vpow2i, LMUL)
    (e_half, vlen);
  VFLOAT p_rest = RVV_LIBMVEC_IMPL_NAME (expm1f_source_vpow2i, LMUL)
    (e_rest, vlen);
  return VFMUL_VV (VFMUL_VV (d, p_half, vlen), p_rest, vlen);
}

static inline void
RVV_LIBMVEC_IMPL_NAME (expm1f_source_dfmul_pair_pair, LMUL)
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
   SLEEF expk2f helper for extended precision exponential */
static inline void
RVV_LIBMVEC_IMPL_NAME (expm1f_expk2f_vf2, LMUL)
  (VFLOAT d_hi, VFLOAT d_lo, VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  /* Argument reduction: exp(d) = 2^q * exp(s) where |s| <= ln(2)/2 */
  VFLOAT u = VFMUL_VF (VFADD_VV (d_hi, d_lo, vlen),
                       RVV_SLEEF_EXPM1F_R_LN2,
                       vlen);
  VINT q = __riscv_vfcvt_x (u, vlen);
  VFLOAT q_flt = __riscv_vfcvt_f (q, vlen);

  /* s = d - q * ln(2) using double-float arithmetic */
  VFLOAT s_hi, s_lo;
  VFLOAT tmp = VFMUL_VF (q_flt, -RVV_SLEEF_EXPM1F_L2U, vlen);
  DFADD2_VF2_VF (d_hi, d_lo, tmp, s_hi, s_lo, vlen);
  tmp = VFMUL_VF (q_flt, -RVV_SLEEF_EXPM1F_L2L, vlen);
  DFADD2_VF2_VF (s_hi, s_lo, tmp, s_hi, s_lo, vlen);

  /* Polynomial approximation for exp(s) - 1 */
  u = VFMV_VF (RVV_SLEEF_EXPM1F_POLY0, vlen);
  u = VFMADD_VV (u, s_hi, VFMV_VF (RVV_SLEEF_EXPM1F_POLY1, vlen), vlen);
  u = VFMADD_VV (u, s_hi, VFMV_VF (RVV_SLEEF_EXPM1F_POLY2, vlen), vlen);
  u = VFMADD_VV (u, s_hi, VFMV_VF (RVV_SLEEF_EXPM1F_POLY3, vlen), vlen);

  /* t = s * u + 0.166666659414234244790680580464 */
  VFLOAT t_hi, t_lo;
  DFMUL_VF2_VF2_VF (s_hi, s_lo, u, t_hi, t_lo, vlen);
  DFADD2_VF2_VF_S (t_hi, t_lo, RVV_SLEEF_EXPM1F_TAYLOR_C,
                   t_hi, t_lo, vlen);

  /* t = s * t + 0.5 */
  VFLOAT tmp_hi, tmp_lo;
  RVV_LIBMVEC_IMPL_NAME (expm1f_source_dfmul_pair_pair, LMUL)
    (s_hi, s_lo, t_hi, t_lo, &tmp_hi, &tmp_lo, vlen);
  DFADD2_VF2_VF_S (tmp_hi, tmp_lo, 0.5f, t_hi, t_lo, vlen);

  /* t = s + s^2 * t */
  VFLOAT s2_hi, s2_lo;
  DFSQU_VF2_VF2 (s_hi, s_lo, s2_hi, s2_lo, vlen);
  RVV_LIBMVEC_IMPL_NAME (expm1f_source_dfmul_pair_pair, LMUL)
    (s2_hi, s2_lo, t_hi, t_lo, &tmp_hi, &tmp_lo, vlen);
  DFADD2_VF2_VF2_VF2 (s_hi, s_lo, tmp_hi, tmp_lo, t_hi, t_lo, vlen);

  /* t = 1 + t */
  DFADD_VF2_VF_VF2 (VFMV_VF (1.0f, vlen), t_hi, t_lo,
                    t_hi, t_lo, vlen);

  /* Scale by 2^q */
  t_hi = RVV_LIBMVEC_IMPL_NAME (expm1f_source_vldexp2, LMUL)
    (t_hi, q, vlen);
  t_lo = RVV_LIBMVEC_IMPL_NAME (expm1f_source_vldexp2, LMUL)
    (t_lo, q, vlen);

  /* Handle underflow: if d.x < -104, return 0 */
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

  VFLOAT vx_orig = vx;
  vy_special = VFMV_VF (0.0f, vlen);
  VUINT vclass = __riscv_vfclass (vx, vlen);
  IDENTIFY (vclass, class_NaN | class_Inf, special_args, vlen);
  VFLOAT vx_cmp = __riscv_vfmerge (vx_orig, fp_posZero, special_args, vlen);
  overflow_x = __riscv_vmfgt (vx_cmp, RVV_SLEEF_EXPM1F_MAX, vlen);
  underflow_x = __riscv_vmflt (vx_cmp, RVV_SLEEF_EXPM1F_MIN, vlen);
  VBOOL special_total
      = __riscv_vmor (special_args, __riscv_vmor (overflow_x, underflow_x, vlen),
                      vlen);

  UINT nb_special_total = __riscv_vcpop (special_total, vlen);
  if (nb_special_total > 0)
    {
      /* Substitute -Inf with -1 */
      VBOOL id_mask;
      IDENTIFY (vclass, class_negInf, id_mask, vlen);
      vx = __riscv_vfmerge (vx, fp_negOne, id_mask, vlen);
      vy_special = __riscv_vfmul (special_args, vx, fp_posOne, vlen);

      /* Finite positive overflow maps to +Inf, finite large negative inputs
         that round to -1 stay on the special path.  */
      vy_special = __riscv_vfmerge (vy_special, fp_posInf, overflow_x, vlen);
      vy_special = __riscv_vfmerge (vy_special, fp_negOne, underflow_x, vlen);

      vx = __riscv_vfmerge (vx, fp_posZero, special_total, vlen);
    }

  /* expm1(x) = expk2f(x) - 1 using double-float arithmetic */
  VFLOAT exp_hi, exp_lo;
  RVV_LIBMVEC_IMPL_NAME (expm1f_expk2f_vf2, LMUL)
    (vx, VFMV_VF (0.0f, vlen), &exp_hi, &exp_lo, vlen);
  VFLOAT y_hi, y_lo;
  DFADD2_VF2_VF_S (exp_hi, exp_lo, -1.0f, y_hi, y_lo, vlen);
  vy = __riscv_vfadd (y_hi, y_lo, vlen);

  /* Handle -0.0. */
  VBOOL negzero = __riscv_vmseq (F_AS_U (vx_orig), VMVU_VX (0x80000000, vlen), vlen);
  vy = __riscv_vfmerge (vy, fp_negZero, negzero, vlen);

  /* Keep the finite range exceptions visible even though their unsafe input
     lanes were sanitized before the finite core.  The probe values are
     discarded by the final special merge below; only their fflags effects
     are architecturally relevant.  */
  if (nb_special_total > 0)
    {
      VFLOAT overflow_probe = VFMV_VF (RVV_SLEEF_EXPM1F_MAX_FINITE, vlen);
      overflow_probe
          = __riscv_vfmul (overflow_x, overflow_probe, 2.0f, vlen);
      vy = VMERGE_VVM_F (vy, overflow_probe, overflow_x, vlen);

      VFLOAT underflow_probe = VFMV_VF (fp_negOne, vlen);
      underflow_probe
          = __riscv_vfadd (underflow_x, underflow_probe,
                           RVV_SLEEF_EXPM1F_MIN_SUBNORMAL, vlen);
      vy = VMERGE_VVM_F (vy, underflow_probe, underflow_x, vlen);
    }

  vy = __riscv_vmerge (vy, vy_special, special_total, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
