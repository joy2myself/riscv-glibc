/* Single-precision vector erfc function.
   erfcf(x) = 1 - erf(x) computes the complementary error function.

   Based on SLEEF xerfcf_u15 (u15 precision).

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

/* Shared float-float helpers live in v_math.h.  */

/* dfsqu: (a_hi, a_lo)^2 -> (r_hi, r_lo).  */

/* Helpers must be uniquely named since this skeleton is included multiple
   times with different LMUL/SIMDLEN in v_f_importer.c.  */
/* Keep the source SLEEF finite graph local to the clean implementation.  The
   public EH envelope below remains clean-specific, but these helpers retain
   the source operation association where the generic helpers in v_math.h
   intentionally use a different graph.  */
static inline void
RVV_ERFCF_HELPER (source_dfmul_pair_pair)
  (VFLOAT a_hi, VFLOAT a_lo, VFLOAT b_hi, VFLOAT b_lo,
   VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT s = VFMUL_VV (a_hi, b_hi, vlen);
  VFLOAT residual = VFMSUB_VV (a_hi, b_hi, s, vlen);
  VFLOAT cross = VFMADD_VV (a_lo, b_hi, residual, vlen);
  *r_hi = s;
  *r_lo = VFMADD_VV (a_hi, b_lo, cross, vlen);
}

static inline void
RVV_ERFCF_HELPER (source_dfdiv_pair_pair)
  (VFLOAT n_hi, VFLOAT n_lo, VFLOAT d_hi, VFLOAT d_lo,
   VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT one = VFMV_VF_F (1.0f, vlen);
  VFLOAT t = VFDIV_VV (one, d_hi, vlen);
  VFLOAT s = VFMUL_VV (n_hi, t, vlen);
  VFLOAT u = VFMSUB_VV (t, n_hi, s, vlen);
  VFLOAT v = VFNMSUB_VV (d_hi, t, one, vlen);
  v = VFNMSUB_VV (d_lo, t, v, vlen);
  VFLOAT nlo_t = VFMADD_VV (n_lo, t, u, vlen);
  *r_hi = s;
  *r_lo = VFMADD_VV (s, v, nlo_t, vlen);
}

static inline void
RVV_ERFCF_HELPER (source_dfadd_scalar_pair)
  (VFLOAT x, VFLOAT y_hi, VFLOAT y_lo,
   VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT s = VFADD_VV (x, y_hi, vlen);
  VFLOAT t = VFADD_VV (VFSUB_VV (x, s, vlen), y_hi, vlen);
  *r_hi = s;
  *r_lo = VFADD_VV (t, y_lo, vlen);
}

static inline VFLOAT
RVV_ERFCF_HELPER (source_vpow2i) (VINT e, size_t vlen)
{
  VINT biased = __riscv_vadd (e, EXP_BIAS_F, vlen);
  return I_AS_F (__riscv_vsll (biased, MAN_LEN_F, vlen));
}

static inline VFLOAT
RVV_ERFCF_HELPER (source_vldexp2) (VFLOAT d, VINT e, size_t vlen)
{
  VINT e_half = __riscv_vsra (e, 1, vlen);
  VINT e_rest = __riscv_vsub (e, e_half, vlen);
  VFLOAT p_half = RVV_ERFCF_HELPER (source_vpow2i) (e_half, vlen);
  VFLOAT p_rest = RVV_ERFCF_HELPER (source_vpow2i) (e_rest, vlen);
  return VFMUL_VV (VFMUL_VV (d, p_half, vlen), p_rest, vlen);
}

/* expk2f: exp((d_hi, d_lo)) -> (r_hi, r_lo)
   SLEEF expk2f helper for extended precision exponential.  */
static inline void
RVV_ERFCF_HELPER (expk2f_vf2) (VFLOAT d_hi, VFLOAT d_lo,
                               VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT u = VFMUL_VF
    (VFADD_VV (d_hi, d_lo, vlen), RVV_ERFCF_R_LN2_F, vlen);
  VINT q = __riscv_vfcvt_x (u, vlen);
  VFLOAT q_flt = __riscv_vfcvt_f (q, vlen);

  /* s = d - q * ln(2) using double-float arithmetic.  */
  VFLOAT s_hi, s_lo;
  VFLOAT tmp = VFMUL_VF (q_flt, RVV_ERFCF_NEG_L2U_F, vlen);
  DFADD2_VF2_VF (d_hi, d_lo, tmp, s_hi, s_lo, vlen);
  tmp = VFMUL_VF (q_flt, RVV_ERFCF_NEG_L2L_F, vlen);
  DFADD2_VF2_VF (s_hi, s_lo, tmp, s_hi, s_lo, vlen);

  /* Polynomial approximation for exp(s) - 1.  */
  u = VFMV_VF_F (RVV_ERFCF_EXP_POLY_0_F, vlen);
  u = VFMADD_VV (u, s_hi, VFMV_VF_F (RVV_ERFCF_EXP_POLY_1_F, vlen), vlen);
  u = VFMADD_VV (u, s_hi, VFMV_VF_F (RVV_ERFCF_EXP_POLY_2_F, vlen), vlen);
  u = VFMADD_VV (u, s_hi, VFMV_VF_F (RVV_ERFCF_EXP_POLY_3_F, vlen), vlen);

  VFLOAT t_hi, t_lo;
  DFMUL_VF2_VF2_VF (s_hi, s_lo, u, t_hi, t_lo, vlen);
  DFADD2_VF2_VF_S (t_hi, t_lo, RVV_ERFCF_EXP_TAYLOR_C_F,
                   t_hi, t_lo, vlen);

  VFLOAT tmp_hi, tmp_lo;
  RVV_ERFCF_HELPER (source_dfmul_pair_pair)
    (s_hi, s_lo, t_hi, t_lo, &tmp_hi, &tmp_lo, vlen);
  DFADD2_VF2_VF_S (tmp_hi, tmp_lo, RVV_ERFCF_EXP_HALF_F,
                   t_hi, t_lo, vlen);

  VFLOAT s2_hi, s2_lo;
  DFSQU_VF2_VF2 (s_hi, s_lo, s2_hi, s2_lo, vlen);
  RVV_ERFCF_HELPER (source_dfmul_pair_pair)
    (s2_hi, s2_lo, t_hi, t_lo, &tmp_hi, &tmp_lo, vlen);
  DFADD2_VF2_VF2_VF2 (s_hi, s_lo, tmp_hi, tmp_lo, t_hi, t_lo, vlen);

  RVV_ERFCF_HELPER (source_dfadd_scalar_pair)
    (VFMV_VF_F (1.0f, vlen), t_hi, t_lo, &t_hi, &t_lo, vlen);

  t_hi = RVV_ERFCF_HELPER (source_vldexp2) (t_hi, q, vlen);
  t_lo = RVV_ERFCF_HELPER (source_vldexp2) (t_lo, q, vlen);

  /* Handle underflow: if d.x < -104, return 0.  */
  VBOOL underflow = __riscv_vmflt (d_hi, RVV_ERFCF_SOURCE_UNDERFLOW_F,
                                   vlen);
  t_hi = VFMERGE_VFM (t_hi, fp_posZero_f, underflow, vlen);
  t_lo = VFMERGE_VFM (t_lo, fp_posZero_f, underflow, vlen);

  *r_hi = t_hi;
  *r_lo = t_lo;
}

/* Select among 4 scalar float constants based on nested masks:
   if o0 then c0 else if o1 then c1 else if o2 then c2 else c3.  */
static inline VFLOAT
RVV_ERFCF_HELPER (select4f) (VBOOL o0, VBOOL o1, VBOOL o2,
                              float k0, float k1, float k2, float k3,
                             size_t vlen)
{
  VFLOAT v = VFMV_VF_F (k3, vlen);
  v = VFMERGE_VFM (v, k2, o2, vlen);
  v = VFMERGE_VFM (v, k1, o1, vlen);
  v = VFMERGE_VFM (v, k0, o0, vlen);
  return v;
}

/* Select among 4 scalar double constants converted to double-float (hi, lo).  */
static inline void
RVV_ERFCF_HELPER (select4df) (VBOOL o0, VBOOL o1, VBOOL o2,
                              double k0, double k1, double k2, double k3,
                              VFLOAT *out_hi, VFLOAT *out_lo, size_t vlen)
{
  /* Keep the high/low casts at run time, as with the data-backed values.  */
  __asm__ ("" : "+f" (k0), "+f" (k1), "+f" (k2), "+f" (k3));
  double d3 = k3;
  float d3_hi = (float) d3;
  float d3_lo = (float) (d3 - (double) d3_hi);
  VFLOAT hi = VFMV_VF_F (d3_hi, vlen);
  VFLOAT lo = VFMV_VF_F (d3_lo, vlen);

  double d2 = k2;
  float d2_hi = (float) d2;
  float d2_lo = (float) (d2 - (double) d2_hi);
  hi = VFMERGE_VFM (hi, d2_hi, o2, vlen);
  lo = VFMERGE_VFM (lo, d2_lo, o2, vlen);

  double d1 = k1;
  float d1_hi = (float) d1;
  float d1_lo = (float) (d1 - (double) d1_hi);
  hi = VFMERGE_VFM (hi, d1_hi, o1, vlen);
  lo = VFMERGE_VFM (lo, d1_lo, o1, vlen);

  double d0 = k0;
  float d0_hi = (float) d0;
  float d0_lo = (float) (d0 - (double) d0_hi);
  hi = VFMERGE_VFM (hi, d0_hi, o0, vlen);
  lo = VFMERGE_VFM (lo, d0_lo, o0, vlen);

  *out_hi = hi;
  *out_lo = lo;
}

/* Port of SLEEF's xerfcf_u15.  */

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT a_orig = x;
  VFLOAT a, vy, vy_special;
  VBOOL special_args, special_total;

  SET_ROUNDTONEAREST;

  /* Handle NaN and Inf like v_f_logf_skeleton.c (classify + mask merge).  */
  VUINT vclass = __riscv_vfclass (a_orig, vlen);
  VBOOL subnormal_args;
  IDENTIFY (vclass, class_Subnormal, subnormal_args, vlen);
  IDENTIFY (vclass, class_NaN | class_Inf, special_args, vlen);
  /* Do not feed NaN lanes to the range classifier.  The classifier is part
     of the EH prelude, so quieting its input here keeps discarded lanes from
     acquiring an implementation-dependent invalid flag.  */
  VFLOAT class_input
      = VFMERGE_VFM (a_orig, fp_posZero_f, special_args, vlen);
  VFLOAT abs_orig = __riscv_vfabs (class_input, vlen);
  VBOOL large_args = __riscv_vmfge (abs_orig, RVV_ERFCF_INTERVAL_3_F, vlen);
  large_args = __riscv_vmand (large_args, __riscv_vmnot (special_args, vlen),
                              vlen);
  special_total = __riscv_vmor (special_args, large_args, vlen);
  special_total = __riscv_vmor (special_total, subnormal_args, vlen);

  VUINT sign_bits = __riscv_vand (F_AS_U (a_orig), 0x80000000u, vlen);
  VBOOL signbit = __riscv_vmsne (sign_bits, 0, vlen);

  UINT nb_special_args = __riscv_vcpop (special_total, vlen);
  vy_special = VFMV_VF_F (fp_posZero_f, vlen);
  if (nb_special_args > 0)
    {
      VBOOL id_mask;

      /* erfcf(+Inf) = 0, erfcf(-Inf) = 2.  */
      IDENTIFY (vclass, class_negInf, id_mask, vlen);
      vy_special = VFMERGE_VFM (vy_special, 2.0f, id_mask, vlen);

      /* Propagate NaNs.  */
      IDENTIFY (vclass, class_NaN, id_mask, vlen);
      VFLOAT nan_res = VFADD_VV_M (id_mask, a_orig, a_orig, vlen);
      vy_special = VMERGE_VVM_F (vy_special, nan_res, id_mask, vlen);

      /* erfcf(-large) = 2, erfcf(+large) = 0.  */
      VBOOL large_neg = __riscv_vmand (large_args, signbit, vlen);
      vy_special = VFMERGE_VFM (vy_special, 2.0f, large_neg, vlen);

      /* erfcf(+-tiny) rounds to 1.  Use a normal-result operation to raise
         only the inexact flag for nonzero subnormals; running the rational
         core on them would also raise underflow in discarded intermediates.  */
      VFLOAT tiny_pos = VFRSUB_VF_M (subnormal_args, abs_orig, 1.0f, vlen);
      VFLOAT tiny_neg = VFADD_VF_M (subnormal_args, abs_orig, 1.0f, vlen);
      VFLOAT tiny_res = VMERGE_VVM_F (tiny_pos, tiny_neg, signbit, vlen);
      vy_special = VMERGE_VVM_F (vy_special, tiny_res, subnormal_args, vlen);
    }

  /* Substitute special and saturated args with 0 to keep intermediates
     finite.  */
  a = VFMERGE_VFM (a_orig, fp_posZero_f, special_total, vlen);

  a = __riscv_vfabs (a, vlen);

  VBOOL o0 = __riscv_vmflt (a, RVV_ERFCF_INTERVAL_0_F, vlen);
  VBOOL o1 = __riscv_vmflt (a, RVV_ERFCF_INTERVAL_1_F, vlen);
  VBOOL o2 = __riscv_vmflt (a, RVV_ERFCF_INTERVAL_2_F, vlen);
  VBOOL o3 = __riscv_vmflt (a, RVV_ERFCF_INTERVAL_3_F, vlen);

  /* u = o1 ? (a,0) : (1,0)/(a,0).  Avoid div-by-zero for a==0 lanes.  */
  VBOOL a_zero = __riscv_vmfeq (a, fp_posZero_f, vlen);
  VFLOAT a_safe = VFMERGE_VFM (a, fp_posOne_f, a_zero, vlen);
  /* The reciprocal is discarded in the o1 branch.  Keep subnormal and zero
     lanes out of that helper as well; otherwise its unmasked 1/a operation
     can raise overflow/underflow flags before u is selected from a.  */
  VFLOAT reciprocal_input = VFMERGE_VFM (a_safe, fp_posOne_f, o1, vlen);

  VFLOAT inv_hi, inv_lo;
  RVV_ERFCF_HELPER (source_dfdiv_pair_pair)
    (VFMV_VF_F (1.0f, vlen), VFMV_VF_F (0.0f, vlen), reciprocal_input,
     VFMV_VF_F (0.0f, vlen), &inv_hi, &inv_lo, vlen);

  VFLOAT u_hi = VMERGE_VVM_F (inv_hi, a, o1, vlen);
  VFLOAT u_lo
      = VMERGE_VVM_F (inv_lo, VFMV_VF_F (0.0f, vlen), o1, vlen);

  VFLOAT t = RVV_ERFCF_HELPER (select4f) (o0, o1, o2,
     RVV_ERFCF_POLY_0_0_F, RVV_ERFCF_POLY_0_1_F,
     RVV_ERFCF_POLY_0_2_F, RVV_ERFCF_POLY_0_3_F, vlen);
  t = VFMADD_VV (
      t, u_hi,
      RVV_ERFCF_HELPER (select4f) (o0, o1, o2,
        RVV_ERFCF_POLY_1_0_F, RVV_ERFCF_POLY_1_1_F,
        RVV_ERFCF_POLY_1_2_F, RVV_ERFCF_POLY_1_3_F, vlen),
      vlen);
  t = VFMADD_VV (
      t, u_hi,
      RVV_ERFCF_HELPER (select4f) (o0, o1, o2,
        RVV_ERFCF_POLY_2_0_F, RVV_ERFCF_POLY_2_1_F,
        RVV_ERFCF_POLY_2_2_F, RVV_ERFCF_POLY_2_3_F, vlen),
      vlen);
  t = VFMADD_VV (
      t, u_hi,
      RVV_ERFCF_HELPER (select4f) (o0, o1, o2,
        RVV_ERFCF_POLY_3_0_F, RVV_ERFCF_POLY_3_1_F,
        RVV_ERFCF_POLY_3_2_F, RVV_ERFCF_POLY_3_3_F, vlen),
      vlen);
  t = VFMADD_VV (
      t, u_hi,
      RVV_ERFCF_HELPER (select4f) (o0, o1, o2,
        RVV_ERFCF_POLY_4_0_F, RVV_ERFCF_POLY_4_1_F,
        RVV_ERFCF_POLY_4_2_F, RVV_ERFCF_POLY_4_3_F, vlen),
      vlen);

  VFLOAT d_hi, d_lo;
  DFMUL_VF2_VF2_VF (u_hi, u_lo, t, d_hi, d_lo, vlen);

  VFLOAT c_hi, c_lo;
  RVV_ERFCF_HELPER (select4df)
    (o0, o1, o2,
     RVV_ERFCF_DF_0_0_D, RVV_ERFCF_DF_0_1_D,
     RVV_ERFCF_DF_0_2_D, RVV_ERFCF_DF_0_3_D,
     &c_hi, &c_lo, vlen);
  DFADD2_VF2_VF2_VF2 (d_hi, d_lo, c_hi, c_lo, d_hi, d_lo, vlen);

  VFLOAT mul_hi, mul_lo;
  RVV_ERFCF_HELPER (source_dfmul_pair_pair)
    (d_hi, d_lo, u_hi, u_lo, &mul_hi, &mul_lo, vlen);
  d_hi = mul_hi;
  d_lo = mul_lo;
  RVV_ERFCF_HELPER (select4df)
    (o0, o1, o2,
     RVV_ERFCF_DF_1_0_D, RVV_ERFCF_DF_1_1_D,
     RVV_ERFCF_DF_1_2_D, RVV_ERFCF_DF_1_3_D,
     &c_hi, &c_lo, vlen);
  DFADD2_VF2_VF2_VF2 (d_hi, d_lo, c_hi, c_lo, d_hi, d_lo, vlen);

  RVV_ERFCF_HELPER (source_dfmul_pair_pair)
    (d_hi, d_lo, u_hi, u_lo, &mul_hi, &mul_lo, vlen);
  d_hi = mul_hi;
  d_lo = mul_lo;
  RVV_ERFCF_HELPER (select4df)
    (o0, o1, o2,
     RVV_ERFCF_DF_2_0_D, RVV_ERFCF_DF_2_1_D,
     RVV_ERFCF_DF_2_2_D, RVV_ERFCF_DF_2_3_D,
     &c_hi, &c_lo, vlen);
  DFADD2_VF2_VF2_VF2 (d_hi, d_lo, c_hi, c_lo, d_hi, d_lo, vlen);

  /* x = o1 ? d*a : (-a)*a, then x = o1 ? x : x + d.  */
  VFLOAT sel_hi = VFNEG_V (a, vlen);
  VFLOAT sel_lo = VFMV_VF_F (0.0f, vlen);
  sel_hi = VMERGE_VVM_F (sel_hi, d_hi, o1, vlen);
  sel_lo = VMERGE_VVM_F (sel_lo, d_lo, o1, vlen);

  VFLOAT x_hi, x_lo;
  DFMUL_VF2_VF2_VF (sel_hi, sel_lo, a, x_hi, x_lo, vlen);

  VFLOAT x_plus_d_hi, x_plus_d_lo;
  DFADD2_VF2_VF2_VF2 (x_hi, x_lo, d_hi, d_lo, x_plus_d_hi, x_plus_d_lo, vlen);
  x_hi = VMERGE_VVM_F (x_plus_d_hi, x_hi, o1, vlen);
  x_lo = VMERGE_VVM_F (x_plus_d_lo, x_lo, o1, vlen);

  VFLOAT exp_hi, exp_lo;
  RVV_ERFCF_HELPER (expk2f_vf2)
    (x_hi, x_lo, &exp_hi, &exp_lo, vlen);
  x_hi = exp_hi;
  x_lo = exp_lo;

  VFLOAT xu_hi, xu_lo;
  RVV_ERFCF_HELPER (source_dfmul_pair_pair)
    (x_hi, x_lo, u_hi, u_lo, &xu_hi, &xu_lo, vlen);
  x_hi = VMERGE_VVM_F (xu_hi, x_hi, o1, vlen);
  x_lo = VMERGE_VVM_F (xu_lo, x_lo, o1, vlen);

  VFLOAT sum = VFADD_VV (x_hi, x_lo, vlen);
  vy = VMERGE_VVM_F (VFMV_VF_F (0.0f, vlen), sum, o3, vlen);

  /* erfc(-|x|) = 2 - erfc(|x|).  */
  VFLOAT neg_res = VFSUB_VV (VFMV_VF_F (2.0f, vlen), vy, vlen);
  vy = VMERGE_VVM_F (vy, neg_res, signbit, vlen);

  /* Incorporate exceptional inputs.  */
  vy = VMERGE_VVM_F (vy, vy_special, special_total, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
