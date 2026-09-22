/* Single-precision source-original SLEEF erfcf skeleton.

   This is a target-owned translation of SLEEF xerfcf_u15 and its expk2f
   helper from the pinned source tree.  The four source intervals, selected
   polynomial, double-float operation order, source vldexp2 scaling, and
   source-native late exceptional selection are intentionally kept separate
   from the clean erfcf implementation.

   Only the RVV representation, LMUL/VL plumbing, and glibc importer ABI are
   adapted.  No SLEEF source file is included, compiled, or linked.

   SLEEF's source build uses the FMA branch with an explicit
   -ffp-contract=off contract; the wrapper's build rule preserves that
   contract for this comparison implementation.

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

/* SLEEF vpow2i_vf_vi2.  The exponent range used by expk2f is within the
   source helper's direct bit-construction range.  */
static inline VFLOAT
RVV_SLEEF_ORIG_ERFCF_HELPER (vpow2i) (VINT e, size_t vlen)
{
  VINT biased = __riscv_vadd (e, 0x7f, vlen);
  return I_AS_F (__riscv_vsll (biased, 23, vlen));
}

/* SLEEF vldexp2_vf_vf_vi2: two source-ordered power-of-two multiplies.  */
static inline VFLOAT
RVV_SLEEF_ORIG_ERFCF_HELPER (vldexp2) (VFLOAT d, VINT e, size_t vlen)
{
  VINT e_half = __riscv_vsra (e, 1, vlen);
  VINT e_rest = __riscv_vsub (e, e_half, vlen);
  VFLOAT p_half = RVV_SLEEF_ORIG_ERFCF_HELPER (vpow2i) (e_half, vlen);
  VFLOAT p_rest = RVV_SLEEF_ORIG_ERFCF_HELPER (vpow2i) (e_rest, vlen);
  return VFMUL_VV (VFMUL_VV (d, p_half, vlen), p_rest, vlen);
}

/* Source expk2f, including the source non-compensated 1+t reconstruction and
   the late underflow clearing operation.  */
static inline void
RVV_SLEEF_ORIG_ERFCF_HELPER (expk2f) (VFLOAT d_hi, VFLOAT d_lo,
                                      VFLOAT *r_hi, VFLOAT *r_lo,
                                      size_t vlen)
{
  VFLOAT u = VFMUL_VF
    (VFADD_VV (d_hi, d_lo, vlen), RVV_ERFCF_R_LN2_F, vlen);
  VINT q = RVV_SLEEF_ORIG_ERFCF_VFCVT_X_RNE (u, vlen);

  VFLOAT q_f = __riscv_vfcvt_f (q, vlen);
  VFLOAT s_hi, s_lo;
  VFLOAT q_l2u = VFMUL_VF
    (q_f, RVV_ERFCF_NEG_L2U_F, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF
    (d_hi, d_lo, q_l2u, s_hi, s_lo, vlen);
  VFLOAT q_l2l = VFMUL_VF
    (q_f, RVV_ERFCF_NEG_L2L_F, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF
    (s_hi, s_lo, q_l2l, s_hi, s_lo, vlen);

  u = VFMV_VF (RVV_ERFCF_EXP_POLY_0_F, vlen);
  u = VFMADD_VV (u, s_hi, VFMV_VF (RVV_ERFCF_EXP_POLY_1_F, vlen), vlen);
  u = VFMADD_VV (u, s_hi, VFMV_VF (RVV_ERFCF_EXP_POLY_2_F, vlen), vlen);
  u = VFMADD_VV (u, s_hi, VFMV_VF (RVV_ERFCF_EXP_POLY_3_F, vlen), vlen);

  VFLOAT t_hi, t_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (s_hi, s_lo, u, t_hi, t_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF_S
    (t_hi, t_lo, RVV_ERFCF_EXP_TAYLOR_C_F,
     t_hi, t_lo, vlen);

  VFLOAT tmp_hi, tmp_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (s_hi, s_lo, t_hi, t_lo, tmp_hi, tmp_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF_S
    (tmp_hi, tmp_lo, RVV_ERFCF_EXP_HALF_F, t_hi, t_lo, vlen);

  VFLOAT s2_hi, s2_lo;
  RVV_SLEEF_ORIG_DFSQU_VF2
    (s_hi, s_lo, s2_hi, s2_lo, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (s2_hi, s2_lo, t_hi, t_lo, tmp_hi, tmp_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF2
    (s_hi, s_lo, tmp_hi, tmp_lo, t_hi, t_lo, vlen);

  /* Source dfadd_vf2_vf_vf2(1,t), deliberately not dfadd2.  */
  RVV_SLEEF_ORIG_DFADD_VF_VF2
    (VFMV_VF (1.0f, vlen), t_hi, t_lo, t_hi, t_lo, vlen);

  t_hi = RVV_SLEEF_ORIG_ERFCF_HELPER (vldexp2) (t_hi, q, vlen);
  t_lo = RVV_SLEEF_ORIG_ERFCF_HELPER (vldexp2) (t_lo, q, vlen);

  VBOOL underflow = __riscv_vmflt (d_hi, RVV_ERFCF_SOURCE_UNDERFLOW_F,
                                   vlen);
  t_hi = VFMERGE_VFM (t_hi, fp_posZero_f, underflow, vlen);
  t_lo = VFMERGE_VFM (t_lo, fp_posZero_f, underflow, vlen);

  *r_hi = t_hi;
  *r_lo = t_lo;
}

/* Source vsel_vf_vo_vo_vo_f_f_f_f for the four interval constants.  */
static inline VFLOAT
RVV_SLEEF_ORIG_ERFCF_HELPER (select4f) (VBOOL o0, VBOOL o1, VBOOL o2,
                              float k0, float k1, float k2, float k3, size_t vlen)
{
  VFLOAT v = VFMV_VF (k3, vlen);
  v = VFMERGE_VFM (v, k2, o2, vlen);
  v = VFMERGE_VFM (v, k1, o1, vlen);
  v = VFMERGE_VFM (v, k0, o0, vlen);
  return v;
}

/* Source vsel_vf2_vo_vo_vo_d_d_d_d.  The source double constants are split
   into the same binary32 high/low pair used by the target-owned DF helpers.  */
static inline void
RVV_SLEEF_ORIG_ERFCF_HELPER (select4df) (VBOOL o0, VBOOL o1, VBOOL o2,
                              double k0, double k1, double k2, double k3, VFLOAT *out_hi,
                                         VFLOAT *out_lo, size_t vlen)
{
  /* Keep the high/low casts at run time, as with the data-backed values.  */
  __asm__ ("" : "+f" (k0), "+f" (k1), "+f" (k2), "+f" (k3));
  double c3 = k3;
  float c3_hi = (float) c3;
  float c3_lo = (float) (c3 - (double) c3_hi);
  VFLOAT hi = VFMV_VF (c3_hi, vlen);
  VFLOAT lo = VFMV_VF (c3_lo, vlen);

  double c2 = k2;
  float c2_hi = (float) c2;
  float c2_lo = (float) (c2 - (double) c2_hi);
  hi = VFMERGE_VFM (hi, c2_hi, o2, vlen);
  lo = VFMERGE_VFM (lo, c2_lo, o2, vlen);

  double c1 = k1;
  float c1_hi = (float) c1;
  float c1_lo = (float) (c1 - (double) c1_hi);
  hi = VFMERGE_VFM (hi, c1_hi, o1, vlen);
  lo = VFMERGE_VFM (lo, c1_lo, o1, vlen);

  double c0 = k0;
  float c0_hi = (float) c0;
  float c0_lo = (float) (c0 - (double) c0_hi);
  hi = VFMERGE_VFM (hi, c0_hi, o0, vlen);
  lo = VFMERGE_VFM (lo, c0_lo, o0, vlen);

  *out_hi = hi;
  *out_lo = lo;
}

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  /* The source SLEEF entry uses its source-native late special selection;
     there is intentionally no clean-side quarantine or early merge here.  */
  VFLOAT s = x;
  VFLOAT a = __riscv_vfabs (x, vlen);
  VFLOAT zero = VFMV_VF (0.0f, vlen);

  VBOOL o0 = __riscv_vmflt (a, RVV_ERFCF_INTERVAL_0_F, vlen);
  VBOOL o1 = __riscv_vmflt (a, RVV_ERFCF_INTERVAL_1_F, vlen);
  VBOOL o2 = __riscv_vmflt (a, RVV_ERFCF_INTERVAL_2_F, vlen);
  VBOOL o3 = __riscv_vmflt (a, RVV_ERFCF_INTERVAL_3_F, vlen);

  VFLOAT inv_hi, inv_lo;
  RVV_SLEEF_ORIG_DFDIV_VF2_VF2
    (VFMV_VF (1.0f, vlen), zero, a, zero, inv_hi, inv_lo, vlen);

  VFLOAT u_hi = VMERGE_VVM_F (inv_hi, a, o1, vlen);
  VFLOAT u_lo = VMERGE_VVM_F (inv_lo, zero, o1, vlen);

  VFLOAT t = RVV_SLEEF_ORIG_ERFCF_HELPER (select4f)
    (o0, o1, o2,
     RVV_ERFCF_POLY_0_0_F, RVV_ERFCF_POLY_0_1_F,
     RVV_ERFCF_POLY_0_2_F, RVV_ERFCF_POLY_0_3_F, vlen);
  t = VFMADD_VV
    (t, u_hi, RVV_SLEEF_ORIG_ERFCF_HELPER (select4f)
       (o0, o1, o2,
     RVV_ERFCF_POLY_1_0_F, RVV_ERFCF_POLY_1_1_F,
     RVV_ERFCF_POLY_1_2_F, RVV_ERFCF_POLY_1_3_F, vlen), vlen);
  t = VFMADD_VV
    (t, u_hi, RVV_SLEEF_ORIG_ERFCF_HELPER (select4f)
       (o0, o1, o2,
     RVV_ERFCF_POLY_2_0_F, RVV_ERFCF_POLY_2_1_F,
     RVV_ERFCF_POLY_2_2_F, RVV_ERFCF_POLY_2_3_F, vlen), vlen);
  t = VFMADD_VV
    (t, u_hi, RVV_SLEEF_ORIG_ERFCF_HELPER (select4f)
       (o0, o1, o2,
     RVV_ERFCF_POLY_3_0_F, RVV_ERFCF_POLY_3_1_F,
     RVV_ERFCF_POLY_3_2_F, RVV_ERFCF_POLY_3_3_F, vlen), vlen);
  t = VFMADD_VV
    (t, u_hi, RVV_SLEEF_ORIG_ERFCF_HELPER (select4f)
       (o0, o1, o2,
     RVV_ERFCF_POLY_4_0_F, RVV_ERFCF_POLY_4_1_F,
     RVV_ERFCF_POLY_4_2_F, RVV_ERFCF_POLY_4_3_F, vlen), vlen);

  VFLOAT d_hi, d_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (u_hi, u_lo, t, d_hi, d_lo, vlen);

  VFLOAT c_hi, c_lo;
  RVV_SLEEF_ORIG_ERFCF_HELPER (select4df)
    (o0, o1, o2,
     RVV_ERFCF_DF_0_0_D, RVV_ERFCF_DF_0_1_D,
     RVV_ERFCF_DF_0_2_D, RVV_ERFCF_DF_0_3_D,
     &c_hi, &c_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF2
    (d_hi, d_lo, c_hi, c_lo, d_hi, d_lo, vlen);

  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (d_hi, d_lo, u_hi, u_lo, d_hi, d_lo, vlen);
  RVV_SLEEF_ORIG_ERFCF_HELPER (select4df)
    (o0, o1, o2,
     RVV_ERFCF_DF_1_0_D, RVV_ERFCF_DF_1_1_D,
     RVV_ERFCF_DF_1_2_D, RVV_ERFCF_DF_1_3_D,
     &c_hi, &c_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF2
    (d_hi, d_lo, c_hi, c_lo, d_hi, d_lo, vlen);

  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (d_hi, d_lo, u_hi, u_lo, d_hi, d_lo, vlen);
  RVV_SLEEF_ORIG_ERFCF_HELPER (select4df)
    (o0, o1, o2,
     RVV_ERFCF_DF_2_0_D, RVV_ERFCF_DF_2_1_D,
     RVV_ERFCF_DF_2_2_D, RVV_ERFCF_DF_2_3_D,
     &c_hi, &c_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF2
    (d_hi, d_lo, c_hi, c_lo, d_hi, d_lo, vlen);

  VFLOAT selected_hi = VFNEG_V (a, vlen);
  VFLOAT selected_lo = zero;
  selected_hi = VMERGE_VVM_F (selected_hi, d_hi, o1, vlen);
  selected_lo = VMERGE_VVM_F (selected_lo, d_lo, o1, vlen);

  VFLOAT x_hi, x_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (selected_hi, selected_lo, a, x_hi, x_lo, vlen);

  VFLOAT x_plus_d_hi, x_plus_d_lo;
  RVV_SLEEF_ORIG_DFADD2_VF2_VF2
    (x_hi, x_lo, d_hi, d_lo, x_plus_d_hi, x_plus_d_lo, vlen);
  x_hi = VMERGE_VVM_F (x_plus_d_hi, x_hi, o1, vlen);
  x_lo = VMERGE_VVM_F (x_plus_d_lo, x_lo, o1, vlen);

  VFLOAT exp_hi, exp_lo;
  RVV_SLEEF_ORIG_ERFCF_HELPER (expk2f)
    (x_hi, x_lo, &exp_hi, &exp_lo, vlen);
  x_hi = exp_hi;
  x_lo = exp_lo;

  VFLOAT xu_hi, xu_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (x_hi, x_lo, u_hi, u_lo, xu_hi, xu_lo, vlen);
  x_hi = VMERGE_VVM_F (xu_hi, x_hi, o1, vlen);
  x_lo = VMERGE_VVM_F (xu_lo, x_lo, o1, vlen);

  VFLOAT sum = VFADD_VV (x_hi, x_lo, vlen);
  VFLOAT result = VMERGE_VVM_F (zero, sum, o3, vlen);
  VBOOL signbit = __riscv_vmsne
    (__riscv_vand (F_AS_U (s), 0x80000000u, vlen), 0, vlen);
  VFLOAT neg_result = VFSUB_VV (VFMV_VF (2.0f, vlen), result, vlen);
  result = VMERGE_VVM_F (result, neg_result, signbit, vlen);

  VBOOL input_nan = __riscv_vmfne (s, s, vlen);
  result = VFMERGE_VFM (result, fp_qNaN_f, input_nan, vlen);

  return result;
}

#undef LMUL
#undef SIMDLEN
