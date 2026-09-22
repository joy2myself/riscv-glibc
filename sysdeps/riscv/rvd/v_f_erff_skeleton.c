/* Single-precision vector erf function.
   erff(x) computes the error function

   Based on SLEEF xerff_u1 implementation.
   Uses two-range polynomial approximation with formula:
   erf(x) = 1 - 1/(1 + p(x))^16

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

static inline void
RVV_ERFF_POLY2DF_VF2 (VFLOAT x, VFLOAT c1, VFLOAT c0_hi, VFLOAT c0_lo,
                  VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT t_hi, t_lo;
  DFMUL_VF2_VF2_VF (c1, VFMV_VF (0.0f, vlen), x, t_hi, t_lo, vlen);
  DFADD_VF2_VF2_VF2 (c0_hi, c0_lo, t_hi, t_lo, *r_hi, *r_lo, vlen);
}

static inline void
RVV_ERFF_POLY2DF_B_VF2 (VFLOAT x, VFLOAT c1_hi, VFLOAT c1_lo,
                    VFLOAT c0_hi, VFLOAT c0_lo,
                    VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT t_hi, t_lo;
  DFMUL_VF2_VF2_VF (c1_hi, c1_lo, x, t_hi, t_lo, vlen);
  DFADD_VF2_VF2_VF2 (c0_hi, c0_lo, t_hi, t_lo, *r_hi, *r_lo, vlen);
}

static inline void
RVV_ERFF_POLY4DF_VF2 (VFLOAT x, VFLOAT c3, VFLOAT c2_hi, VFLOAT c2_lo,
                  VFLOAT c1_hi, VFLOAT c1_lo,
                  VFLOAT c0_hi, VFLOAT c0_lo,
                  VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT x2 = VFMUL_VV (x, x, vlen);
  VFLOAT p1_hi, p1_lo;
  VFLOAT p0_hi, p0_lo;
  RVV_ERFF_POLY2DF_VF2 (x, c3, c2_hi, c2_lo, &p1_hi, &p1_lo, vlen);
  RVV_ERFF_POLY2DF_B_VF2 (x, c1_hi, c1_lo, c0_hi, c0_lo,
                          &p0_hi, &p0_lo, vlen);
  VFLOAT t_hi, t_lo;
  DFMUL_VF2_VF2_VF (p1_hi, p1_lo, x2, t_hi, t_lo, vlen);
  DFADD_VF2_VF2_VF2 (p0_hi, p0_lo, t_hi, t_lo, *r_hi, *r_lo, vlen);
}

/* Source vldexp_vf_vf_vi2: retain the split exponent and four coarse-scale
   multiplies used by SLEEF.  FAST_LDEXP is a clean-only finite-core
   optimization and is intentionally not used for this source-shaped path.  */
static inline VFLOAT
RVV_ERFF_VLDEXP_VF (VFLOAT x, VINT q, size_t vlen)
{
  VINT m = __riscv_vsra (q, 31, vlen);
  VINT tmp = __riscv_vadd (m, q, vlen);
  tmp = __riscv_vsra (tmp, 6, vlen);
  tmp = __riscv_vsub (tmp, m, vlen);
  m = __riscv_vsll (tmp, 4, vlen);
  q = __riscv_vsub (q, __riscv_vsll (m, 2, vlen), vlen);

  m = __riscv_vadd (m, EXP_BIAS_F, vlen);
  VBOOL m_positive = __riscv_vmsgt (m, 0, vlen);
  m = __riscv_vmerge (VMVI_VX (0, vlen), m, m_positive, vlen);
  VBOOL m_too_large = __riscv_vmsgt (m, 0xff, vlen);
  m = __riscv_vmerge (m, 0xff, m_too_large, vlen);

  VFLOAT scale = I_AS_F (__riscv_vsll (m, MAN_LEN_F, vlen));
  x = VFMUL_VV (x, scale, vlen);
  x = VFMUL_VV (x, scale, vlen);
  x = VFMUL_VV (x, scale, vlen);
  x = VFMUL_VV (x, scale, vlen);

  scale = I_AS_F (__riscv_vsll
                  (__riscv_vadd (q, EXP_BIAS_F, vlen), MAN_LEN_F, vlen));
  return VFMUL_VV (x, scale, vlen);
}

static inline VFLOAT
RVV_ERFF_EXPKF_VF2 (VFLOAT d_hi, VFLOAT d_lo,
                size_t vlen)
{
  VFLOAT u = VFMUL_VV (VFADD_VV (d_hi, d_lo, vlen),
                       VFMV_VF (RVV_ERFF_R_LN2_F, vlen), vlen);
  VINT q = __riscv_vfcvt_x (u, vlen);
  VFLOAT q_flt = __riscv_vfcvt_f (q, vlen);

  VFLOAT s_hi, s_lo;
  VFLOAT tmp = VFMUL_VV (q_flt, VFMV_VF (-RVV_ERFF_L2U_F, vlen), vlen);
  DFADD2_VF2_VF2_VF (d_hi, d_lo, tmp, s_hi, s_lo, vlen);
  tmp = VFMUL_VV (q_flt, VFMV_VF (-RVV_ERFF_L2L_F, vlen), vlen);
  VFLOAT s_norm_hi, s_norm_lo;
  DFADD2_VF2_VF2_VF (s_hi, s_lo, tmp, s_norm_hi, s_norm_lo, vlen);
  DFNORMALIZE_VF2_VF2 (s_norm_hi, s_norm_lo, s_hi, s_lo, vlen);

  u = VFMV_VF (RVV_ERFF_EXP_POLY_0_F, vlen);
  u = VFMADD_VV (u, s_hi, VFMV_VF (RVV_ERFF_EXP_POLY_1_F, vlen), vlen);
  u = VFMADD_VV (u, s_hi, VFMV_VF (RVV_ERFF_EXP_POLY_2_F, vlen), vlen);
  u = VFMADD_VV (u, s_hi, VFMV_VF (RVV_ERFF_EXP_POLY_3_F, vlen), vlen);
  u = VFMADD_VV (u, s_hi, VFMV_VF (RVV_ERFF_EXP_POLY_4_F, vlen), vlen);

  VFLOAT s2_hi, s2_lo;
  DFSQU_VF2_VF2 (s_hi, s_lo, s2_hi, s2_lo, vlen);
  VFLOAT t_hi, t_lo;
  DFMUL_VF2_VF2_VF (s2_hi, s2_lo, u, t_hi, t_lo, vlen);
  VFLOAT sum_hi, sum_lo;
  DFADD_VF2_VF2_VF2 (s_hi, s_lo, t_hi, t_lo, sum_hi, sum_lo, vlen);
  t_hi = sum_hi;
  t_lo = sum_lo;
  VFLOAT one_hi, one_lo;
  DFADD_VF2_VF_VF2 (VFMV_VF (fp_posOne_f, vlen), t_hi, t_lo,
                    one_hi, one_lo, vlen);

  u = VFADD_VV (one_hi, one_lo, vlen);
  u = RVV_ERFF_VLDEXP_VF (u, q, vlen);

  VBOOL underflow = __riscv_vmflt (d_hi, RVV_ERFF_EXP_UNDERFLOW_F, vlen);
  u = VFMERGE_VFM (u, fp_posZero_f, underflow, vlen);
  return u;
}

/* The error function erf is an odd function: erf(-x) = -erf(x)
   For x >= 6, erf(x) rounds to 1.0 in single precision
   Uses SLEEF's two-range polynomial approximation. */

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT vx_orig = x;
  VFLOAT vx, vy;

  SET_ROUNDTONEAREST;

  /* Handle Inf and NaN */
  VUINT vclass = __riscv_vfclass (vx_orig, vlen);
  VBOOL is_nan;
  VBOOL is_inf;
  IDENTIFY (vclass, class_NaN, is_nan, vlen);
  IDENTIFY (vclass, class_Inf, is_inf, vlen);
  VBOOL special_args = VMOR (is_nan, is_inf, vlen);
  VFLOAT inf_val = __riscv_vfsgnj (VFMV_VF (fp_posOne_f, vlen), vx_orig, vlen);
  VFLOAT nan_val = VFADD_VV_M (is_nan, vx_orig, vx_orig, vlen);

  /* Keep NaN and Inf lanes out of the range comparisons.  This is the
     classify -> prepare -> sanitize boundary of the vector EH path.  */
  VFLOAT class_input
      = VFMERGE_VFM (vx_orig, fp_posZero_f, special_args, vlen);
  VFLOAT abs_orig = __riscv_vfabs (class_input, vlen);
  VBOOL eq_0 = __riscv_vmfeq (abs_orig, fp_posZero_f, vlen);
  VBOOL ge_6 = __riscv_vmfge (abs_orig, RVV_ERFF_SATURATION_THRESHOLD_F,
                              vlen);
  VBOOL core_skip = VMOR (special_args, ge_6, vlen);
  core_skip = VMOR (core_skip, eq_0, vlen);
  VFLOAT vx_safe = VFMERGE_VFM (vx_orig, fp_posZero_f, core_skip, vlen);

  /* Work with |x| */
  vx = __riscv_vfabs (vx_safe, vlen);

  /* Determine which range: x <= 2.5 or x > 2.5 */
  VBOOL le_25 = __riscv_vmfle (vx, RVV_ERFF_RANGE_SPLIT_F, vlen);

  /* Compute x^2 */
  VFLOAT x2 = __riscv_vfmul (vx, vx, vlen);
  VFLOAT x4 = __riscv_vfmul (x2, x2, vlen);

  /* Use conditional selection for coefficients based on range */
  VFLOAT p5 = VFMERGE_VFM (VFMV_VF (RVV_ERFF_POLY_GT_0_F, vlen),
                           RVV_ERFF_POLY25_0_F, le_25, vlen);
  VFLOAT p4 = VFMERGE_VFM (VFMV_VF (RVV_ERFF_POLY_GT_1_F, vlen),
                           RVV_ERFF_POLY25_1_F, le_25, vlen);
  VFLOAT p3 = VFMERGE_VFM (VFMV_VF (RVV_ERFF_POLY_GT_2_F, vlen),
                           RVV_ERFF_POLY25_2_F, le_25, vlen);
  VFLOAT p2 = VFMERGE_VFM (VFMV_VF (RVV_ERFF_POLY_GT_3_F, vlen),
                           RVV_ERFF_POLY25_3_F, le_25, vlen);
  VFLOAT p1 = VFMERGE_VFM (VFMV_VF (RVV_ERFF_POLY_GT_4_F, vlen),
                           RVV_ERFF_POLY25_4_F, le_25, vlen);
  VFLOAT p0 = VFMERGE_VFM (VFMV_VF (RVV_ERFF_POLY_GT_5_F, vlen),
                           RVV_ERFF_POLY25_5_F, le_25, vlen);

  /* Estrin-style polynomial for t */
  VFLOAT poly2 = __riscv_vfmadd (vx, p5, p4, vlen);
  VFLOAT poly4 = __riscv_vfmadd (vx, p3, p2, vlen);
  poly4 = __riscv_vfmadd (x2, poly4, __riscv_vfmadd (vx, p1, p0, vlen), vlen);
  VFLOAT t = __riscv_vfmadd (x4, poly2, poly4, vlen);

  VFLOAT c2_hi = VFMERGE_VFM (VFMV_VF (RVV_ERFF_DF_GT_0_0_F, vlen),
                              RVV_ERFF_DF25_0_0_F, le_25, vlen);
  VFLOAT c2_lo = VFMERGE_VFM (VFMV_VF (RVV_ERFF_DF_GT_0_1_F, vlen),
                              RVV_ERFF_DF25_0_1_F, le_25, vlen);
  VFLOAT c1_hi = VFMERGE_VFM (VFMV_VF (RVV_ERFF_DF_GT_1_0_F, vlen),
                              RVV_ERFF_DF25_1_0_F, le_25, vlen);
  VFLOAT c1_lo = VFMERGE_VFM (VFMV_VF (RVV_ERFF_DF_GT_1_1_F, vlen),
                              RVV_ERFF_DF25_1_1_F, le_25, vlen);
  VFLOAT c0_hi = VFMERGE_VFM (VFMV_VF (RVV_ERFF_DF_GT_2_0_F, vlen),
                              RVV_ERFF_DF25_2_0_F, le_25, vlen);
  VFLOAT c0_lo = VFMERGE_VFM (VFMV_VF (RVV_ERFF_DF_GT_2_1_F, vlen),
                              RVV_ERFF_DF25_2_1_F, le_25, vlen);

  VFLOAT t2_hi, t2_lo;
  RVV_ERFF_POLY4DF_VF2 (vx, t, c2_hi, c2_lo, c1_hi, c1_lo, c0_hi, c0_lo,
                        &t2_hi, &t2_lo, vlen);

  DFMUL_VF2_VF2_VF (t2_hi, t2_lo, vx, t2_hi, t2_lo, vlen);

  VFLOAT s2_hi, s2_lo;
  VFLOAT one = VFMV_VF (fp_posOne_f, vlen);
  DFADD_VF2_VF_VF2 (one, t2_hi, t2_lo, s2_hi, s2_lo, vlen);
  VFLOAT sq_hi, sq_lo;
  DFSQU_VF2_VF2 (s2_hi, s2_lo, sq_hi, sq_lo, vlen);
  s2_hi = sq_hi;
  s2_lo = sq_lo;
  DFSQU_VF2_VF2 (s2_hi, s2_lo, sq_hi, sq_lo, vlen);
  s2_hi = sq_hi;
  s2_lo = sq_lo;
  DFSQU_VF2_VF2 (s2_hi, s2_lo, sq_hi, sq_lo, vlen);
  s2_hi = sq_hi;
  s2_lo = sq_lo;
  DFSQU_VF2_VF2 (s2_hi, s2_lo, sq_hi, sq_lo, vlen);
  s2_hi = sq_hi;
  s2_lo = sq_lo;
  VFLOAT rec_hi, rec_lo;
  DFREC_VF2_VF2 (s2_hi, s2_lo, rec_hi, rec_lo, vlen);
  s2_hi = rec_hi;
  s2_lo = rec_lo;

  UINT nb_le_25 = __riscv_vcpop (le_25, vlen);
  if (nb_le_25 == vlen)
    {
      t2_hi = s2_hi;
      t2_lo = s2_lo;
    }
  else
    {
      VFLOAT exp_val = RVV_ERFF_EXPKF_VF2 (t2_hi, t2_lo, vlen);
      t2_hi = VMERGE_VVM_F (exp_val, s2_hi, le_25, vlen);
      t2_lo = VMERGE_VVM_F (VFMV_VF (0.0f, vlen), s2_lo, le_25, vlen);
    }

  VFLOAT t3_hi, t3_lo;
  DFADD2_VF2_VF_S (t2_hi, t2_lo, -1.0f, t3_hi, t3_lo, vlen);
  t2_hi = t3_hi;
  t2_lo = t3_lo;

  /* For very small |x| < 1e-4, use erf(x) ≈ (2/sqrt(pi)) * x */
  VBOOL very_small = __riscv_vmflt (vx, RVV_ERFF_TINY_THRESHOLD_F, vlen);
  VFLOAT small_hi, small_lo;
  DFMUL_VF2_VF2_VF (VFMV_VF (RVV_ERFF_TWO_OVER_SQRT_PI_F, vlen),
                    VFMV_VF (RVV_ERFF_TWO_OVER_SQRT_PI_LO_F, vlen),
                    vx, small_hi, small_lo, vlen);
  t2_hi = VMERGE_VVM_F (t2_hi, small_hi, very_small, vlen);
  t2_lo = VMERGE_VVM_F (t2_lo, small_lo, very_small, vlen);

  vy = VFNEG_V (VFADD_VV (t2_hi, t2_lo, vlen), vlen);

  /* For x >= 6, return 1.0 */
  vy = VFMERGE_VFM (vy, fp_posOne_f, ge_6, vlen);

  /* For x == 0, return 0.0 */
  vy = VFMERGE_VFM (vy, fp_posZero_f, eq_0, vlen);

  /* Apply sign from original x */
  vy = __riscv_vfsgnj (vy, vx_orig, vlen);

  vy = VMERGE_VVM_F (vy, inf_val, is_inf, vlen);
  vy = VMERGE_VVM_F (vy, nan_val, is_nan, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
