/* Single-precision vector acos function with ~1 ULP accuracy.
   Based on SLEEF xacosf_u1 implementation.

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

/* The pinned RVV SLEEF build selects ENABLE_FMA_SP.  Keep its
   dfsqrt_vf2_vf operation graph local to this skeleton; the generic clean
   helper uses a different residual/reconstruction order.  */
static inline void
RVV_ACOSF_SOURCE_MANGLE (acosf_source_dfsqrt)
  (VFLOAT d, VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT t = __riscv_vfsqrt (d, vlen);

  VFLOAT tsq_hi, tsq_lo;
  RVV_SLEEF_ORIG_DFMUL_VF_VF (t, t, tsq_hi, tsq_lo, vlen);

  VFLOAT sum_hi, sum_lo;
  RVV_SLEEF_ORIG_DFADD2_VF_VF2
    (d, tsq_hi, tsq_lo, sum_hi, sum_lo, vlen);

  VFLOAT rec_hi = VFDIV_VV (VFMV_VF (1.0f, vlen), t, vlen);
  VFLOAT rec_lo = VFMUL_VV
    (rec_hi, VFNMSUB_VV (t, rec_hi, VFMV_VF (1.0f, vlen), vlen), vlen);

  VFLOAT q_hi, q_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (sum_hi, sum_lo, rec_hi, rec_lo, q_hi, q_lo, vlen);
  RVV_SLEEF_ORIG_DFSCALE_VF2_VF
    (q_hi, q_lo, VFMV_VF (0.5f, vlen), *r_hi, *r_lo, vlen);
}
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT vx, vy, vy_special;
  VBOOL special_args;
  VUINT ux;
  VBOOL x_nan, x_neg;

  SET_ROUNDTONEAREST;
  vx = x;
  ux = F_AS_U (vx);

  VUINT vclass = __riscv_vfclass (vx, vlen);
  IDENTIFY_F (vclass, class_NaN, x_nan, vlen);
  x_neg = __riscv_vmsne (__riscv_vand (ux, 0x80000000u, vlen), 0, vlen);

  /* Handle special cases: |x| > 1 returns NaN, |x| = 1 returns 0 or pi */
  VFLOAT abs_x = __riscv_vfsgnj (vx, fp_posOne, vlen);
  VFLOAT abs_x_cmp = __riscv_vfmerge (abs_x, fp_posZero, x_nan, vlen);
  VBOOL abs_x_gt1 = __riscv_vmfgt (abs_x_cmp, fp_posOne, vlen);
  VBOOL abs_x_eq1 = __riscv_vmfeq (abs_x_cmp, fp_posOne, vlen);
  special_args = __riscv_vmor (abs_x_gt1, abs_x_eq1, vlen);
  special_args = __riscv_vmor (special_args, x_nan, vlen);
  UINT nb_special_args = __riscv_vcpop (special_args, vlen);

  vy_special = vx;
  if (nb_special_args > 0)
    {
      /* Quiet invalid lanes with x+x after substituting sNaN.  */
      VFLOAT invalid_result = __riscv_vfmerge (abs_x_cmp, fp_sNaN, abs_x_gt1,
                                               vlen);
      invalid_result = __riscv_vfadd (abs_x_gt1, invalid_result, invalid_result,
                                      vlen);
      vy_special = __riscv_vmerge (vy_special, invalid_result, abs_x_gt1, vlen);

      VFLOAT vy_nan = __riscv_vfadd (x_nan, vx, vx, vlen);
      vy_special = __riscv_vmerge (vy_special, vy_nan, x_nan, vlen);

      /* For x = 1, acos(1) = 0 */
      VBOOL x_eq_pos1 = __riscv_vmandn (abs_x_eq1, x_neg, vlen);
      vy_special = __riscv_vfmerge (vy_special, fp_posZero, x_eq_pos1, vlen);

      /* For x = -1, acos(-1) = pi */
      VBOOL x_eq_neg1 = __riscv_vmand (abs_x_eq1, x_neg, vlen);
      vy_special = __riscv_vfmerge
        (vy_special, RVV_ACOSF_PI_HI_F, x_eq_neg1, vlen);

      /* Replace special args with safe values for main computation */
      vx = __riscv_vfmerge (vx, fp_posZero, special_args, vlen);
      abs_x = __riscv_vfsgnj (vx, fp_posOne, vlen);
    }

  VBOOL main_args = __riscv_vmnot (special_args, vlen);
  /* Do not run the discarded finite approximation when every lane is
     exceptional.  Its sanitized arithmetic can still set FE_INEXACT, which
     the final special-result merge cannot remove.  */
  if (__riscv_vcpop (main_args, vlen) == 0)
    {
      RESTORE_FRM;
      return vy_special;
    }

  /* Algorithm:
     For |x| <= 0.5: acos(x) = pi/2 - asin(x)
                             = pi/2 - (x + x*x2*P(x2))
     For |x| > 0.5:  acos(x) = 2*asin(sqrt((1-|x|)/2))
                             = 2*(sqrt(z) + sqrt(z)*z*P(z)) for x > 0
                             = pi - 2*(sqrt(z) + sqrt(z)*z*P(z)) for x < 0
     where z = (1-|x|)/2

     For _u1 accuracy, we use double-float arithmetic for key operations. */

  VBOOL x_lt_half = __riscv_vmflt (abs_x, 0.5f, vlen);
  VBOOL x_ge_half = __riscv_vmnot (x_lt_half, vlen);

  /* Compute x2 = x*x for |x| <= 0.5, or z = (1-|x|)/2 for |x| > 0.5.
     For |x| < 2^-63, x*x is below the normal range and the correction
     x2*x is far below the result precision.  Quarantine that intermediate
     product so subnormal inputs do not manufacture FE_UNDERFLOW; the
     source x itself remains in x_hi and therefore still participates in the
     leading pi/2 - x term.  */
  VBOOL x_tiny = __riscv_vmflt (abs_x, 0x1.0p-63f, vlen);
  VFLOAT x_for_square = __riscv_vfmerge (vx, fp_posZero, x_tiny, vlen);
  VFLOAT x2 = __riscv_vfmul (x_for_square, x_for_square, vlen);
  VFLOAT one_minus_abs_x
      = __riscv_vfsub (VFMV_VF (fp_posOne, vlen), abs_x, vlen);
  VFLOAT z = __riscv_vfmul (one_minus_abs_x, 0x1.0p-1f, vlen);
  x2 = __riscv_vmerge (x2, z, x_ge_half, vlen);

  /* Keep the source selection while quarantining the unused sqrt input on
     the small branch and on sanitized special lanes.  */
  VFLOAT x2_sqrt = __riscv_vfmerge (x2, fp_posOne, x_lt_half, vlen);
  VFLOAT sqrt_z_hi, sqrt_z_lo;
  RVV_ACOSF_SOURCE_MANGLE (acosf_source_dfsqrt)
    (x2_sqrt, &sqrt_z_hi, &sqrt_z_lo, vlen);

  /* Select base value: (abs_x, 0) for |x| <= 0.5, (sqrt_z_hi, sqrt_z_lo) for
   * |x| > 0.5 */
  VFLOAT x_hi = abs_x;
  VFLOAT x_lo = VFMV_VF (0.0f, vlen);
  x_hi = __riscv_vmerge (x_hi, sqrt_z_hi, x_ge_half, vlen);
  x_lo = __riscv_vmerge (x_lo, sqrt_z_lo, x_ge_half, vlen);

  /* For |x| = 1, set x_hi = x_lo = 0 to avoid issues */
  x_hi = __riscv_vfmerge (x_hi, fp_posZero,
                          __riscv_vmfeq (abs_x, fp_posOne, vlen), vlen);
  x_lo = __riscv_vfmerge (x_lo, fp_posZero,
                          __riscv_vmfeq (abs_x, fp_posOne, vlen), vlen);

  /* Polynomial approximation: P(x2) = C4 + x2*(C3 + x2*(C2 + x2*(C1 +
   * x2*C0))) */
  VFLOAT u = VFMV_VF (RVV_ACOSF_C4, vlen);
  u = __riscv_vfmacc (VFMV_VF (RVV_ACOSF_C3, vlen), u, x2, vlen);
  u = __riscv_vfmacc (VFMV_VF (RVV_ACOSF_C2, vlen), u, x2, vlen);
  u = __riscv_vfmacc (VFMV_VF (RVV_ACOSF_C1, vlen), u, x2, vlen);
  u = __riscv_vfmacc (VFMV_VF (RVV_ACOSF_C0, vlen), u, x2, vlen);

  /* u = x2 * x_hi * u (polynomial correction term) */
  u = __riscv_vfmul (u, __riscv_vfmul (x2, x_hi, vlen), vlen);

  /* For |x| <= 0.5:
     Compute pi/2 - (sign(x)*x_hi + sign(x)*u) using double-float subtraction
     y_hi = pi/2_hi - sign(x)*(x_hi + u)
     y_lo = pi/2_lo (absorbed into final result) */
  VFLOAT signed_x_hi = __riscv_vfsgnj (x_hi, vx, vlen);
  VFLOAT signed_u = __riscv_vfsgnj (u, vx, vlen);

  VFLOAT small_term_hi, small_term_lo;
  RVV_SLEEF_ORIG_DFADD_VF_VF
    (signed_x_hi, signed_u, small_term_hi, small_term_lo, vlen);
  VFLOAT y_hi_small, y_lo_small;
  RVV_SLEEF_ORIG_DFSUB_VF2_VF2
    (VFMV_VF (RVV_ACOSF_PIBY2_HI_F, vlen),
     VFMV_VF (RVV_ACOSF_PIBY2_LO_F, vlen),
     small_term_hi, small_term_lo, y_hi_small, y_lo_small, vlen);
  VFLOAT y_small = __riscv_vfadd (y_hi_small, y_lo_small, vlen);

  /* For |x| > 0.5:
     Compute 2*(x_hi + x_lo + u) for x > 0
     Compute pi - 2*(x_hi + x_lo + u) for x < 0

     Using double-float: (x_hi, x_lo) + u, then scale by 2 */
  VFLOAT sum_hi, sum_lo;
  RVV_SLEEF_ORIG_DFADD_VF2_VF2_VF
    (x_hi, x_lo, u, sum_hi, sum_lo, vlen);

  VFLOAT two_sum_hi, two_sum_lo;
  RVV_SLEEF_ORIG_DFSCALE_VF2_VF
    (sum_hi, sum_lo, VFMV_VF (2.0f, vlen), two_sum_hi, two_sum_lo, vlen);

  VBOOL x_lt_zero = __riscv_vmflt (vx, fp_posZero, vlen);
  VBOOL x_ge_half_neg = __riscv_vmand (x_ge_half, x_lt_zero, vlen);

  /* For x > 0.5: y = 2*(x_hi + x_lo + u) */
  VFLOAT y_large = __riscv_vfadd (two_sum_hi, two_sum_lo, vlen);

  /* For x < -0.5: y = (pi_hi, pi_lo) - (two_sum_hi, two_sum_lo) */
  VFLOAT y_neg_hi, y_neg_lo;
  RVV_SLEEF_ORIG_DFSUB_VF2_VF2
    (VFMV_VF (RVV_ACOSF_PI_HI_F, vlen),
     VFMV_VF (RVV_ACOSF_PI_LO_F, vlen), two_sum_hi,
     two_sum_lo, y_neg_hi, y_neg_lo, vlen);
  VFLOAT y_neg = __riscv_vfadd (y_neg_hi, y_neg_lo, vlen);
  y_large = __riscv_vmerge (y_large, y_neg, x_ge_half_neg, vlen);

  /* Select result based on |x| <= 0.5 or |x| > 0.5 */
  vy = __riscv_vmerge (y_large, y_small, x_lt_half, vlen);

  /* Merge special cases */
  vy = __riscv_vmerge (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
