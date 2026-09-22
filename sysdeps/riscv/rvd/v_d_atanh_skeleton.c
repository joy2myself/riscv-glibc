/* Double-precision vector atanh function.
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

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT vx = x;
  VFLOAT vx_orig, vy, vy_special;
  VBOOL special_args;

  SET_ROUNDTONEAREST;
  vx_orig = vx;

  /* Exception handling: Inf, NaN, |x| >= 1, and |x| < 2^(-30).  */
  VUINT vclass = __riscv_vfclass (vx, vlen);
  VBOOL x_nan;
  IDENTIFY (vclass, class_NaN, x_nan, vlen);
  VFLOAT abs_x = __riscv_vfabs (vx, vlen);
  VFLOAT abs_x_cmp = VFMERGE_VFM (abs_x, fp_posZero, x_nan, vlen);
  VBOOL x_ge_1 = __riscv_vmfge (abs_x_cmp, fp_posOne, vlen);
  VBOOL x_eq_1 = __riscv_vmfeq (abs_x_cmp, fp_posOne, vlen);
  VBOOL x_gt_1 = __riscv_vmfgt (abs_x_cmp, fp_posOne, vlen);

  VUINT expo_x = __riscv_vsrl (F_AS_U (vx), MAN_LEN, vlen);
  expo_x = __riscv_vand (expo_x, 0x7FF, vlen);
  VBOOL x_small = __riscv_vmsleu (expo_x, EXP_BIAS - 31, vlen);

  special_args = __riscv_vmor (x_ge_1, x_small, vlen);
  special_args = __riscv_vmor (special_args, x_nan, vlen);
  UINT nb_special_args = __riscv_vcpop (special_args, vlen);

  vy_special = vx;
  if (nb_special_args > 0)
    {
      /* For |x| = 1: atanh(+/-1) = +/-Inf with div-by-zero signal.  */
      VFLOAT at_one = VFMV_VF (fp_posZero, vlen);
      at_one = __riscv_vfsgnj (at_one, vx, vlen);
      at_one = __riscv_vfrec7 (x_eq_1, at_one, vlen);
      vy_special = __riscv_vmerge (vy_special, at_one, x_eq_1, vlen);

      /* For |x| > 1: return NaN.  */
      VFLOAT x_tmp = __riscv_vfmerge (abs_x, fp_sNaN, x_gt_1, vlen);
      x_tmp = VFADD_VV_M (x_gt_1, x_tmp, x_tmp, vlen);
      vy_special = __riscv_vmerge (vy_special, x_tmp, x_gt_1, vlen);

      /* For small |x|: atanh(x) ~= x.  */
      VFLOAT near_zero = __riscv_vfmadd (x_small, vx, 0x1.0p-60, vx, vlen);
      vy_special = __riscv_vmerge (vy_special, near_zero, x_small, vlen);

      /* Propagate NaNs via x + x.  */
      VFLOAT x_nan_val = VFADD_VV_M (x_nan, vx, vx, vlen);
      vy_special = __riscv_vmerge (vy_special, x_nan_val, x_nan, vlen);

      vx = __riscv_vfmerge (vx, fp_posZero, special_args, vlen);
    }

  /* Work with |x|.  */
  vx = __riscv_vfsgnj (vx, fp_posOne, vlen);

  /* Get n so that 2^(-n) * (1+x)/(1-x) is roughly in [1/sqrt(2), sqrt(2)].  */
  VFLOAT one_plus_x = __riscv_vfadd (vx, fp_posOne, vlen);
  VFLOAT one_minus_x = __riscv_vfrsub (vx, fp_posOne, vlen);
  VFLOAT ratio
    = __riscv_vfmul (one_plus_x, __riscv_vfrec7 (one_minus_x, vlen), vlen);
  VUINT n = __riscv_vadd (__riscv_vsrl (F_AS_U (ratio), MAN_LEN - 8, vlen),
                          0x96, vlen);
  n = __riscv_vsub (__riscv_vsrl (n, 8, vlen), EXP_BIAS, vlen);

  /* Convert to fixed-point with scale 60.  */
  VINT x_fixed = __riscv_vfcvt_x (__riscv_vfmul (vx, 0x1.0p60, vlen), vlen);
  VINT one_plus_x_fixed = __riscv_vadd (x_fixed, RVV_ONE_Q60, vlen);
  VINT one_minus_x_fixed = __riscv_vrsub (x_fixed, RVV_ONE_Q60, vlen);
  one_minus_x_fixed = __riscv_vsll (one_minus_x_fixed, n, vlen);
  VINT numer_fixed = __riscv_vsub (one_plus_x_fixed, one_minus_x_fixed, vlen);
  VINT denom_fixed = __riscv_vadd (one_plus_x_fixed, one_minus_x_fixed, vlen);

  /* Convert back to floating-point with extra precision.  */
  VFLOAT numer = __riscv_vfcvt_f (numer_fixed, vlen);
  VINT tail = __riscv_vsub (numer_fixed, __riscv_vfcvt_x (numer, vlen), vlen);
  VFLOAT delta_numer = __riscv_vfcvt_f (tail, vlen);
  VFLOAT denom = __riscv_vfcvt_f (denom_fixed, vlen);
  tail = __riscv_vsub (denom_fixed, __riscv_vfcvt_x (denom, vlen), vlen);
  VFLOAT delta_denom = __riscv_vfcvt_f (tail, vlen);

  /* Compute r = numer/denom with extra precision.  */
  VFLOAT r_hi, r_lo, r;
  DIV2_N2D2 (numer, delta_numer, denom, delta_denom, r_hi, r_lo, vlen);

  /* For small |x| < 2^(-8), use x directly.  */
  VBOOL x_in_range = __riscv_vmflt (vx, 0x1.0p-8, vlen);
  r_hi = __riscv_vmerge (r_hi, vx, x_in_range, vlen);
  r_lo = __riscv_vfmerge (r_lo, fp_posZero, x_in_range, vlen);
  n = __riscv_vmerge (n, 0, x_in_range, vlen);

  r = __riscv_vfadd (r_hi, r_lo, vlen);
  VFLOAT rsq = __riscv_vfmul (r, r, vlen);
  VFLOAT rcube = __riscv_vfmul (rsq, r, vlen);
  VFLOAT r6 = __riscv_vfmul (rcube, rcube, vlen);

  /* Polynomial approximation for atanh.  */
  VFLOAT poly_right = PSTEP (
      0x1.c71c4a9aa397dp-4, rsq,
      PSTEP (0x1.7467d1711e0d8p-4, rsq,
             PSTEP (0x1.397813e4ac2d0p-4, 0x1.30b2960ceaa62p-4, rsq, vlen),
             vlen),
      vlen);

  VFLOAT poly_left = PSTEP (
      0x1.55555555555aep-2, rsq,
      PSTEP (0x1.999999997646fp-3, 0x1.2492494ac4a16p-3, rsq, vlen), vlen);

  VFLOAT poly = __riscv_vfmadd (poly_right, r6, poly_left, vlen);
  poly = __riscv_vfmadd (poly, rcube, r_lo, vlen);

  /* Compose final answer: (n/2)*log(2) + atanh(r).  */
  VFLOAT n_flt = __riscv_vfcvt_f (n, vlen);
  VFLOAT a = __riscv_vfmul (n_flt, RVV_LOG2_BY2_HI_D, vlen);
  VFLOAT s_hi, s_lo;
  FAST2SUM (a, r_hi, s_hi, s_lo, vlen);
  s_lo = __riscv_vfmacc (s_lo, RVV_LOG2_BY2_LO_D, n_flt, vlen);
  s_lo = __riscv_vfadd (s_lo, poly, vlen);
  vy = __riscv_vfadd (s_hi, s_lo, vlen);

  /* Restore sign.  */
  vy = __riscv_vfsgnj (vy, vx_orig, vlen);
  vy = __riscv_vmerge (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
