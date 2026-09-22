/* Double-precision vector acospi function.
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
  VFLOAT vx_orig = x;
  VFLOAT vx, vy, vy_special;
  VBOOL special_args;

  SET_ROUNDTONEAREST;

  /* Exception handling for special arguments:
     |x| >= 1, Inf, NaN, and |x| < 2^(-30).  */
  VUINT expo_x = __riscv_vsrl (F_AS_U (vx_orig), MAN_LEN, vlen);
  expo_x = __riscv_vand (expo_x, 0x7FF, vlen);
  VBOOL expo_ge_bias = __riscv_vmsgeu (expo_x, EXP_BIAS, vlen);
  VBOOL expo_le_bias_m31 = __riscv_vmsleu (expo_x, EXP_BIAS - 31, vlen);
  special_args = __riscv_vmor (expo_ge_bias, expo_le_bias_m31, vlen);
  UINT nb_special_args = __riscv_vcpop (special_args, vlen);

  vy_special = vx_orig;
  if (nb_special_args > 0)
    {
      VFLOAT x_tmp = __riscv_vfsgnj (vx_orig, fp_posOne, vlen);
      VBOOL abs_x_1 = __riscv_vmand (
          special_args, __riscv_vmfeq (x_tmp, fp_posOne, vlen), vlen);
      VBOOL abs_x_gt1 = __riscv_vmand (
          special_args, __riscv_vmfgt (x_tmp, fp_posOne, vlen), vlen);

      /* For |x| > 1, return NaN.  */
      x_tmp = __riscv_vfmerge (x_tmp, fp_sNaN, abs_x_gt1, vlen);
      x_tmp = __riscv_vfadd (special_args, x_tmp, x_tmp, vlen);
      vy_special = __riscv_vmerge (vy_special, x_tmp, special_args, vlen);

      /* For |x| = 1: acospi(1) = 0, acospi(-1) = 1.  */
      VFLOAT at_one = __riscv_vfsgnjn (VFMV_VF (0x1.0p-1, vlen), vx_orig, vlen);
      at_one = __riscv_vfadd (at_one, VFMV_VF (0x1.0p-1, vlen), vlen);
      vy_special = __riscv_vmerge (vy_special, at_one, abs_x_1, vlen);

      /* For |x| < 2^(-30): acospi(x) ~= 1/2 - x/pi.  */
      VFLOAT near_zero
          = __riscv_vfnmsub (expo_le_bias_m31, vx_orig, RVV_ONE_OVER_PI_HI_D,
                             VFMV_VF (0x1.0p-1, vlen), vlen);
      vy_special
        = __riscv_vmerge (vy_special, near_zero, expo_le_bias_m31, vlen);

      vx_orig = __riscv_vfmerge (vx_orig, fp_posZero, special_args, vlen);
    }

  vx = __riscv_vfsgnj (vx_orig, fp_posOne, vlen);

  /* acospi(x) = acos(x) / pi.
     Use the same algorithm as acos and multiply the fixed-point result by
     1/pi at the end.  */
  VBOOL x_le_half = __riscv_vmfle (vx, 0x1.0p-1, vlen);
  VBOOL x_gt_half = __riscv_vmnot (x_le_half, vlen);
  VBOOL x_orig_le_half = __riscv_vmfle (vx_orig, 0x1.0p-1, vlen);
  VBOOL x_orig_lt_neghalf = __riscv_vmflt (vx_orig, 0x1.0p-1, vlen);

  VFLOAT alpha = vx;
  VFLOAT beta = VFMV_VF (0.0, vlen);
  alpha = __riscv_vfmerge (alpha, -0x1.0p-1, x_gt_half, vlen);
  beta = __riscv_vfmerge (beta, 0x1.0p-1, x_gt_half, vlen);
  VFLOAT s = __riscv_vfmadd (alpha, vx, beta, vlen);

  /* s is x*x or (1-x)/2.  */
  VINT s_fixed = __riscv_vfcvt_x (__riscv_vfmul (s, 0x1.0p63, vlen), vlen);
  VINT s_sq = __riscv_vsmul (s_fixed, s_fixed, 1, vlen);

  /* Only the x > 1/2 path needs sqrt/reciprocal work.  */
  VFLOAT sqrt_s = __riscv_vfsqrt (x_gt_half, s, vlen);
  VFLOAT delta = __riscv_vfnmsub (x_gt_half, sqrt_s, sqrt_s, s, vlen);
  delta = __riscv_vfmul (x_gt_half, delta, __riscv_vfrec7 (s, vlen), vlen);
  delta = __riscv_vfmul (x_gt_half, delta, sqrt_s, vlen);
  delta = __riscv_vfmul (x_gt_half, delta, 0x1.0p-1, vlen);
  sqrt_s = __riscv_vmerge (VFMV_VF (0.0, vlen), sqrt_s, x_gt_half, vlen);
  delta = __riscv_vmerge (VFMV_VF (0.0, vlen), delta, x_gt_half, vlen);

  /* Polynomial approximation for asin.  */
  VINT p_even = PSTEP_I (
      0x15555555555390ddLL, s_sq,
      PSTEP_I (
          0x5b6db6d09b27a82LL, s_sq,
          PSTEP_I (
              0x2dd13e6dd791f29LL, s_sq,
              PSTEP_I (0x1c6fc7fedf424bbLL, s_sq,
                       PSTEP_I (0xd5bd98b325786cLL, -0x21470ca28feec71LL, s_sq,
                                vlen),
                       vlen),
              vlen),
          vlen),
      vlen);

  VINT p_odd = PSTEP_I (
      0x99999999b7428adLL, s_sq,
      PSTEP_I (
          0x3e38e587fad54b2LL, s_sq,
          PSTEP_I (
              0x238d7e0436a1c30LL, s_sq,
              PSTEP_I (0x18ecc06c5a390e3LL, s_sq,
                       PSTEP_I (0x28063c8b4b6a072LL, 0x41646ebd6edd35eLL, s_sq,
                                vlen),
                       vlen),
              vlen),
          vlen),
      vlen);

  p_odd = __riscv_vsmul (p_odd, s_fixed, 1, vlen);
  VINT poly = __riscv_vsadd (p_odd, p_even, vlen);
  poly = __riscv_vsmul (poly, s_fixed, 1, vlen);

  VFLOAT r = __riscv_vmerge (vx, sqrt_s, x_gt_half, vlen);

  VINT m = U_AS_I (__riscv_vrsub (__riscv_vsrl (F_AS_U (r), MAN_LEN, vlen),
                                  EXP_BIAS, vlen));
  m = __riscv_vmin (m, 60, vlen);
  VINT q = __riscv_vadd (m, 60, vlen);

  q = __riscv_vmerge (q, 60, x_orig_le_half, vlen);
  r = __riscv_vfsgnjx (r, vx_orig, vlen);
  delta = __riscv_vfsgnjx (delta, vx_orig, vlen);

  VFLOAT scale_r = U_AS_F (__riscv_vsll (
      I_AS_U (__riscv_vadd (q, EXP_BIAS, vlen)), MAN_LEN, vlen));
  VINT r_fixed = __riscv_vfcvt_x (__riscv_vfmul (r, scale_r, vlen), vlen);
  r_fixed = __riscv_vsadd (
      __riscv_vfcvt_x (__riscv_vfmul (delta, scale_r, vlen), vlen), r_fixed,
      vlen);
  poly = __riscv_vsadd (r_fixed, __riscv_vsmul (poly, r_fixed, 1, vlen), vlen);
  VINT poly_prime = __riscv_vsadd (poly, poly, vlen);

  poly = __riscv_vrsub (poly, 0, vlen);
  poly = __riscv_vmerge (poly, poly_prime, x_gt_half, vlen);

  /* Add constant based on input range.  */
  VINT c = VMVI_VX (0, vlen);
  c = __riscv_vmerge (c, RVV_PI_Q60, x_orig_lt_neghalf, vlen);
  c = __riscv_vmerge (c, RVV_PIBY2_Q60, x_le_half, vlen);
  poly = __riscv_vsadd (c, poly, vlen);

  /* Multiply by 1/pi to get acospi.  */
  poly = __riscv_vsmul (poly, RVV_ONE_OVER_PI_Q63, 1, vlen);

  VFLOAT inv_scale_r = U_AS_F (__riscv_vsll (
      I_AS_U (__riscv_vrsub (q, EXP_BIAS, vlen)), MAN_LEN, vlen));
  vy = __riscv_vfmul (inv_scale_r, __riscv_vfcvt_f (poly, vlen), vlen);
  vy = __riscv_vmerge (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
