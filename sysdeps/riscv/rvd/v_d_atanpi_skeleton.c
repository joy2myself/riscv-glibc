/* Double-precision vector atanpi function.
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
  VFLOAT vy, vy_special;
  VBOOL special_args;

  SET_ROUNDTONEAREST;

  /* Classify special lanes up front, following the log/logf skeleton style:
     NaN propagates, Inf and very large finite inputs round to
     sign(x) * (0.5 - 1/(pi*x)), and tiny inputs use the linearized result to
     preserve the expected flags.  */
  VUINT expo_x = __riscv_vsrl (F_AS_U (vx_orig), MAN_LEN, vlen);
  expo_x = __riscv_vand (expo_x, 0x7FF, vlen);

  VBOOL expo_x_large = __riscv_vmsgeu (expo_x, EXP_BIAS + 60, vlen);
  VBOOL x_small = __riscv_vmsleu (expo_x, EXP_BIAS - 31, vlen);
  special_args = __riscv_vmor (expo_x_large, x_small, vlen);
  UINT nb_special_args = __riscv_vcpop (special_args, vlen);

  vy_special = vx_orig;
  if (nb_special_args > 0)
    {
      /* For small |x|: atanpi(x) ~= x / pi.  */
      VFLOAT near_zero
          = __riscv_vfmul (x_small, vx_orig, RVV_ONE_OVER_PI_LO_D, vlen);
      near_zero = __riscv_vfmacc (x_small, near_zero,
                                  RVV_ONE_OVER_PI_HI_D, vx_orig, vlen);
      vy_special = __riscv_vmerge (vy_special, near_zero, x_small, vlen);

      /* Preserve signed zero exactly for the current harness contract.  */
      VBOOL x_zero = __riscv_vmfeq (vx_orig, fp_posZero, vlen);
      vy_special = __riscv_vmerge (vy_special, vx_orig, x_zero, vlen);

      /* For large |x|, match veclibm's exception-generation path.  */
      VFLOAT large_sign
        = __riscv_vfsgnj (VFMV_VF (0x1.0p-1, vlen), vx_orig, vlen);
      VFLOAT large_result
        = __riscv_vfsub (expo_x_large, large_sign,
                         __riscv_vfrec7 (expo_x_large, vx_orig, vlen), vlen);
      vy_special = __riscv_vmerge (vy_special, large_result, expo_x_large,
                                   vlen);

      vx_orig = __riscv_vfmerge (vx_orig, fp_posZero, special_args, vlen);
    }

  /* atanpi(-x) = -atanpi(x), so compute sign(x)*atanpi(|x|).  */
  VFLOAT a = __riscv_vfsgnj (vx_orig, fp_posOne, vlen);
  VBOOL a_ge_one = __riscv_vmfge (a, fp_posOne, vlen);
  VBOOL a_lt_one = __riscv_vmnot (a_ge_one, vlen);

  /* For |x| >= 1: z = 1/x, else z = x.  */
  VFLOAT z = __riscv_vfrdiv (a_ge_one, a, fp_posOne, vlen);
  z = __riscv_vmerge (z, a, a_lt_one, vlen);

  /* Compute 1/a to extra precision for |x| >= 1.  */
  VFLOAT delta = VFMV_VF (fp_posOne, vlen);
  delta = __riscv_vfnmsac (a_ge_one, delta, z, a, vlen);
  delta = __riscv_vfmul (a_ge_one, delta, z, vlen);
  delta = __riscv_vfmerge (delta, fp_posZero, a_lt_one, vlen);

  /* Convert z to fixed point with scale 61+m where 2^(-m) <= |z| < 2^(-m+1).  */
  VUINT expo_61pm = __riscv_vsrl (F_AS_U (z), MAN_LEN, vlen);
  expo_61pm = __riscv_vmaxu (expo_61pm, EXP_BIAS - 60, vlen);
  expo_61pm = __riscv_vrsub (expo_61pm, 2 * EXP_BIAS + 61, vlen);
  VFLOAT scale_61pm = U_AS_F (__riscv_vsll (expo_61pm, MAN_LEN, vlen));
  VINT z_fixed = __riscv_vfcvt_x (__riscv_vfmul (z, scale_61pm, vlen), vlen);
  VINT delta_fixed = __riscv_vfcvt_x (__riscv_vfmul (delta, scale_61pm, vlen),
                                      vlen);
  delta_fixed = __riscv_vsadd (delta_fixed, z_fixed, vlen);
  z_fixed = __riscv_vmerge (z_fixed, delta_fixed, a_ge_one, vlen);

  /* V = z*z with scale 60 + 2m.  */
  VINT v = __riscv_vsmul (z_fixed, __riscv_vsll (z_fixed, 1, vlen), 1, vlen);
  /* VV = -z*z for polynomial evaluation.  */
  VINT vv = __riscv_vrsub (v, 0, vlen);

  VUINT m = __riscv_vsub (expo_61pm, EXP_BIAS + 61, vlen);
  VUINT two_m = __riscv_vsll (m, 1, vlen);
  VBOOL left_shift = __riscv_vmsltu (two_m, 3, vlen);
  VBOOL right_shift = __riscv_vmnot (left_shift, vlen);

  /* Adjust V and VV to scale 62 and 63 respectively.  */
  VINT i_tmp = __riscv_vsll (v, __riscv_vrsub (two_m, 2, vlen), vlen);
  v = __riscv_vmerge (v, i_tmp, left_shift, vlen);
  i_tmp = __riscv_vsll (vv, __riscv_vrsub (two_m, 3, vlen), vlen);
  vv = __riscv_vmerge (vv, i_tmp, left_shift, vlen);

  i_tmp = __riscv_vsra (v, __riscv_vsub (two_m, 2, vlen), vlen);
  v = __riscv_vmerge (v, i_tmp, right_shift, vlen);
  i_tmp = __riscv_vsra (vv, __riscv_vsub (two_m, 3, vlen), vlen);
  vv = __riscv_vmerge (vv, i_tmp, right_shift, vlen);

  /* WW = -z^4 in scale 63.  */
  VINT ww = __riscv_vsll (__riscv_vsmul (v, vv, 1, vlen), 1, vlen);

  /* Polynomial evaluation using fixed-point arithmetic.  */
  VINT p_even = PSTEPN_I (
      -0x56629d839b68685LL, ww,
      PSTEPN_I (
          -0x3d2984d0a6f836aLL, ww,
          PSTEPN_I (
              -0x1c5e8b5228f9fe4LL, ww,
              PSTEPN_I (-0x05deca0ae3a1a5dLL, -0x004efe42fda24d7LL, ww, vlen),
              vlen),
          vlen),
      vlen);

  p_even = PSTEPN_I (
      -0x2aaaaaaaaaaa49d3LL, ww,
      PSTEPN_I (
          -0x12492492378aaf69LL, ww,
          PSTEPN_I (
              -0xba2e88c805cbaf8LL, ww,
              PSTEPN_I (-0x888722719d1260aLL, ww,
                        PSTEPN_I (-0x6b96ef57ce79cc3LL, ww, p_even, vlen),
                        vlen),
              vlen),
          vlen),
      vlen);

  VINT p_odd = PSTEPN_I (
      0x04afe3b1345b489bLL, ww,
      PSTEPN_I (
          0x02cec355111c7439LL, ww,
          PSTEPN_I (
              0x00eaa9acebf3963eLL, ww,
              PSTEPN_I (0x001b053368ecfa14LL, 0x00006da7bb4399ddLL, ww, vlen),
              vlen),
          vlen),
      vlen);

  p_odd = PSTEPN_I (
      0x1999999999540349LL, ww,
      PSTEPN_I (
          0x0e38e38bf1671f42LL, ww,
          PSTEPN_I (
              0x09d89b293ef5f4d9LL, ww,
              PSTEPN_I (0x0786ec3df324db61LL, ww,
                        PSTEPN_I (0x060b457b3c56e750LL, ww, p_odd, vlen),
                        vlen),
              vlen),
          vlen),
      vlen);

  p_odd = __riscv_vsmul (vv, p_odd, 1, vlen);
  VINT p = __riscv_vsub (p_even, p_odd, vlen);

  p = __riscv_vsmul (vv, p, 1, vlen);
  p = __riscv_vsmul (z_fixed, p, 1, vlen);
  p = __riscv_vsub (z_fixed, p, vlen);

  /* For |x| >= 1: atanpi(x) = 0.5 - atanpi(1/x).  */
  VINT p_a_ge_one = __riscv_vsra (p, m, vlen);
  p_a_ge_one = __riscv_vrsub (p_a_ge_one, RVV_PIBY2_Q61, vlen);
  p = __riscv_vmerge (p, p_a_ge_one, a_ge_one, vlen);

  /* Multiply by 1/pi to get atanpi.  */
  p = __riscv_vsmul (p, RVV_ONE_OVER_PI_Q63, 1, vlen);

  /* Scale result back to floating point.  */
  VUINT expo_scale = __riscv_vrsub (expo_61pm, 2 * EXP_BIAS, vlen);
  expo_scale = __riscv_vmerge (expo_scale, EXP_BIAS - 61, a_ge_one, vlen);
  VFLOAT scale_result = U_AS_F (__riscv_vsll (expo_scale, MAN_LEN, vlen));
  vy = __riscv_vfcvt_f (p, vlen);
  vy = __riscv_vfmul (vy, scale_result, vlen);

  /* For very large |x|, atanpi(x) ~= sign(x) * (0.5 - 1/(pi*x)).  */
  VBOOL a_huge = __riscv_vmfgt (a, 0x1.0p27, vlen);
  VFLOAT inv = __riscv_vfrdiv (a_huge, a, fp_posOne, vlen);
  inv = VFMERGE_VFM (inv, fp_posZero, __riscv_vmnot (a_huge, vlen), vlen);
  VFLOAT inv_pi_hi
    = __riscv_vfmul (a_huge, inv, RVV_ONE_OVER_PI_HI_D, vlen);
  VFLOAT inv_pi_lo
    = __riscv_vfmul (a_huge, inv, RVV_ONE_OVER_PI_LO_D, vlen);
  VFLOAT inv_pi = __riscv_vfadd (a_huge, inv_pi_hi, inv_pi_lo, vlen);
  VFLOAT huge = __riscv_vfsub (a_huge, VFMV_VF (0x1.0p-1, vlen), inv_pi, vlen);
  huge = __riscv_vfsgnj (huge, vx_orig, vlen);
  vy = __riscv_vmerge (vy, huge, a_huge, vlen);

  /* Restore sign.  */
  vy = __riscv_vfsgnj (vy, vx_orig, vlen);

  vy = __riscv_vmerge (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
