/* Double-precision vector atan2 function.
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
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT vy, VFLOAT vx, size_t vlen)
{
  VFLOAT vw, vw_special;
  VUINT vclass_y, vclass_x;
  UINT stencil, class_of_interest;
  VBOOL special_y, special_x, special_args, id_mask;
  UINT nb_special_args;

  SET_ROUNDTONEAREST;

  stencil = class_NaN | class_Inf | class_Zero;
  vclass_y = __riscv_vfclass (vy, vlen);
  IDENTIFY (vclass_y, stencil, special_y, vlen);
  vclass_x = __riscv_vfclass (vx, vlen);
  IDENTIFY (vclass_x, stencil, special_x, vlen);
  special_args = __riscv_vmor (special_y, special_x, vlen);
  nb_special_args = __riscv_vcpop (special_args, vlen);

  vw_special = VFMV_VF (fp_posZero, vlen);
  if (nb_special_args > 0)
    {
      VFLOAT result_tmp;

      class_of_interest = class_NaN;
      IDENTIFY (vclass_y, class_of_interest, special_y, vlen);
      IDENTIFY (vclass_x, class_of_interest, special_x, vlen);
      VBOOL y_notNaN = __riscv_vmnot (special_y, vlen);
      id_mask = __riscv_vmor (special_y, special_x, vlen);
      result_tmp = VFADD_VV_M (id_mask, vy, vx, vlen);
      vw_special = __riscv_vmerge (vw_special, result_tmp, id_mask, vlen);

      VBOOL x_0, x_neg0, x_pos0;
      class_of_interest = class_posZero;
      IDENTIFY (vclass_x, class_of_interest, x_pos0, vlen);
      class_of_interest = class_negZero;
      IDENTIFY (vclass_x, class_of_interest, x_neg0, vlen);
      x_0 = __riscv_vmor (x_pos0, x_neg0, vlen);

      VBOOL y_0, y_not0;
      class_of_interest = class_Zero;
      IDENTIFY (vclass_y, class_of_interest, y_0, vlen);
      y_not0 = __riscv_vmnot (y_0, vlen);
      y_not0 = __riscv_vmand (y_not0, y_notNaN, vlen);
      id_mask = __riscv_vmand (x_0, y_not0, vlen);
      result_tmp = VFMV_VF (RVV_ATAN2_Y_NONZERO_X_ZERO, vlen);
      result_tmp = __riscv_vfsgnj (result_tmp, vy, vlen);
      vw_special = __riscv_vmerge (vw_special, result_tmp, id_mask, vlen);

      id_mask = __riscv_vmand (y_0, x_pos0, vlen);
      vw_special = __riscv_vmerge (vw_special, vy, id_mask, vlen);

      id_mask = __riscv_vmand (y_0, x_neg0, vlen);
      result_tmp = VFMV_VF (RVV_ATAN2_Y_ZERO_X_NEGATIVE, vlen);
      result_tmp = __riscv_vfsgnj (result_tmp, vy, vlen);
      vw_special = __riscv_vmerge (vw_special, result_tmp, id_mask, vlen);

      VBOOL x_posInf, x_negInf, y_Inf, y_finite;
      class_of_interest = class_Inf;
      IDENTIFY (vclass_y, class_of_interest, y_Inf, vlen);
      y_finite = __riscv_vmandn (y_notNaN, y_Inf, vlen);
      x_posInf = __riscv_vmfeq (vx, fp_posInf, vlen);
      x_negInf = __riscv_vmfeq (vx, fp_negInf, vlen);

      id_mask = __riscv_vmand (x_posInf, y_Inf, vlen);
      result_tmp = VFMV_VF (RVV_ATAN2_POSINF_POSINF, vlen);
      result_tmp = __riscv_vfsgnj (result_tmp, vy, vlen);
      vw_special = __riscv_vmerge (vw_special, result_tmp, id_mask, vlen);

      id_mask = __riscv_vmand (x_posInf, y_finite, vlen);
      result_tmp = VFMV_VF (fp_posZero, vlen);
      result_tmp = __riscv_vfsgnj (result_tmp, vy, vlen);
      vw_special = __riscv_vmerge (vw_special, result_tmp, id_mask, vlen);

      id_mask = __riscv_vmand (x_negInf, y_Inf, vlen);
      result_tmp = VFMV_VF (RVV_ATAN2_NEGINF_INF, vlen);
      result_tmp = __riscv_vfsgnj (result_tmp, vy, vlen);
      vw_special = __riscv_vmerge (vw_special, result_tmp, id_mask, vlen);

      id_mask = __riscv_vmand (x_negInf, y_finite, vlen);
      result_tmp = VFMV_VF (RVV_ATAN2_NEGINF_FINITE, vlen);
      result_tmp = __riscv_vfsgnj (result_tmp, vy, vlen);
      vw_special = __riscv_vmerge (vw_special, result_tmp, id_mask, vlen);

      class_of_interest = class_finite_pos;
      VBOOL x_finite_pos;
      IDENTIFY (vclass_x, class_of_interest, x_finite_pos, vlen);
      id_mask = __riscv_vmand (x_finite_pos, y_0, vlen);
      vw_special = __riscv_vmerge (vw_special, vy, id_mask, vlen);

      class_of_interest = class_finite_neg;
      VBOOL x_finite_neg;
      IDENTIFY (vclass_x, class_of_interest, x_finite_neg, vlen);
      id_mask = __riscv_vmand (x_finite_neg, y_0, vlen);
      result_tmp = VFMV_VF (RVV_ATAN2_FINITE_NEGATIVE_ZERO, vlen);
      result_tmp = __riscv_vfsgnj (result_tmp, vy, vlen);
      vw_special = __riscv_vmerge (vw_special, result_tmp, id_mask, vlen);

      id_mask = __riscv_vmor (x_finite_pos, x_finite_neg, vlen);
      id_mask = __riscv_vmand (id_mask, y_Inf, vlen);
      result_tmp = VFMV_VF (RVV_ATAN2_Y_INF_FINITE, vlen);
      result_tmp = __riscv_vfsgnj (result_tmp, vy, vlen);
      vw_special = __riscv_vmerge (vw_special, result_tmp, id_mask, vlen);

#if RVV_MATH_IS_ATAN2PI (V_MATH_FUNC)
      vy = __riscv_vfmerge (vy, fp_posHalf, special_args, vlen);
      vx = __riscv_vfmerge (vx, fp_posHalf, special_args, vlen);
#else
      vy = __riscv_vfmerge (vy, 0x1.0p-59, special_args, vlen);
      vx = __riscv_vfmerge (vx, 0x1.0p0, special_args, vlen);
#endif
    }

  VUINT expo_y =
      __riscv_vand (__riscv_vsrl (F_AS_U (vy), MAN_LEN, vlen), 0x7FF, vlen);
  VUINT expo_x =
      __riscv_vand (__riscv_vsrl (F_AS_U (vx), MAN_LEN, vlen), 0x7FF, vlen);
  VINT exp_diff = __riscv_vsub (U_AS_I (expo_y), U_AS_I (expo_x), vlen);
  VBOOL exp_diff_large = __riscv_vmsge (exp_diff, 60, vlen);
  exp_diff_large =
      __riscv_vmor (exp_diff_large, __riscv_vmsle (exp_diff, -60, vlen), vlen);

  nb_special_args = __riscv_vcpop (exp_diff_large, vlen);
  special_args = __riscv_vmor (special_args, exp_diff_large, vlen);

  if (nb_special_args > 0)
    {
      VBOOL swap_yx = __riscv_vmsgtu (expo_y, expo_x, vlen);
      VBOOL x_neg = __riscv_vmslt (F_AS_I (vx), 0, vlen);

      VBOOL no_divide = __riscv_vmor (swap_yx, x_neg, vlen);
      no_divide = __riscv_vmand (no_divide, exp_diff_large, vlen);

      VBOOL divide = __riscv_vmnot (swap_yx, vlen);
      divide = __riscv_vmandn (divide, x_neg, vlen);
      divide = __riscv_vmand (divide, exp_diff_large, vlen);

#if RVV_MATH_IS_ATAN2PI (V_MATH_FUNC)
      VFLOAT tmp1 = __riscv_vfmul (divide, vx, 0x1.0p-55, vlen);
      VFLOAT tmp2 = __riscv_vfmul (divide, tmp1, RVV_ATAN2_PI_HI, vlen);
      VFLOAT tmp3 = __riscv_vfmsac (divide, tmp2, RVV_ATAN2_PI_HI, tmp1,
                                    vlen);
      tmp3 = __riscv_vfmacc (divide, tmp3, RVV_ATAN2_PI_LO, tmp1, vlen);
      VFLOAT R = __riscv_vfrdiv (divide, tmp2, fp_posOne, vlen);
      VFLOAT r = VFMV_VF (fp_posOne, vlen);
      r = __riscv_vfnmsac (divide, r, R, tmp2, vlen);
      r = __riscv_vfnmsac (divide, r, R, tmp3, vlen);
      r = __riscv_vfmul (divide, r, R, vlen);
      tmp1 = __riscv_vfmul (divide, vy, 0x1.0p55, vlen);
      tmp2 = __riscv_vfmul (divide, tmp1, r, vlen);
      tmp1 = __riscv_vfmadd (divide, tmp1, R, tmp2, vlen);
      tmp1 = __riscv_vfmul (divide, tmp1, 0x1.0p-110, vlen);
      tmp1 = __riscv_vfmerge (tmp1, 0x1.0p-60, no_divide, vlen);
      tmp1 = __riscv_vfsgnj (tmp1, vx, vlen);

      tmp2 = __riscv_vfsgnj (divide, tmp1, vy, vlen);
      vw_special = __riscv_vmerge (vw_special, tmp2, divide, vlen);

      VBOOL use_half = __riscv_vmand (swap_yx, exp_diff_large, vlen);
      tmp2 = VFMV_VF (RVV_ATAN2_Y_NONZERO_X_ZERO, vlen);
      tmp2 = __riscv_vfsgnj (use_half, tmp2, vy, vlen);
      vw_special = __riscv_vmerge (vw_special, tmp2, use_half, vlen);

      VBOOL use_one = __riscv_vmandn (x_neg, swap_yx, vlen);
      use_one = __riscv_vmand (use_one, exp_diff_large, vlen);
      tmp2 = VFMV_VF (RVV_ATAN2_NEGINF_FINITE, vlen);
      tmp2 = __riscv_vfsgnj (use_one, tmp2, vy, vlen);
      vw_special = __riscv_vmerge (vw_special, tmp2, use_one, vlen);
#else
      VFLOAT abs_y = __riscv_vfsgnj (vy, fp_posOne, vlen);
      VFLOAT tmp1 = __riscv_vfdiv (divide, abs_y, vx, vlen);
      tmp1 = __riscv_vfmerge (tmp1, 0x1.0p-60, no_divide, vlen);
      tmp1 = __riscv_vfsgnj (tmp1, vx, vlen);

      VFLOAT tmp2 = __riscv_vfsgnj (divide, tmp1, vy, vlen);
      vw_special = __riscv_vmerge (vw_special, tmp2, divide, vlen);

      VBOOL use_piby2 = __riscv_vmand (swap_yx, exp_diff_large, vlen);
      tmp2 = VFMV_VF (RVV_ATAN2_PIBY2_HI, vlen);
      tmp2 = __riscv_vfsgnj (use_piby2, tmp2, vy, vlen);
      vw_special = __riscv_vmerge (vw_special, tmp2, use_piby2, vlen);

      VBOOL use_pi = __riscv_vmandn (x_neg, swap_yx, vlen);
      use_pi = __riscv_vmand (use_pi, exp_diff_large, vlen);
      tmp2 = VFMV_VF (RVV_ATAN2_NEGINF_FINITE, vlen);
      tmp2 = __riscv_vfsgnj (use_pi, tmp2, vy, vlen);
      vw_special = __riscv_vmerge (vw_special, tmp2, use_pi, vlen);
#endif
#if RVV_MATH_IS_ATAN2PI (V_MATH_FUNC)
      vy = __riscv_vfmerge (vy, fp_posHalf, special_args, vlen);
      vx = __riscv_vfmerge (vx, fp_posHalf, special_args, vlen);
#else
      vy = __riscv_vfmerge (vy, 0x1.0p-59, special_args, vlen);
      vx = __riscv_vfmerge (vx, 0x1.0p0, special_args, vlen);
#endif
    }

  VBOOL core_args = __riscv_vmnot (special_args, vlen);

  VFLOAT abs_y = __riscv_vfsgnj (vy, fp_posOne, vlen);
  VFLOAT abs_x = __riscv_vfsgnj (vx, fp_posOne, vlen);
  VBOOL swap_yx = __riscv_vmfge (abs_y, abs_x, vlen);
  VFLOAT numer = abs_y;
  VFLOAT denom = abs_x;
  numer = __riscv_vmerge (numer, abs_x, swap_yx, vlen);
  denom = __riscv_vmerge (denom, abs_y, swap_yx, vlen);
  numer = __riscv_vfsgnj (numer, vx, vlen);

  VUINT exp_normalize = __riscv_vsrl (F_AS_U (denom), MAN_LEN, vlen);
  exp_normalize = __riscv_vmaxu (exp_normalize, 10, vlen);
  exp_normalize = __riscv_vminu (exp_normalize, 2036, vlen);
  exp_normalize = __riscv_vrsub (exp_normalize, 2046, vlen);
  VFLOAT scale_normalize = U_AS_F (__riscv_vsll (exp_normalize, MAN_LEN, vlen));
  numer = __riscv_vfmul (core_args, numer, scale_normalize, vlen);
  denom = __riscv_vfmul (core_args, denom, scale_normalize, vlen);

  VFLOAT z = __riscv_vfdiv (core_args, numer, denom, vlen);
  VFLOAT delta = numer;
  delta = __riscv_vfnmsac (core_args, delta, z, denom, vlen);
  VFLOAT numer_recip = VFREC7_M (core_args, numer, vlen);
  delta = __riscv_vfmul (core_args, delta, numer_recip, vlen);
  delta = __riscv_vfmul (core_args, delta, z, vlen);

  VUINT expo_61pm = __riscv_vsrl (F_AS_U (z), MAN_LEN, vlen);
  expo_61pm = __riscv_vand (expo_61pm, 0x7FF, vlen);
  expo_61pm = __riscv_vmaxu (expo_61pm, EXP_BIAS - 60, vlen);
  expo_61pm = __riscv_vrsub (expo_61pm, 2 * EXP_BIAS + 61, vlen);

  VFLOAT scale_61pm = U_AS_F (__riscv_vsll (expo_61pm, MAN_LEN, vlen));
  VINT Z = VFCVT_X_F_M
    (core_args, __riscv_vfmul (core_args, z, scale_61pm, vlen), vlen);
  VINT Delta = VFCVT_X_F_M
    (core_args, __riscv_vfmul (core_args, delta, scale_61pm, vlen), vlen);
  Z = __riscv_vsadd (Z, Delta, vlen);

  VINT V = __riscv_vsmul (Z, __riscv_vsll (Z, 1, vlen), 1, vlen);
  VINT VV = __riscv_vrsub (V, 0, vlen);

  VUINT m = __riscv_vsub (expo_61pm, EXP_BIAS + 61, vlen);
  VUINT two_m = __riscv_vsll (m, 1, vlen);
  VBOOL left_shift = __riscv_vmsltu (two_m, 3, vlen);
  VBOOL right_shift = __riscv_vmnot (left_shift, vlen);

  VINT I_tmp = __riscv_vsll (left_shift, V, __riscv_vrsub (two_m, 2, vlen),
                             vlen);
  V = __riscv_vmerge (V, I_tmp, left_shift, vlen);
  I_tmp = __riscv_vsll (left_shift, VV, __riscv_vrsub (two_m, 3, vlen), vlen);
  VV = __riscv_vmerge (VV, I_tmp, left_shift, vlen);

  /* RVV shifts mask the shift amount (e.g. 64 -> 0).  For large m, the
     normalization shifts can be >= 64 which must yield 0 instead.  */
  VUINT shift_v = __riscv_vsub (two_m, 2, vlen);
  VBOOL shift_v_too_large = __riscv_vmsgeu (shift_v, 64, vlen);
  shift_v = __riscv_vminu (shift_v, 63, vlen);
  I_tmp = __riscv_vsra (right_shift, V, shift_v, vlen);
  I_tmp = __riscv_vmerge (I_tmp, 0, shift_v_too_large, vlen);
  V = __riscv_vmerge (V, I_tmp, right_shift, vlen);

  VUINT shift_vv = __riscv_vsub (two_m, 3, vlen);
  VBOOL shift_vv_too_large = __riscv_vmsgeu (shift_vv, 64, vlen);
  shift_vv = __riscv_vminu (shift_vv, 63, vlen);
  I_tmp = __riscv_vsra (right_shift, VV, shift_vv, vlen);
  I_tmp = __riscv_vmerge (I_tmp, 0, shift_vv_too_large, vlen);
  VV = __riscv_vmerge (VV, I_tmp, right_shift, vlen);

  VINT WW = __riscv_vsll (__riscv_vsmul (V, VV, 1, vlen), 1, vlen);

  VINT P_even = PSTEPN_I (-0x56629d839b68685, WW,
                          PSTEPN_I (-0x3d2984d0a6f836a, WW,
                                   PSTEPN_I (-0x1c5e8b5228f9fe4, WW,
                                            PSTEPN_I (-0x05deca0ae3a1a5d,
                                                      -0x004efe42fda24d7, WW,
                                                      vlen),
                                            vlen),
                                   vlen),
                          vlen);

  P_even = PSTEPN_I (-0x2aaaaaaaaaaa49d3, WW,
                     PSTEPN_I (-0x12492492378aaf69, WW,
                              PSTEPN_I (-0xba2e88c805cbaf8, WW,
                                       PSTEPN_I (-0x888722719d1260a, WW,
                                                PSTEPN_I (-0x6b96ef57ce79cc3,
                                                         WW, P_even, vlen),
                                                vlen),
                                       vlen),
                              vlen),
                     vlen);

  VINT P_odd = PSTEPN_I (0x04afe3b1345b489b, WW,
                         PSTEPN_I (0x02cec355111c7439, WW,
                                  PSTEPN_I (0x00eaa9acebf3963e, WW,
                                           PSTEPN_I (0x001b053368ecfa14,
                                                     0x00006da7bb4399dd, WW,
                                                     vlen),
                                           vlen),
                                  vlen),
                         vlen);

  P_odd = PSTEPN_I (0x1999999999540349, WW,
                    PSTEPN_I (0x0e38e38bf1671f42, WW,
                             PSTEPN_I (0x09d89b293ef5f4d9, WW,
                                      PSTEPN_I (0x0786ec3df324db61, WW,
                                               PSTEPN_I (0x060b457b3c56e750,
                                                        WW, P_odd, vlen),
                                               vlen),
                                      vlen),
                             vlen),
                    vlen);

  P_odd = __riscv_vsmul (VV, P_odd, 1, vlen);
  VINT P = __riscv_vsub (P_even, P_odd, vlen);

  P = __riscv_vsmul (VV, P, 1, vlen);
  P = __riscv_vsmul (Z, P, 1, vlen);
  P = __riscv_vsub (Z, P, vlen);

  VBOOL xneg = __riscv_vmslt (F_AS_I (vx), 0, vlen);
  VBOOL xneg_or_swap = __riscv_vmor (xneg, swap_yx, vlen);
  VBOOL xneg_and_noswap = __riscv_vmandn (xneg, swap_yx, vlen);

  VINT P_tmp = __riscv_vsra (P, m, vlen);
  P = __riscv_vmerge (P, P_tmp, xneg_or_swap, vlen);

  P = RVV_ATAN2_POST_SCALE (P, vlen);

  P_tmp = __riscv_vrsub (P, RVV_ATAN2_FIXED_HALF, vlen);
  P = __riscv_vmerge (P, P_tmp, swap_yx, vlen);

  P_tmp = __riscv_vadd (P, RVV_ATAN2_FIXED_ONE, vlen);
  P = __riscv_vmerge (P, P_tmp, xneg_and_noswap, vlen);

  VUINT expo_scale = __riscv_vrsub (expo_61pm, 2 * EXP_BIAS, vlen);
  expo_scale = __riscv_vmerge (expo_scale, EXP_BIAS - 61, xneg_or_swap, vlen);
  VFLOAT scale_result = U_AS_F (__riscv_vsll (expo_scale, MAN_LEN, vlen));
  vw = VFCVT_F_X_M (core_args, P, vlen);
  vw = __riscv_vfmul (core_args, vw, scale_result, vlen);
  vw = __riscv_vfsgnj (vw, vy, vlen);
  vw = __riscv_vmerge (vw, vw_special, special_args, vlen);

  RESTORE_FRM;
  return vw;
}
