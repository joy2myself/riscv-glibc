/* Double-precision vector erf function.
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
  VFLOAT vx_orig = vx;
  VFLOAT vy, vy_special = VFMV_VF (fp_posZero, vlen);
  VBOOL special_args;

  SET_ROUNDTONEAREST;

  VUINT vclass = __riscv_vfclass (vx, vlen);
  VUINT expo_x = __riscv_vsrl (F_AS_U (vx), MAN_LEN, vlen);
  expo_x = __riscv_vand (expo_x, 0x7FF, vlen);
  IDENTIFY (vclass, class_NaN | class_Inf, special_args, vlen);
  VBOOL expo_small = __riscv_vmsltu (expo_x, EXP_BIAS - 30, vlen);
  special_args = __riscv_vmor (special_args, expo_small, vlen);

  UINT nb_special_args = __riscv_vcpop (special_args, vlen);
  if (nb_special_args > 0)
    {
      VBOOL id_mask;
      IDENTIFY (vclass, class_NaN, id_mask, vlen);
      vy_special = __riscv_vfadd (id_mask, vx, vx, vlen);
      IDENTIFY (vclass, class_posInf, id_mask, vlen);
      vy_special = __riscv_vfmerge (vy_special, fp_posOne, id_mask, vlen);
      IDENTIFY (vclass, class_negInf, id_mask, vlen);
      vy_special = __riscv_vfmerge (vy_special, fp_negOne, id_mask, vlen);

      VFLOAT vy_small = VFMV_VF (fp_posZero, vlen);
      vy_small = __riscv_vfmul (expo_small, vx, RVV_ERF_ALPHA_LO_D, vlen);
      vy_small
        = __riscv_vfmacc (expo_small, vy_small, RVV_ERF_ALPHA_HI_D, vx, vlen);
      vy_special = __riscv_vmerge (vy_special, vy_small, expo_small, vlen);

      vx = __riscv_vfmerge (vx, fp_posZero, special_args, vlen);
    }

  vx = __riscv_vfsgnj (vx, fp_posOne, vlen);
  vx = __riscv_vfmin (vx, 0x1.8p+2, vlen);

  VBOOL x_gt_t = __riscv_vmfgt (vx, RVV_ERF_T_D, vlen);

  VFLOAT xsq = __riscv_vfmul (vx, vx, vlen);
  VFLOAT r = __riscv_vmerge (xsq, vx, x_gt_t, vlen);
  VFLOAT delta_r = __riscv_vfmsub (vx, vx, xsq, vlen);
  delta_r = __riscv_vfmerge (delta_r, fp_posZero, x_gt_t, vlen);

  VFLOAT poly
    = PSTEP (RVV_ERF_Q_TAIL9_D, r,
             PSTEP (RVV_ERF_Q_TAIL10_D, RVV_ERF_Q_TAIL11_D, r, vlen), vlen);

  VFLOAT poly_right = VFMV_VF (fp_posZero, vlen);
  poly_right = __riscv_vmerge (poly_right, poly, x_gt_t, vlen);
  poly_right
    = PSTEP_ab (x_gt_t, RVV_ERF_Q_TAIL8_D, RVV_ERF_P_TAIL8_D, r, poly_right,
                vlen);
  poly_right
    = PSTEP_ab (x_gt_t, RVV_ERF_Q_TAIL7_D, RVV_ERF_P_TAIL7_D, r, poly_right,
                vlen);
  poly_right
    = PSTEP_ab (x_gt_t, RVV_ERF_Q_TAIL6_D, RVV_ERF_P_TAIL6_D, r, poly_right,
                vlen);
  poly_right
    = PSTEP_ab (x_gt_t, RVV_ERF_Q_TAIL5_D, RVV_ERF_P_TAIL5_D, r, poly_right,
                vlen);
  poly_right
    = PSTEP_ab (x_gt_t, RVV_ERF_Q_TAIL4_D, RVV_ERF_P_TAIL4_D, r, poly_right,
                vlen);

  VFLOAT r4 = __riscv_vfmul (r, r, vlen);
  r4 = __riscv_vfmul (r4, r4, vlen);

  VFLOAT poly_left = VFMV_VF (RVV_ERF_P_TAIL3_D, vlen);
  poly_left = __riscv_vfmerge (poly_left, RVV_ERF_Q_TAIL3_D, x_gt_t, vlen);
  poly_left
    = PSTEP_ab (x_gt_t, RVV_ERF_Q_TAIL2_D, RVV_ERF_P_TAIL2_D, r, poly_left,
                vlen);
  poly_left
    = PSTEP_ab (x_gt_t, RVV_ERF_Q_TAIL1_D, RVV_ERF_P_TAIL1_D, r, poly_left,
                vlen);
  poly_left
    = PSTEP_ab (x_gt_t, RVV_ERF_Q_TAIL0_D, RVV_ERF_P_TAIL0_D, r, poly_left,
                vlen);

  poly = __riscv_vfmadd (poly_right, r4, poly_left, vlen);
  VINT poly_fixed = __riscv_vfcvt_x (__riscv_vfmul (poly, RVV_ERF_POLY_SCALE_D,
                                                    vlen),
                                     vlen);

  VINT r_fixed
    = __riscv_vfcvt_x (__riscv_vfmul (r, RVV_ERF_R_SCALE_D, vlen), vlen);
  VINT delta_r_fixed
    = __riscv_vfcvt_x (__riscv_vfmul (delta_r, RVV_ERF_R_SCALE_D, vlen), vlen);
  r_fixed = __riscv_vadd (r_fixed, delta_r_fixed, vlen);

  VINT coeff
    = __riscv_vmerge (VMVI_VX (RVV_ERF_P_FIXED8_I, vlen), RVV_ERF_Q_FIXED8_I,
                      x_gt_t, vlen);
  poly_fixed = __riscv_vsll (__riscv_vsmul (r_fixed, poly_fixed, 1, vlen), 1,
                             vlen);
  poly_fixed = __riscv_vadd (poly_fixed, coeff, vlen);

  coeff
    = __riscv_vmerge (VMVI_VX (RVV_ERF_P_FIXED7_I, vlen), RVV_ERF_Q_FIXED7_I,
                      x_gt_t, vlen);
  poly_fixed = __riscv_vsll (__riscv_vsmul (r_fixed, poly_fixed, 1, vlen), 1,
                             vlen);
  poly_fixed = __riscv_vadd (poly_fixed, coeff, vlen);

  coeff
    = __riscv_vmerge (VMVI_VX (RVV_ERF_P_FIXED6_I, vlen), RVV_ERF_Q_FIXED6_I,
                      x_gt_t, vlen);
  poly_fixed = __riscv_vsll (__riscv_vsmul (r_fixed, poly_fixed, 1, vlen), 1,
                             vlen);
  poly_fixed = __riscv_vadd (poly_fixed, coeff, vlen);

  coeff
    = __riscv_vmerge (VMVI_VX (RVV_ERF_P_FIXED5_I, vlen), RVV_ERF_Q_FIXED5_I,
                      x_gt_t, vlen);
  poly_fixed = __riscv_vsll (__riscv_vsmul (r_fixed, poly_fixed, 1, vlen), 1,
                             vlen);
  poly_fixed = __riscv_vadd (poly_fixed, coeff, vlen);

  coeff
    = __riscv_vmerge (VMVI_VX (RVV_ERF_P_FIXED4_I, vlen), RVV_ERF_Q_FIXED4_I,
                      x_gt_t, vlen);
  poly_fixed = __riscv_vsll (__riscv_vsmul (r_fixed, poly_fixed, 1, vlen), 1,
                             vlen);
  poly_fixed = __riscv_vadd (poly_fixed, coeff, vlen);

  coeff
    = __riscv_vmerge (VMVI_VX (RVV_ERF_P_FIXED3_I, vlen), RVV_ERF_Q_FIXED3_I,
                      x_gt_t, vlen);
  poly_fixed = __riscv_vsll (__riscv_vsmul (r_fixed, poly_fixed, 1, vlen), 1,
                             vlen);
  poly_fixed = __riscv_vadd (poly_fixed, coeff, vlen);

  coeff
    = __riscv_vmerge (VMVI_VX (RVV_ERF_P_FIXED2_I, vlen), RVV_ERF_Q_FIXED2_I,
                      x_gt_t, vlen);
  poly_fixed = __riscv_vsll (__riscv_vsmul (r_fixed, poly_fixed, 1, vlen), 1,
                             vlen);
  poly_fixed = __riscv_vadd (poly_fixed, coeff, vlen);

  coeff
    = __riscv_vmerge (VMVI_VX (RVV_ERF_P_FIXED1_I, vlen), RVV_ERF_Q_FIXED1_I,
                      x_gt_t, vlen);
  poly_fixed = __riscv_vsll (__riscv_vsmul (r_fixed, poly_fixed, 1, vlen), 1,
                             vlen);
  poly_fixed = __riscv_vadd (poly_fixed, coeff, vlen);

  coeff
    = __riscv_vmerge (VMVI_VX (RVV_ERF_P_FIXED0_I, vlen), RVV_ERF_Q_FIXED0_I,
                      x_gt_t, vlen);
  poly_fixed = __riscv_vsll (__riscv_vsmul (r_fixed, poly_fixed, 1, vlen), 1,
                             vlen);
  poly_fixed = __riscv_vadd (poly_fixed, coeff, vlen);

  VINT poly_right_fixed = __riscv_vsll (poly_fixed, 3, vlen);
  poly_right_fixed = __riscv_vsmul (poly_right_fixed, poly_right_fixed, 1,
                                    vlen);
  poly_right_fixed = __riscv_vsmul (poly_right_fixed, poly_right_fixed, 1,
                                    vlen);
  poly_right_fixed = __riscv_vsmul (poly_right_fixed, poly_right_fixed, 1,
                                    vlen);

  VINT e = __riscv_vsra (F_AS_I (vx), MAN_LEN, vlen);
  e = __riscv_vmax (e, EXP_BIAS - 40, vlen);
  e = __riscv_vrsub (e, 2 * EXP_BIAS + 62, vlen);
  VFLOAT scale = I_AS_F (__riscv_vsll (e, MAN_LEN, vlen));
  VINT x_fixed = __riscv_vfcvt_x (__riscv_vfmul (vx, scale, vlen), vlen);
  poly_fixed = __riscv_vsmul (x_fixed, poly_fixed, 1, vlen);

  INT one = (1LL << 62);
  poly_right_fixed = __riscv_vsra (poly_right_fixed, 1, vlen);
  poly_right_fixed = __riscv_vrsub (poly_right_fixed, one, vlen);
  poly_fixed = __riscv_vmerge (poly_fixed, poly_right_fixed, x_gt_t, vlen);

  e = __riscv_vrsub (e, 2 * EXP_BIAS + 3, vlen);
  e = __riscv_vmerge (e, EXP_BIAS - 62, x_gt_t, vlen);
  scale = I_AS_F (__riscv_vsll (e, MAN_LEN, vlen));

  vy = __riscv_vfcvt_f (poly_fixed, vlen);
  vy = __riscv_vfmul (vy, scale, vlen);
  vy = __riscv_vfsgnj (vy, vx_orig, vlen);

  vy = __riscv_vmerge (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
