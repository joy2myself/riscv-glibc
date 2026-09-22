/* Double-precision vector acosh function.
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
  VFLOAT vy, vy_special;
  VBOOL special_args;

  SET_ROUNDTONEAREST;

  /* Exception handling: Inf, NaN, and input <= 1.0.  */
  VFLOAT vxm1 = __riscv_vfsub (vx, fp_posOne, vlen);
  VUINT vclass = __riscv_vfclass (vxm1, vlen);

  VBOOL is_negative;
  IDENTIFY (vclass, class_negative, is_negative, vlen);
  vxm1 = __riscv_vfmerge (vxm1, fp_sNaN, is_negative, vlen);

  IDENTIFY (vclass, class_NaN | class_Inf | class_negative | class_Zero,
            special_args, vlen);
  UINT nb_special_args = __riscv_vcpop (special_args, vlen);

  vy_special = vx;
  if (nb_special_args > 0)
    {
      vy_special = __riscv_vfmul (special_args, vxm1, vxm1, vlen);
      vx = __riscv_vfmerge (vx, fp_posOne, special_args, vlen);
    }

  /* Scale x down by 2^(-550) if x >= 2^500 to avoid overflow.  */
  VUINT expo_x = __riscv_vsrl (F_AS_U (vx), MAN_LEN, vlen);
  VBOOL x_large = __riscv_vmsgeu (expo_x, EXP_BIAS + 500, vlen);
  VINT n = VMVI_VX (0, vlen);
  n = __riscv_vmerge (n, 550, x_large, vlen);
  VFLOAT u = VFMV_VF (-0x1.0p0, vlen);
  u = __riscv_vfmerge (u, fp_posZero, x_large, vlen);
  vx = I_AS_F (__riscv_vsub (F_AS_I (vx), __riscv_vsll (n, MAN_LEN, vlen),
                             vlen));

  /* Compute x*x + u with extra precision.  */
  VFLOAT p_hi, p_lo;
  PROD_X1Y1 (vx, vx, p_hi, p_lo, vlen);
  VFLOAT tmp1, tmp2;
  FAST2SUM (p_hi, u, tmp1, tmp2, vlen);
  tmp2 = __riscv_vfadd (tmp2, p_lo, vlen);
  VFLOAT a_hi, a_lo;
  FAST2SUM (tmp1, tmp2, a_hi, a_lo, vlen);

  /* Compute sqrt(x*x + u) with extra precision.  */
  VFLOAT a_sum = __riscv_vfadd (a_hi, a_lo, vlen);
  VBOOL x_eq_0 = __riscv_vmfeq (a_sum, fp_posZero, vlen);
  a_sum = __riscv_vfmerge (a_sum, fp_posOne, x_eq_0, vlen);
  VFLOAT b_hi = __riscv_vfsqrt (a_sum, vlen);
  VFLOAT b_lo = __riscv_vfnmsub (b_hi, b_hi, a_hi, vlen);
  b_lo = __riscv_vfadd (b_lo, a_lo, vlen);
  b_lo = __riscv_vfmul (b_lo, __riscv_vfrec7 (a_sum, vlen), vlen);
  b_lo = __riscv_vfmul (b_lo, 0x1.0p-1, vlen);
  b_lo = __riscv_vfmul (b_lo, b_hi, vlen);
  b_hi = __riscv_vfmerge (b_hi, fp_posZero, x_eq_0, vlen);
  b_lo = __riscv_vfmerge (b_lo, fp_posZero, x_eq_0, vlen);

  /* x + sqrt(x*x + u) with extra precision.  */
  VFLOAT s_hi, s_lo;
  FAST2SUM (vx, b_hi, s_hi, s_lo, vlen);
  s_lo = __riscv_vfadd (s_lo, b_lo, vlen);

  /* Scale S and s to keep the log core near 1.  */
  VINT expo_s = __riscv_vsra (F_AS_I (s_hi), MAN_LEN - 8, vlen);
  expo_s = __riscv_vadd (expo_s, 0x96, vlen);
  expo_s = __riscv_vsra (expo_s, 8, vlen);
  VINT n_adjust = __riscv_vsub (expo_s, EXP_BIAS, vlen);
  n = __riscv_vadd (n, n_adjust, vlen);
  VINT scale_exp = __riscv_vsll (
      __riscv_vrsub (expo_s, 2 * EXP_BIAS, vlen), MAN_LEN, vlen);
  s_hi = I_AS_F (__riscv_vsub (F_AS_I (s_hi),
                               __riscv_vsll (n_adjust, MAN_LEN, vlen), vlen));
  s_lo = __riscv_vfmul (s_lo, I_AS_F (scale_exp), vlen);

  /* log(y) = 2 * atanh ((y - 1) / (y + 1)), y = s_hi + s_lo.  */
  VFLOAT s_tmp = __riscv_vfsub (s_hi, fp_posOne, vlen);
  VFLOAT numer, delta_numer;
  FAST2SUM (s_tmp, s_lo, numer, delta_numer, vlen);
  numer = __riscv_vfadd (numer, numer, vlen);
  delta_numer = __riscv_vfadd (delta_numer, delta_numer, vlen);

  VFLOAT denom, delta_denom;
  FAST2SUM (VFMV_VF (fp_posOne, vlen), s_hi, denom, delta_denom, vlen);
  delta_denom = __riscv_vfadd (delta_denom, s_lo, vlen);

  VFLOAT r_hi, r_lo, r;
  DIV2_N2D2 (numer, delta_numer, denom, delta_denom, r_hi, r_lo, vlen);
  r = __riscv_vfadd (r_hi, r_lo, vlen);

  VFLOAT rsq = __riscv_vfmul (r, r, vlen);
  VFLOAT rcube = __riscv_vfmul (rsq, r, vlen);
  VFLOAT r6 = __riscv_vfmul (rcube, rcube, vlen);

  VFLOAT poly_right = PSTEP (
      0x1.c71c543983a27p-12, rsq,
      PSTEP (0x1.7465c27ee47d0p-14, rsq,
             PSTEP (0x1.39af2e90a6554p-16, 0x1.2e74f2255e096p-18, rsq, vlen),
             vlen),
      vlen);
  VFLOAT poly_left = PSTEP (
      0x1.555555555558cp-4, rsq,
      PSTEP (0x1.9999999982550p-7, 0x1.2492493f7cc71p-9, rsq, vlen), vlen);
  VFLOAT poly = __riscv_vfmadd (poly_right, r6, poly_left, vlen);
  poly = __riscv_vfmadd (poly, rcube, r_lo, vlen);

  /* Reconstruction: n * log(2) + r_hi + poly.  */
  VFLOAT n_flt = __riscv_vfcvt_f (n, vlen);
  a_hi = __riscv_vfmul (n_flt, RVV_LOG2_HI_D, vlen);
  FAST2SUM (a_hi, r_hi, s_hi, s_lo, vlen);
  s_lo = __riscv_vfmacc (s_lo, RVV_LOG2_LO_D, n_flt, vlen);
  s_lo = __riscv_vfadd (s_lo, poly, vlen);

  vy = __riscv_vfadd (s_hi, s_lo, vlen);
  vy = __riscv_vmerge (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
