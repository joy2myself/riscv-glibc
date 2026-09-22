/* Double-precision vector sin function.
   sin(x) computes the sine of x

   Based on RISC-V Vector Library for Math (rvvlm) implementation.

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
/* This version reduces argument to [-pi/4, pi/4] and computes sin(r) */

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT vx_orig = x;
  VFLOAT vx, vy, vy_special = VFMV_VF (fp_posZero, vlen);
  VBOOL special_args;
  VUINT expo_x;

  SET_ROUNDTONEAREST;

  vx = __riscv_vfsgnj (vx_orig, fp_posOne, vlen);
  expo_x = __riscv_vsrl (F_AS_U (vx), MAN_LEN, vlen);

  /* Set results for input of NaN and Inf and also for |x| very small */
  EXCEPTION_HANDLING_TRIG(vx_orig, expo_x, special_args, vy_special, vlen);
  vx = __riscv_vfmerge (vx, fp_posZero, special_args, vlen);

  VBOOL x_large = __riscv_vmsgeu (expo_x, EXP_BIAS + 24, vlen);
  x_large = __riscv_vmandn (x_large, special_args, vlen);
  VFLOAT vx_copy = vx;
  vx = __riscv_vfmerge (vx, fp_posZero, x_large, vlen);

  /* Argument reduction: reduce x to [-pi/4, pi/4] */
  VFLOAT n_flt = __riscv_vfmul (vx, PIBY2_INV, vlen);
  VINT n = __riscv_vfcvt_x (n_flt, vlen);
  n_flt = __riscv_vfcvt_f (n, vlen);
  VFLOAT r_hi = __riscv_vfnmsac (vx, PIBY2_HI, n_flt, vlen);
  VUINT expo_r = __riscv_vsrl (F_AS_U (r_hi), MAN_LEN, vlen);
  expo_r = __riscv_vand (expo_r, 0x7FF, vlen);
  VBOOL r_small = __riscv_vmsleu (expo_r, EXP_BIAS - 16, vlen);
  VFLOAT r = __riscv_vfnmsac (r_hi, PIBY2_MID, n_flt, vlen);
  VFLOAT r_delta = __riscv_vfsub (r_hi, r, vlen);
  r_delta = __riscv_vfnmsac (r_delta, PIBY2_MID, n_flt, vlen);

  if (__riscv_vcpop (r_small, vlen) > 0)
    {
      VFLOAT A = __riscv_vfmul (n_flt, PIBY2_MID, vlen);
      VFLOAT a = __riscv_vfmsub (n_flt, PIBY2_MID, A, vlen);
      VFLOAT S = __riscv_vfsub (r_hi, A, vlen);
      VFLOAT s = __riscv_vfsub (r_hi, S, vlen);
      s = __riscv_vfsub (s, A, vlen);
      s = __riscv_vfsub (s, a, vlen);
      s = __riscv_vfnmsac (s, PIBY2_LO, n_flt, vlen);
      r = __riscv_vmerge (r, S, r_small, vlen);
      r_delta = __riscv_vmerge (r_delta, s, r_small, vlen);
    }

  if (__riscv_vcpop (x_large, vlen) > 0)
    {
      VFLOAT r_xlarge, r_delta_xlarge;
      VINT n_xlarge;
      LARGE_ARGUMENT_REDUCTION_Piby2 (vx_copy, vlen, x_large, n_xlarge, r_xlarge,
                                      r_delta_xlarge);
      r = __riscv_vmerge (r, r_xlarge, x_large, vlen);
      r_delta = __riscv_vmerge (r_delta, r_delta_xlarge, x_large, vlen);
      n = __riscv_vmerge (n, n_xlarge, x_large, vlen);
    }

  VUINT n_lsb = __riscv_vand (I_AS_U (n), 0x1, vlen);
  VBOOL pick_c = __riscv_vmsne (n_lsb, 0, vlen);

  /* Compute sin(r) or cos(r) based on quadrant */
  VFLOAT rsq = __riscv_vfmul (r, r, vlen);
  VFLOAT rcube = __riscv_vfmul (rsq, r, vlen);
  VFLOAT r_to_6 = __riscv_vfmul (rcube, rcube, vlen);

  VFLOAT r_prime = __riscv_vfmul (r, -0x1.0p-1, vlen);
  VFLOAT One = VFMV_VF (fp_posOne, vlen);
  VFLOAT C = __riscv_vfmacc (One, r_prime, r, vlen);
  VFLOAT s_corr = __riscv_vfmul (r_delta, C, vlen);

  VFLOAT c_corr = __riscv_vfsub (One, C, vlen);
  c_corr = __riscv_vfmacc (c_corr, r, r_prime, vlen);
  c_corr = __riscv_vfnmsac (c_corr, r, r_delta, vlen);

  VFLOAT poly_right = VFMV_VF (0x1.5d8b5ae12066ap-33, vlen);
  poly_right = __riscv_vfmerge (poly_right, -0x1.8f5dd75850673p-37, pick_c, vlen);
  poly_right =
      PSTEP_ab (pick_c, -0x1.27e4f72551e3dp-22, 0x1.71de35553ddb6p-19, rsq,
                PSTEP_ab (pick_c, 0x1.1ee950032f74cp-29, -0x1.ae5e4b94836f8p-26,
                          rsq, poly_right, vlen),
                vlen);

  VFLOAT poly_left = VFMV_VF (-0x1.a01a019be932ap-13, vlen);
  poly_left = __riscv_vfmerge (poly_left, 0x1.a01a019b77545p-16, pick_c, vlen);
  poly_left =
      PSTEP_ab (pick_c, 0x1.5555555555546p-5, -0x1.5555555555548p-3, rsq,
                PSTEP_ab (pick_c, -0x1.6c16c16c1450cp-10, 0x1.111111110f730p-7,
                          rsq, poly_left, vlen),
                vlen);

  poly_right = __riscv_vfmadd (poly_right, r_to_6, poly_left, vlen);

  VFLOAT t = __riscv_vfmul (rsq, rsq, vlen);
  t = __riscv_vmerge (rcube, t, pick_c, vlen);

  VFLOAT A = __riscv_vmerge (r, C, pick_c, vlen);
  VFLOAT a = __riscv_vmerge (s_corr, c_corr, pick_c, vlen);
  vy = __riscv_vfmadd (poly_right, t, a, vlen);
  vy = __riscv_vfadd (A, vy, vlen);

  /* Apply sign based on quadrant */
  n = __riscv_vsll (n, BIT_WIDTH - 2, vlen);
  vy = __riscv_vfsgnjx (vy, I_AS_F (n), vlen);

  /* Apply sign from original input */
  vy = __riscv_vfsgnjx (vy, vx_orig, vlen);

  vy = __riscv_vmerge (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
