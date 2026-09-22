/* Double-precision vector tan function.
   tan(x) computes the tangent of x

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
/* This version reduces argument to [-pi/4, pi/4] and computes tan(r) */

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
  VBOOL numer_pick_c = __riscv_vmsne (n_lsb, 0, vlen);
  VBOOL denom_pick_c = __riscv_vmnot (numer_pick_c, vlen);

  /* Compute high-precision sin and cos for tan = sin/cos */
  VFLOAT rsq = __riscv_vfmul (r, r, vlen);

  UINT mask_r_head = 1;
  mask_r_head = ~((mask_r_head << 36) - 1);
  VFLOAT r_head = U_AS_F (__riscv_vand (F_AS_U (r), mask_r_head, vlen));
  VFLOAT r_tail = __riscv_vfsub (r, r_head, vlen);

  UINT exp_m1 = 1;
  exp_m1 = (exp_m1 << 52);
  VFLOAT r_prime = U_AS_F (__riscv_vsub (F_AS_U (r), exp_m1, vlen));

  VFLOAT ONE = VFMV_VF (fp_posOne, vlen);
  VFLOAT cos_hi = __riscv_vfnmsac (ONE, r, r_prime, vlen);
  VFLOAT cos_corr = __riscv_vfsub (ONE, cos_hi, vlen);
  cos_corr = __riscv_vfnmsac (cos_corr, r, r_prime, vlen);
  cos_corr = __riscv_vfnmsac (cos_corr, r_delta, r, vlen);

  double coeff = -0x1.5555555555555p-3;
  VFLOAT r_head_cube = __riscv_vfmul (r_head, r_head, vlen);
  r_head_cube = __riscv_vfmul (r_head_cube, r_head, vlen);
  VFLOAT sin_hi = __riscv_vfmadd (r_head_cube, coeff, r, vlen);
  VFLOAT sin_corr = __riscv_vfsub (r, sin_hi, vlen);
  sin_corr = __riscv_vfmacc (sin_corr, coeff, r_head_cube, vlen);
  VFLOAT tmp = __riscv_vfmadd (r_head, r_head, rsq, vlen);
  VFLOAT tmp2 = __riscv_vfmul (r_tail, coeff, vlen);
  tmp = __riscv_vfmacc (tmp, r_head, r, vlen);
  sin_corr = __riscv_vfmacc (sin_corr, tmp, tmp2, vlen);
  sin_corr = __riscv_vfmacc (sin_corr, r_delta, cos_hi, vlen);

  VFLOAT poly_s =
      PSTEP (0x1.1111111111069p-7, rsq,
            PSTEP (-0x1.a01a019ffe527p-13, rsq,
                  PSTEP (0x1.71de3a33a62c6p-19, rsq,
                        PSTEP (-0x1.ae642c52fc493p-26, rsq,
                              PSTEP (0x1.6109be886e15cp-33,
                                    -0x1.9ffe1dd295e78p-41, rsq, vlen),
                              vlen),
                        vlen),
                  vlen),
            vlen);

  VFLOAT poly_c =
      PSTEP (0x1.5555555555546p-5, rsq,
            PSTEP (-0x1.6c16c16c1450cp-10, rsq,
                  PSTEP (0x1.a01a019b77545p-16, rsq,
                        PSTEP (-0x1.27e4f72551e3dp-22, rsq,
                              PSTEP (0x1.1ee950032f74cp-29,
                                    -0x1.8f5dd75850673p-37, rsq, vlen),
                              vlen),
                        vlen),
                  vlen),
            vlen);

  VFLOAT r_to_4 = __riscv_vfmul (rsq, rsq, vlen);
  VFLOAT r_to_5 = __riscv_vfmul (r_to_4, r, vlen);

  poly_c = __riscv_vfmadd (poly_c, r_to_4, cos_corr, vlen);
  poly_s = __riscv_vfmadd (poly_s, r_to_5, sin_corr, vlen);

  VFLOAT S, s, C, c;
  FAST2SUM (sin_hi, poly_s, S, s, vlen);
  FAST2SUM (cos_hi, poly_c, C, c, vlen);

  VFLOAT numer_hi, numer_lo, denom_hi, denom_lo;
  numer_hi = S;
  numer_hi = __riscv_vmerge (numer_hi, C, numer_pick_c, vlen);
  numer_lo = s;
  numer_lo = __riscv_vmerge (numer_lo, c, numer_pick_c, vlen);

  denom_hi = S;
  denom_hi = __riscv_vmerge (denom_hi, C, denom_pick_c, vlen);
  denom_lo = s;
  denom_lo = __riscv_vmerge (denom_lo, c, denom_pick_c, vlen);

  DIV_N2D2 (numer_hi, numer_lo, denom_hi, denom_lo, vy, vlen);

  /* Apply correct sign */
  n = __riscv_vsll (n, BIT_WIDTH - 1, vlen);
  vy = __riscv_vfsgnjx (vy, I_AS_F (n), vlen);
  vy = __riscv_vfsgnjx (vy, vx_orig, vlen);

  vy = __riscv_vmerge (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
