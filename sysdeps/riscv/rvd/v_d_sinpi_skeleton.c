/* Double-precision vector sinpi function.
   sinpi(x) computes sin(pi * x)

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

/* This version uses trigpi reduction:
   N = rint(2|x|); rem = 2|x| - N with |rem| <= 1/2.
   Then pi*x = N*(pi/2) + rem*(pi/2), so reduced |rem*pi/2| <= pi/4.
   Evaluate sin/cos polynomials on the reduced argument.  */

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT vx_orig = x;
  VFLOAT vx, vy, vy_special;
  VBOOL special_args;
  VUINT expo_x;

  SET_ROUNDTONEAREST;

  vx = __riscv_vfsgnj (vx_orig, fp_posOne, vlen);
  expo_x = __riscv_vsrl (F_AS_U (vx), MAN_LEN, vlen);

  /* Special cases: NaN/Inf and |x| very small.  */
  VUINT vclass = __riscv_vfclass (vx_orig, vlen);
  VBOOL NaN_Inf;
  IDENTIFY (vclass, class_NaN | class_Inf, NaN_Inf, vlen);
  VBOOL small_x = __riscv_vmsleu (expo_x, EXP_BIAS - MAN_LEN - 5, vlen);
  special_args = __riscv_vmor (NaN_Inf, small_x, vlen);

  vy_special = VFMV_VF (fp_posZero, vlen);
  UINT nb_special_args = __riscv_vcpop (special_args, vlen);
  if (nb_special_args > 0)
    {
      /* Substitute Inf with sNaN and propagate NaNs via x+x.  */
      VBOOL id_mask;
      IDENTIFY (vclass, class_Inf, id_mask, vlen);

      /* For very small |x|, sin(pi*x) ~= pi*x.  */
      VFLOAT sin_near_zero
        = __riscv_vfmul (small_x, vx_orig, RVV_PI_MID_D, vlen);
      sin_near_zero = __riscv_vfmadd (small_x, vx_orig, RVV_PI_HI_D,
                                      sin_near_zero, vlen);
      vy_special = sin_near_zero;

      VFLOAT vy_NaN_Inf = __riscv_vfmerge (vx_orig, fp_sNaN, id_mask, vlen);
      vy_NaN_Inf = __riscv_vfadd (NaN_Inf, vy_NaN_Inf, vy_NaN_Inf, vlen);
      vy_special = __riscv_vmerge (vy_special, vy_NaN_Inf, NaN_Inf, vlen);

      /* Avoid spurious exceptions from the core path.  */
      vx = __riscv_vfmerge (vx, fp_posZero, special_args, vlen);
    }

  /* For |x| >= 2^53, the input is an even integer. Hence sinpi(x) is a signed
     zero. Mask these lanes out from argument reduction.  */
  VBOOL x_large = __riscv_vmsgeu (expo_x, EXP_BIAS + 53, vlen);
  vx = __riscv_vfmerge (vx, fp_posZero, x_large, vlen);

  /* Argument reduction.  */
  VFLOAT two_x = __riscv_vfadd (vx, vx, vlen);
  VINT n = __riscv_vfcvt_x (two_x, vlen);
  VFLOAT n_flt = __riscv_vfcvt_f (n, vlen);
  VFLOAT rem = __riscv_vfsub (two_x, n_flt, vlen);
  VBOOL x_is_n_piby2 = __riscv_vmseq (F_AS_U (rem), 0, vlen);

  /* rem * (pi/2) as r + r_delta.  */
  VFLOAT r = __riscv_vfmul (rem, PIBY2_HI, vlen);
  VFLOAT r_delta = __riscv_vfmsac (r, PIBY2_HI, rem, vlen);
  r_delta = __riscv_vfmacc (r_delta, PIBY2_MID, rem, vlen);

  VUINT n_lsb = __riscv_vand (I_AS_U (n), 0x1, vlen);
  VBOOL sin_use_c = __riscv_vmsne (n_lsb, 0, vlen);

  /* Exact zeros for integers: rem == 0 and N even.  */
  VBOOL sinX_zero = __riscv_vmandn (x_is_n_piby2, sin_use_c, vlen);

  /* Compute both sin(r+r_delta) and cos(r+r_delta) for |r| <= pi/4.  */
  VFLOAT One = VFMV_VF (fp_posOne, vlen);
  VFLOAT rsq = __riscv_vfmul (r, r, vlen);
  VFLOAT rcube = __riscv_vfmul (rsq, r, vlen);

  VFLOAT r_prime = __riscv_vfmul (r, -0x1.0p-1, vlen);
  VFLOAT C = __riscv_vfmacc (One, r_prime, r, vlen);
  VFLOAT s_corr = __riscv_vfmul (r_delta, C, vlen);

  VFLOAT c_corr = __riscv_vfsub (One, C, vlen);
  c_corr = __riscv_vfmacc (c_corr, r, r_prime, vlen);
  c_corr = __riscv_vfnmsac (c_corr, r, r_delta, vlen);

  VFLOAT poly_sin
      = PSTEP (-0x1.5555555555548p-3, rsq,
               PSTEP (0x1.111111110f730p-7, rsq,
                      PSTEP (-0x1.a01a019be932ap-13, rsq,
                             PSTEP (0x1.71de35553ddb6p-19, rsq,
                                    PSTEP (-0x1.ae5e4b94836f8p-26,
                                           0x1.5d8b5ae12066ap-33, rsq, vlen),
                                    vlen),
                             vlen),
                      vlen),
               vlen);

  VFLOAT poly_cos
      = PSTEP (0x1.5555555555546p-5, rsq,
               PSTEP (-0x1.6c16c16c1450cp-10, rsq,
                      PSTEP (0x1.a01a019b77545p-16, rsq,
                             PSTEP (-0x1.27e4f72551e3dp-22, rsq,
                                    PSTEP (0x1.1ee950032f74cp-29,
                                           -0x1.8f5dd75850673p-37, rsq, vlen),
                                    vlen),
                             vlen),
                      vlen),
               vlen);

  VFLOAT r_4 = __riscv_vfmul (rsq, rsq, vlen);

  poly_sin = __riscv_vfmadd (poly_sin, rcube, s_corr, vlen);
  poly_sin = __riscv_vfadd (r, poly_sin, vlen);

  poly_cos = __riscv_vfmadd (poly_cos, r_4, c_corr, vlen);
  poly_cos = __riscv_vfadd (C, poly_cos, vlen);

  /* sin(pi*x) from (sin, cos) based on quadrant N mod 4.  */
  vy = __riscv_vmerge (poly_sin, poly_cos, sin_use_c, vlen);

  /* Sign adjustment from quadrant N mod 4 (flip when bit1 is set).  */
  n = __riscv_vsll (n, BIT_WIDTH - 2, vlen);
  vy = __riscv_vfsgnjx (vy, I_AS_F (n), vlen);

  /* Force exact zeros at integer points.  */
  vy = __riscv_vmerge (vy, VFMV_VF (fp_posZero, vlen), sinX_zero, vlen);

  /* sin(pi*x) is odd.  */
  vy = __riscv_vfsgnjx (vy, vx_orig, vlen);

  vy = __riscv_vmerge (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
