/* Double-precision vector sincos function.
   sincos(x) computes both sin(x) and cos(x)

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

/* This version computes both sin and cos simultaneously.  */

static void
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT vx_orig, VFLOAT *vsin_out,
                                           VFLOAT *vcos_out, size_t vlen)
{
  VFLOAT vx, vy, vy_special, vz, vz_special;
  VBOOL special_args;
  VUINT expo_x;

  SET_ROUNDTONEAREST;

  vx = __riscv_vfsgnj (vx_orig, fp_posOne, vlen);
  expo_x = __riscv_vsrl (F_AS_U (vx), MAN_LEN, vlen);

  /* Set results for input of NaN/Inf and also for |x| very small.  */
  VUINT vclass = __riscv_vfclass (vx_orig, vlen);
  VBOOL NaN_Inf;
  IDENTIFY (vclass, class_NaN | class_Inf, NaN_Inf, vlen);
  VBOOL small_x = __riscv_vmsleu (expo_x, EXP_BIAS - MAN_LEN - 5, vlen);
  special_args = __riscv_vmor (NaN_Inf, small_x, vlen);

  vy_special = VFMV_VF (fp_posZero, vlen);
  vz_special = VFMV_VF (fp_posZero, vlen);
  UINT nb_special_args = __riscv_vcpop (special_args, vlen);
  if (nb_special_args > 0)
    {
      /* Substitute Inf with sNaN.  */
      VBOOL id_mask;
      IDENTIFY (vclass, class_Inf, id_mask, vlen);

      vy_special
        = __riscv_vfmadd (small_x, vx_orig, 0x1.0p-60, vx_orig, vlen);
      vz_special = __riscv_vfadd (small_x, vx_orig, fp_posOne, vlen);

      VFLOAT vy_NaN_Inf = __riscv_vfmerge (vx_orig, fp_sNaN, id_mask, vlen);
      vy_NaN_Inf = __riscv_vfadd (NaN_Inf, vy_NaN_Inf, vy_NaN_Inf, vlen);
      vy_special = __riscv_vmerge (vy_special, vy_NaN_Inf, NaN_Inf, vlen);
      vz_special = __riscv_vmerge (vz_special, vy_NaN_Inf, NaN_Inf, vlen);

      /* Avoid spurious exceptions from the core path.  */
      vx = __riscv_vfmerge (vx, fp_posZero, special_args, vlen);
    }
  vx = __riscv_vfmerge (vx, fp_posZero, special_args, vlen);

  /* Argument reduction: reduce x to [-pi/4, pi/4].  */
  VBOOL x_large = __riscv_vmsgeu (expo_x, EXP_BIAS + 24, vlen);
  x_large = __riscv_vmandn (x_large, special_args, vlen);
  VFLOAT vx_copy = vx;
  vx = __riscv_vfmerge (vx, fp_posZero, x_large, vlen);

  VFLOAT n_flt = __riscv_vfmul (vx, PIBY2_INV, vlen);
  VINT n = __riscv_vfcvt_x (n_flt, vlen);
  n_flt = __riscv_vfcvt_f (n, vlen);
  VFLOAT r_hi = __riscv_vfnmsac (vx, PIBY2_HI, n_flt, vlen);
  VUINT expo_r = __riscv_vsrl (F_AS_U (r_hi), MAN_LEN, vlen);
  expo_r = __riscv_vand (expo_r, 0x7FF, vlen);
  VBOOL r_small = __riscv_vmsleu (expo_r, EXP_BIAS - 16, vlen);
  UINT nb_r_small = __riscv_vcpop (r_small, vlen);
  VFLOAT r = __riscv_vfnmsac (r_hi, PIBY2_MID, n_flt, vlen);
  VFLOAT r_delta = __riscv_vfsub (r_hi, r, vlen);
  r_delta = __riscv_vfnmsac (r_delta, PIBY2_MID, n_flt, vlen);
  /* At this point, r + r_delta is an accurate reduced argument PROVIDED
     |r_hi| >= 2^(-15).  */
  if (nb_r_small > 0)
    {
      VFLOAT A = __riscv_vfmul (n_flt, PIBY2_MID, vlen);
      VFLOAT a = __riscv_vfmsub (n_flt, PIBY2_MID, A, vlen);
      /* A + a is n * piby2_mid exactly.  */
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
      LARGE_ARGUMENT_REDUCTION_Piby2 (vx_copy, vlen, x_large, n_xlarge,
                                      r_xlarge, r_delta_xlarge);
      r = __riscv_vmerge (r, r_xlarge, x_large, vlen);
      r_delta = __riscv_vmerge (r_delta, r_delta_xlarge, x_large, vlen);
      n = __riscv_vmerge (n, n_xlarge, x_large, vlen);
    }

  VUINT n_lsb = __riscv_vand (I_AS_U (n), 0x1, vlen);
  VBOOL sin_use_c = __riscv_vmsne (n_lsb, 0, vlen);

  /* We compute both sin(r) and cos(r).  */
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

  vy = __riscv_vmerge (poly_sin, poly_cos, sin_use_c, vlen);
  vz = __riscv_vmerge (poly_cos, poly_sin, sin_use_c, vlen);

  n = __riscv_vsll (n, BIT_WIDTH - 2, vlen);
  vy = __riscv_vfsgnjx (vy, I_AS_F (n), vlen);
  vz = __riscv_vfsgnjx (vz, I_AS_F (n), vlen);

  n_lsb = __riscv_vsll (n_lsb, 63, vlen);
  vz = __riscv_vfsgnjx (vz, U_AS_F (n_lsb), vlen);

  vy = __riscv_vfsgnjx (vy, vx_orig, vlen);

  vy = __riscv_vmerge (vy, vy_special, special_args, vlen);
  vz = __riscv_vmerge (vz, vz_special, special_args, vlen);

  *vsin_out = vy;
  *vcos_out = vz;

  RESTORE_FRM;
}

#undef LMUL
#undef SIMDLEN
