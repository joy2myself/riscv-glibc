/* Double-precision vector cospi function.
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
  VFLOAT vx, vy, vy_special = VFMV_VF (fp_posZero, vlen);
  VBOOL special_args;
  VINT n;

  SET_ROUNDTONEAREST;

  vx = __riscv_vfsgnj (vx_orig, fp_posOne, vlen);
  VUINT expo_x = __riscv_vsrl (F_AS_U (vx), MAN_LEN, vlen);

  VUINT vclass = __riscv_vfclass (vx_orig, vlen);
  VBOOL nan_inf;
  IDENTIFY (vclass, class_NaN | class_Inf, nan_inf, vlen);
  VBOOL small_x
    = __riscv_vmsleu (expo_x, EXP_BIAS - MAN_LEN - RVV_COSPI_SMALL_EXP_OFFSET,
                      vlen);
  special_args = __riscv_vmor (nan_inf, small_x, vlen);

  UINT nb_special_args = __riscv_vcpop (special_args, vlen);
  if (nb_special_args > 0)
    {
      VBOOL id_mask;
      IDENTIFY (vclass, class_Inf, id_mask, vlen);
      vy_special = __riscv_vfadd (small_x, vx_orig, fp_posOne, vlen);
      VFLOAT vy_nan_inf = __riscv_vfmerge (vx_orig, fp_sNaN, id_mask, vlen);
      vy_nan_inf = __riscv_vfadd (nan_inf, vy_nan_inf, vy_nan_inf, vlen);
      vy_special = __riscv_vmerge (vy_special, vy_nan_inf, nan_inf, vlen);
      vx = __riscv_vfmerge (vx, fp_posZero, special_args, vlen);
    }

  VBOOL x_large
    = __riscv_vmsgeu (expo_x, EXP_BIAS + RVV_COSPI_LARGE_EXP_OFFSET, vlen);
  vx = __riscv_vfmerge (vx, fp_posZero, x_large, vlen);

  VFLOAT two_x = __riscv_vfadd (vx, vx, vlen);
  n = __riscv_vfcvt_x (two_x, vlen);
  VFLOAT n_flt = __riscv_vfcvt_f (n, vlen);
  VFLOAT rem = __riscv_vfsub (two_x, n_flt, vlen);
  VBOOL x_is_n_piby2 = __riscv_vmseq (F_AS_U (rem), 0, vlen);

  VFLOAT r = __riscv_vfmul (rem, PIBY2_HI, vlen);
  VFLOAT r_delta = __riscv_vfmsac (r, PIBY2_HI, rem, vlen);
  r_delta = __riscv_vfmacc (r_delta, PIBY2_MID, rem, vlen);

  VUINT n_lsb = __riscv_vand (I_AS_U (n), 0x1, vlen);
  VBOOL pick_c = __riscv_vmseq (n_lsb, 0, vlen);
  VBOOL exact_zero = __riscv_vmandn (x_is_n_piby2, pick_c, vlen);

  VFLOAT one = VFMV_VF (fp_posOne, vlen);
  VFLOAT rsq = __riscv_vfmul (r, r, vlen);
  VFLOAT rcube = __riscv_vfmul (rsq, r, vlen);
  VFLOAT r_to_6 = __riscv_vfmul (rcube, rcube, vlen);

  VFLOAT r_prime = __riscv_vfmul (r, RVV_COSPI_R_PRIME_SCALE, vlen);
  VFLOAT c = __riscv_vfmacc (one, r_prime, r, vlen);
  VFLOAT s_corr = __riscv_vfmul (r_delta, c, vlen);

  VFLOAT c_corr = __riscv_vfsub (one, c, vlen);
  c_corr = __riscv_vfmacc (c_corr, r, r_prime, vlen);
  c_corr = __riscv_vfnmsac (c_corr, r, r_delta, vlen);

  VFLOAT poly_right = VFMV_VF (RVV_COSPI_POLY_RIGHT_SEED_B, vlen);
  poly_right = __riscv_vfmerge (poly_right, RVV_COSPI_POLY_RIGHT_SEED_A, pick_c,
                                vlen);
  poly_right = PSTEP_ab (pick_c, RVV_COSPI_POLY_RIGHT_OUTER_A,
                         RVV_COSPI_POLY_RIGHT_OUTER_B, rsq,
                         PSTEP_ab (pick_c, RVV_COSPI_POLY_RIGHT_INNER_A,
                                   RVV_COSPI_POLY_RIGHT_INNER_B, rsq, poly_right,
                                   vlen),
                         vlen);

  VFLOAT poly_left = VFMV_VF (RVV_COSPI_POLY_LEFT_SEED_B, vlen);
  poly_left = __riscv_vfmerge (poly_left, RVV_COSPI_POLY_LEFT_SEED_A, pick_c,
                               vlen);
  poly_left = PSTEP_ab (pick_c, RVV_COSPI_POLY_LEFT_OUTER_A,
                        RVV_COSPI_POLY_LEFT_OUTER_B, rsq,
                        PSTEP_ab (pick_c, RVV_COSPI_POLY_LEFT_INNER_A,
                                  RVV_COSPI_POLY_LEFT_INNER_B, rsq, poly_left,
                                  vlen),
                        vlen);

  poly_right = __riscv_vfmadd (poly_right, r_to_6, poly_left, vlen);

  VFLOAT t = __riscv_vfmul (rsq, rsq, vlen);
  t = __riscv_vmerge (rcube, t, pick_c, vlen);

  VFLOAT a_hi = __riscv_vmerge (r, c, pick_c, vlen);
  VFLOAT a_lo = __riscv_vmerge (s_corr, c_corr, pick_c, vlen);
  vy = __riscv_vfmadd (poly_right, t, a_lo, vlen);
  vy = __riscv_vfadd (a_hi, vy, vlen);

  n = __riscv_vsll (n, BIT_WIDTH - 2, vlen);
  vy = __riscv_vfsgnjx (vy, I_AS_F (n), vlen);
  n_lsb = __riscv_vsll (n_lsb, 63, vlen);
  vy = __riscv_vfsgnjx (vy, U_AS_F (n_lsb), vlen);

  vy = __riscv_vmerge (vy, VFMV_VF (fp_posZero, vlen), exact_zero, vlen);
  vy = __riscv_vmerge (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
