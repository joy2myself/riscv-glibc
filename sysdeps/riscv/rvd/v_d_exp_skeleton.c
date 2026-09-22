/* Double-precision vector exp function.
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
  VFLOAT vy, vy_special = VFMV_VF (fp_posZero, vlen);
  VBOOL special_args;

  SET_ROUNDTONEAREST;

  VUINT vclass = __riscv_vfclass (vx, vlen);
  IDENTIFY (vclass, class_NaN | class_Inf, special_args, vlen);

  UINT nb_special_args = __riscv_vcpop (special_args, vlen);
  if (nb_special_args > 0)
    {
      VBOOL id_mask;
      IDENTIFY (vclass, class_negInf, id_mask, vlen);
      vx = __riscv_vfmerge (vx, fp_posZero, id_mask, vlen);
      vy_special = __riscv_vfadd (vx, vx, vlen);
      vx = __riscv_vfmerge (vx, fp_posZero, special_args, vlen);
    }

  vx = FCLIP (vx, RVV_EXP_X_MIN_D, RVV_EXP_X_MAX_D, vlen);

  VFLOAT flt_n = __riscv_vfmul (vx, RVV_LN2_INV_D, vlen);
  VINT n = __riscv_vfcvt_x (flt_n, vlen);
  flt_n = __riscv_vfcvt_f (n, vlen);
  VFLOAT r = __riscv_vfnmsac (vx, RVV_LN2_HI_D, flt_n, vlen);
  r = __riscv_vfnmsac (r, RVV_LN2_LO_D, flt_n, vlen);

  VFLOAT poly_right = PSTEP (
      RVV_EXP_POLY_RIGHT0_D, r,
      PSTEP (RVV_EXP_POLY_RIGHT1_D, RVV_EXP_POLY_RIGHT2_D, r, vlen), vlen);

  VFLOAT poly_mid = PSTEP (
      RVV_EXP_POLY_MID0_D, r,
      PSTEP (RVV_EXP_POLY_MID1_D, RVV_EXP_POLY_MID2_D, r, vlen), vlen);

  VFLOAT poly_left = PSTEP (
      RVV_EXP_POLY_LEFT0_D, r,
      PSTEP (RVV_EXP_POLY_LEFT1_D, RVV_EXP_POLY_LEFT2_D, r, vlen), vlen);

  VFLOAT r_sq = __riscv_vfmul (r, r, vlen);
  VFLOAT r_cube = __riscv_vfmul (r_sq, r, vlen);

  VFLOAT poly = __riscv_vfmadd (poly_right, r_cube, poly_mid, vlen);
  poly = __riscv_vfmadd (poly, r_cube, poly_left, vlen);
  poly = PSTEP (RVV_EXP_POLY_FINAL_D, r, poly, vlen);

  r = __riscv_vfmacc (r, r_sq, poly, vlen);
  vy = __riscv_vfadd (r, fp_posOne, vlen);

  FAST_LDEXP (vy, n, vlen);

  vy = __riscv_vmerge (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
