/* Double-precision vector expm1 function.
   expm1(x) computes e^x - 1

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

/* This version uses EPsim version of exp to compute expm1. */

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT vx = x;
  VFLOAT vy, vy_special;
  VBOOL special_args;

  SET_ROUNDTONEAREST;

  VUINT vclass = __riscv_vfclass (vx, vlen);
  IDENTIFY (vclass, class_NaN | class_Inf, special_args, vlen);
  vy_special = vx;

  UINT nb_special_args = __riscv_vcpop (special_args, vlen);
  if (nb_special_args > 0)
    {
      /* Substitute -Inf with -1 */
      VBOOL id_mask;
      IDENTIFY (vclass, class_negInf, id_mask, vlen);
      vx = __riscv_vfmerge (vx, fp_negOne, id_mask, vlen);
      vy_special = __riscv_vfmul (special_args, vx, fp_posOne, vlen);
      vx = __riscv_vfmerge (vx, fp_posZero, special_args, vlen);
    }

  /* Clip */
  vx = FCLIP (vx, RVV_EXPM1_X_MIN_D, RVV_EXP_X_MAX_D, vlen);

  /* Argument reduction */
  VFLOAT flt_n = __riscv_vfmul (vx, RVV_LN2_INV_D, vlen);
  VINT n = __riscv_vfcvt_x (flt_n, vlen);
  flt_n = __riscv_vfcvt_f (n, vlen);

  VFLOAT r_tmp = __riscv_vfnmsac (vx, RVV_LN2_HI_D, flt_n, vlen);
  VFLOAT r = __riscv_vfnmsub (flt_n, RVV_LN2_LO_D, r_tmp, vlen);
  VFLOAT r_lo = __riscv_vfsub (r_tmp, r, vlen);
  r_lo = __riscv_vfnmsac (r_lo, RVV_LN2_LO_D, flt_n, vlen);
  /* r is the reduced argument in working precision; but r + r_lo is extra
     precise exp(r+r_lo) is 1 + (r+r_lo) + (r+r_lo)^2/2 + r^3 * polynomial(r)
     1 + r + r^2/2 + (r_lo + r * r_lo) + r^3 * polynomial(r) */
  r_lo = __riscv_vfmacc (r_lo, r, r_lo, vlen);
  /* 1 + r + r^2/2 + r_lo + r^3 * polynomial(r) */

  /* Compute P_head + P_tail as r + r^2/2 accurately */
  VFLOAT r_prime = __riscv_vfmul (r, RVV_HALF_D, vlen);
  VFLOAT P_head = __riscv_vfmadd (r, r_prime, r, vlen);
  VFLOAT P_tail = __riscv_vfsub (r, P_head, vlen);
  P_tail = __riscv_vfmacc (P_tail, r, r_prime, vlen);

  /* Polynomial computation, we have a degree 11 polynomial */
  VFLOAT poly_right = PSTEP (
      RVV_EXPM1_POLY_RIGHT0_D, r,
      PSTEP (RVV_EXPM1_POLY_RIGHT1_D, r,
            PSTEP (RVV_EXPM1_POLY_RIGHT2_D, RVV_EXPM1_POLY_RIGHT3_D, r, vlen),
            vlen),
      vlen);

  VFLOAT poly_mid = PSTEP (
      RVV_EXPM1_POLY_MID0_D, r,
      PSTEP (RVV_EXPM1_POLY_MID1_D, RVV_EXPM1_POLY_MID2_D, r, vlen), vlen);

  VFLOAT poly_left =
      PSTEP (RVV_EXPM1_POLY_LEFT0_D, r,
            PSTEP (RVV_EXPM1_POLY_LEFT1_D, RVV_EXPM1_POLY_LEFT2_D, r, vlen),
            vlen);

  VFLOAT r_sq = __riscv_vfmul (r, r, vlen);
  VFLOAT r_cube = __riscv_vfmul (r_sq, r, vlen);

  VFLOAT poly = __riscv_vfmadd (poly_right, r_cube, poly_mid, vlen);
  poly = __riscv_vfmadd (poly, r_cube, poly_left, vlen);
  /* At this point, exp(r) is 1 + P_head + P_tail + r_lo + r^3 * poly
     expm1(x) = 2^n ( 1 - 2^(-n) + P_head + P_tail + r_lo + r^3 * poly ) */

  /* Compute 1 - 2^(-n) accurately */
  VFLOAT One = VFMV_VF (fp_posOne, vlen);
  VINT n_clip = __riscv_vmin (n, 64, vlen);
  VBOOL n_le53 = __riscv_vmsle (n_clip, 53, vlen);
  VBOOL n_le_neg54 = __riscv_vmsle (n_clip, -54, vlen);
  VINT I_tail = __riscv_vrsub (n_clip, -(EXP_BIAS + 2), vlen);
  VFLOAT u = I_AS_F (__riscv_vsll (I_tail, MAN_LEN, vlen)); /* u = -2^(-n_clip) */
  I_tail = __riscv_vmerge (I_tail, 0, n_le53, vlen);
  I_tail = __riscv_vmerge (I_tail, EXP_BIAS, n_le_neg54, vlen);
  VFLOAT a = I_AS_F (__riscv_vsll (I_tail, MAN_LEN, vlen));
  VFLOAT A = __riscv_vfadd (One, u, vlen);

  /* Compute A + a + P_head + P_tail + r_lo + r^3 * poly */
  P_tail = __riscv_vfadd (P_tail, a, vlen);
  P_tail = __riscv_vfadd (P_tail, r_lo, vlen);
  poly = __riscv_vfmadd (poly, r_cube, P_tail, vlen);
  P_head = __riscv_vfadd (P_head, poly, vlen);

  vy = __riscv_vfadd (P_head, A, vlen);
  /* vy is now exp(r) - 1 accurately. */

  /* Need to compute 2^n * exp(r). */
  FAST_LDEXP (vy, n, vlen);

  /* Special lanes were sanitized in VX above.  Use that value for the
     tiny-input predicate so discarded NaN/Inf lanes cannot raise an
     unrelated comparison exception; the final merge still restores the
     original special result from VY_SPECIAL.  */
  VFLOAT abs_x = __riscv_vfabs (vx, vlen);
  VBOOL tiny_x = __riscv_vmflt (abs_x, RVV_EXPM1_TINY_BOUND_D, vlen);
  vy = VMERGE_VVM_F (vy, x, tiny_x, vlen);
  vy = __riscv_vmerge (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
