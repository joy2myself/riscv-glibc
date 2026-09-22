/* Double-precision vector log2 function.
   log2(x) computes the base-2 logarithm of x

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

/* This version uses a 128-entry lookup table for high precision */

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  const struct v_log_data *d = ptr_barrier (&__v_log_data);
  VFLOAT vx = x;
  VFLOAT vy, vy_special = vx;
  VBOOL special_args;
  VINT n_adjust;

  SET_ROUNDTONEAREST;

  /* Handle special cases: NaN, Inf, negative, zero, and denormals */
  VUINT vclass = __riscv_vfclass (vx, vlen);
  /* special handling except positive normal number */
  IDENTIFY (vclass, 0x3BF, special_args, vlen);

  UINT nb_special_args = __riscv_vcpop (special_args, vlen);
  n_adjust = VMVI_VX (0, vlen);
  if (nb_special_args > 0)
    {
      VBOOL id_mask;
      /* substitute negative arguments with sNaN */
      IDENTIFY (vclass, class_negative, id_mask, vlen);
      vx = __riscv_vfmerge (vx, fp_sNaN, id_mask, vlen);
      /* substitute +0 argument with -0 */
      IDENTIFY (vclass, class_posZero, id_mask, vlen);
      vx = __riscv_vfmerge (vx, fp_negZero, id_mask, vlen);
      /* eliminate positive denorm input from special_args */
      IDENTIFY (vclass, 0x39F, special_args, vlen);
      /* for narrowed set of special arguments, compute vx+vfrec7(vx) */
      vy_special = __riscv_vfrec7 (special_args, vx, vlen);
      vy_special = __riscv_vfadd (special_args, vy_special, vx, vlen);
      vx = __riscv_vfmerge (vx, fp_posOne, special_args, vlen);
      /* scale up input for positive denormals */
      IDENTIFY (vclass, class_posDenorm, id_mask, vlen);
      n_adjust = __riscv_vmerge (n_adjust, 64, id_mask, vlen);
      VFLOAT vx_normalized = __riscv_vfmul (id_mask, vx, 0x1.0p64, vlen);
      vx = __riscv_vmerge (vx, vx_normalized, id_mask, vlen);
    }

  /* Decompose in_arg into n, B, r: in_arg = 2^n (1/B) (1 + r)
     B is equivalently defined by ind, 0 <= ind < 128 */
  VINT n = U_AS_I (
      __riscv_vadd (__riscv_vsrl (F_AS_U (vx), MAN_LEN - 1, vlen), 1, vlen));
  n = __riscv_vsra (n, 1, vlen);
  n = __riscv_vsub (n, EXP_BIAS, vlen);
  vx = U_AS_F (__riscv_vsrl (__riscv_vsll (F_AS_U (vx), BIT_WIDTH - MAN_LEN, vlen),
                            BIT_WIDTH - MAN_LEN, vlen));
  vx = U_AS_F (__riscv_vadd (F_AS_U (vx), (uint64_t)EXP_BIAS << MAN_LEN, vlen));
  n = __riscv_vsub (n, n_adjust, vlen);
  VFLOAT n_flt = __riscv_vfcvt_f (n, vlen);
  VFLOAT B = __riscv_vfrec7 (vx, vlen);
  /* get 7 msb of mantissa, and left shift by 3 to get address */
  VUINT ind =
      __riscv_vand (__riscv_vsrl (F_AS_U (vx), MAN_LEN - 10, vlen), 0x3F8, vlen);
  /* adjust B to be 1.0 if ind == 0 */
  VBOOL adjust_B = __riscv_vmseq (ind, 0, vlen);
  B = __riscv_vfmerge (B, fp_posOne, adjust_B, vlen);
  /* finally get r = B * in_arg - 1.0 */
  VFLOAT r = VFMV_VF (fp_posOne, vlen);
  r = __riscv_vfmsac (r, vx, B, vlen);

  /* Base-2 log is log2(in_arg) = log2(2^n * 1/B) + log2(1 + r)
     (n + log2(1/B))*log2(2) + log(1+r)*log2(e)
     log2(1/B) is stored in a table
     and log(1+r) is approximated by r + poly
     poly is a polynomial in r in the form r^2 * (p0 + p1 r + ... )
     To deliver this result accurately, one uses log2(2) and log2(e)
     with extra precision and sums the various terms in an appropriate order */
  VFLOAT rsq = __riscv_vfmul (r, r, vlen);
  VFLOAT rcube = __riscv_vfmul (rsq, r, vlen);

  VFLOAT poly_right = PSTEP (
      0x1.9999998877038p-3, r,
      PSTEP (-0x1.555c54f8b7c6cp-3, 0x1.2499765b3c27ap-3, r, vlen), vlen);

  VFLOAT poly_left = PSTEP (
      -0x1.000000000001cp-1, r,
      PSTEP (0x1.55555555555a9p-2, -0x1.fffffff2018cfp-3, r, vlen), vlen);

  VFLOAT poly = __riscv_vfmadd (poly_right, rcube, poly_left, vlen);
  poly = __riscv_vfmul (rsq, poly, vlen);
  /* log_e(1+r) is r + poly */

  /* Load table values and get n_flt + T to be A + a */
  VINT T = __riscv_vluxei64 (d->table, ind, vlen);
  VINT T_hi = __riscv_vsll (__riscv_vsra (T, 24, vlen), 24, vlen);
  VINT T_lo = __riscv_vsub (T, T_hi, vlen);
  VFLOAT T_hi_flt = __riscv_vfcvt_f (T_hi, vlen);
  VFLOAT A = __riscv_vfmadd (T_hi_flt, 0x1.0p-63, n_flt, vlen);
  VFLOAT a = __riscv_vfcvt_f (T_lo, vlen);
  a = __riscv_vfmul (a, 0x1.0p-63, vlen);

  /* Compute (A + a) * (log2_2_hi + log2_2_lo) + (r + P) * (log2_e_hi +
     log2_e_lo) where log2_2 = 1.0 and log2_e = 1/ln(2) */
  VFLOAT delta_1 = __riscv_vfmul (A, 0.0, vlen);
  delta_1 = __riscv_vfmadd (a, 0x1.0p0, delta_1, vlen);
  delta_1 = __riscv_vfmacc (delta_1, RVV_LN2_INV_LO_D, r, vlen);
  poly = __riscv_vfmadd (poly, RVV_LN2_INV_D, delta_1, vlen);

  poly = __riscv_vfmacc (poly, RVV_LN2_INV_D, r, vlen);

  vy = __riscv_vfadd (A, poly, vlen);
  vy = __riscv_vmerge (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
