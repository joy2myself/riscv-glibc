/* Double-precision vector log1p function.
   log1p(x) computes log(1 + x) accurately for small x

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

/* This version uses a table-based approach for log1p */

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT vx = x;
  VFLOAT vxp1, vy, vy_special = VFMV_VF (fp_posZero, vlen);
  VBOOL special_args;
  VINT n_adjust;
  const struct v_log_data *d = ptr_barrier (&__v_log_data);

  SET_ROUNDTONEAREST;

  /* Compute 1 + x */
  vxp1 = __riscv_vfadd (vx, fp_posOne, vlen);

  /* Handle special cases: NaN, Inf, x < -1, x = -1 */
  VUINT vclass = __riscv_vfclass (vxp1, vlen);
  IDENTIFY (vclass, 0x3BF, special_args, vlen);

  UINT nb_special_args = __riscv_vcpop (special_args, vlen);
  if (nb_special_args > 0)
    {
      VBOOL id_mask;
      /* log1p(x < -1) = NaN */
      IDENTIFY (vclass, class_negative, id_mask, vlen);
      vxp1 = __riscv_vfmerge (vxp1, fp_sNaN, id_mask, vlen);
      /* log1p(-1) = -Inf */
      IDENTIFY (vclass, class_Zero, id_mask, vlen);
      vxp1 = __riscv_vfmerge (vxp1, fp_negZero, id_mask, vlen);
      /* Eliminate positive denorm input from special_args */
      IDENTIFY (vclass, 0x39F, special_args, vlen);
      vy_special = __riscv_vfrec7 (special_args, vxp1, vlen);
      vy_special = __riscv_vfadd (special_args, vy_special, vxp1, vlen);
      vx = __riscv_vfmerge (vx, fp_posOne, special_args, vlen);
      vxp1 = __riscv_vfmerge (vxp1, 0x1.0p1, special_args, vlen);
    }

  /* Compute 1+x as vxp1 + u for better precision */
  VFLOAT first = __riscv_vfmax (vx, fp_posOne, vlen);
  VFLOAT second = __riscv_vfmin (vx, fp_posOne, vlen);
  VFLOAT u = __riscv_vfsub (first, vxp1, vlen);
  u = __riscv_vfadd (u, second, vlen);

  /* Decompose 1+x = 2^n * (1/B) * (1 + r) */
  n_adjust = VMVI_VX (0, vlen);
  VFLOAT scale = VFMV_VF (fp_posOne, vlen);
  VBOOL adjust = __riscv_vmfge (vxp1, 0x1.0p1020, vlen);
  n_adjust = __riscv_vmerge (n_adjust, 4, adjust, vlen);
  scale = __riscv_vfmerge (scale, 0x1.0p-4, adjust, vlen);
  vxp1 = __riscv_vfmul (vxp1, scale, vlen);

  VINT n = U_AS_I (
      __riscv_vadd (__riscv_vsrl (F_AS_U (vxp1), MAN_LEN - 1, vlen), 1, vlen));
  n = __riscv_vsra (n, 1, vlen);
  n = __riscv_vsub (n, EXP_BIAS, vlen);
  n = __riscv_vadd (n, n_adjust, vlen);
  VFLOAT n_flt = __riscv_vfcvt_f (n, vlen);

  VUINT ind = __riscv_vand (__riscv_vsrl (F_AS_U (vxp1), MAN_LEN - 10, vlen),
                            0x3F8, vlen);

  /* Normalize exponent of vxp1. */
  VFLOAT B = __riscv_vfrec7 (vxp1, vlen);
  /* Adjust B if ind == 0, rounding up at the 7th bit. */
  VBOOL adjust_B = __riscv_vmseq (ind, 0, vlen);
  VINT B_bits = F_AS_I (B);
  VINT B_bits_adj = __riscv_vadd (B_bits, (int64_t) 0x0000200000000000, vlen);
  B_bits = __riscv_vmerge (B_bits, B_bits_adj, adjust_B, vlen);
  B = I_AS_F (B_bits);

  VFLOAT r_hi = VFMV_VF (fp_posOne, vlen);
  r_hi = __riscv_vfmsac (r_hi, vxp1, B, vlen);
  VFLOAT r = __riscv_vfmacc (r_hi, u, B, vlen);
  VFLOAT delta_r = __riscv_vfsub (r_hi, r, vlen);
  delta_r = __riscv_vfmacc (delta_r, u, B, vlen);

  VFLOAT rsq = __riscv_vfmul (r, r, vlen);
  VFLOAT rcube = __riscv_vfmul (rsq, r, vlen);

  VFLOAT poly_right = PSTEP (
      0x1.9999998877038p-3, r,
      PSTEP (-0x1.555c54f8b7c6cp-3, 0x1.2499765b3c27ap-3, r, vlen), vlen);

  VFLOAT poly_left = PSTEP (
      -0x1.000000000001cp-1, r,
      PSTEP (0x1.55555555555a9p-2, -0x1.fffffff2018cfp-3, r, vlen), vlen);

  VFLOAT poly = __riscv_vfmadd (poly_right, rcube, poly_left, vlen);
  poly = __riscv_vfmadd (rsq, poly, delta_r, vlen);

  /* Load table values and get n_flt + T to be A + a. */
  VINT T = __riscv_vluxei64 (d->table, ind, vlen);
  VINT T_hi = __riscv_vsll (__riscv_vsra (T, 24, vlen), 24, vlen);
  VINT T_lo = __riscv_vsub (T, T_hi, vlen);
  VFLOAT T_hi_flt = __riscv_vfcvt_f (T_hi, vlen);
  VFLOAT A = __riscv_vfmadd (T_hi_flt, 0x1.0p-63, n_flt, vlen);
  VFLOAT a = __riscv_vfcvt_f (T_lo, vlen);
  a = __riscv_vfmul (a, 0x1.0p-63, vlen);

  /* Compute (A + a) * (log2_hi + log2_lo) + (r + poly). */
  VFLOAT delta_1 = __riscv_vfmul (A, RVV_LN2_LO_D, vlen);
  delta_1 = __riscv_vfmadd (a, RVV_LN2_HI_D, delta_1, vlen);
  poly = __riscv_vfadd (poly, delta_1, vlen);
  poly = __riscv_vfadd (poly, r, vlen);
  vy = __riscv_vfmadd (A, RVV_LN2_HI_D, poly, vlen);

  VFLOAT abs_x = __riscv_vfabs (x, vlen);
  VBOOL tiny_x = __riscv_vmflt (abs_x, 0x1.0p-54, vlen);
  vy = VMERGE_VVM_F (vy, x, tiny_x, vlen);
  vy = __riscv_vmerge (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
