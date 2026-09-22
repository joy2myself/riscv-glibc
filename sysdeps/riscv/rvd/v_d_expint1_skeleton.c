/* Double-precision vector expint1 function.

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

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT vx = x;
  VFLOAT vy, vy_special;
  VBOOL special_args;
  VINT n_adjust;

  SET_ROUNDTONEAREST;

  vy_special = VFMV_VF (fp_posZero, vlen);
  EXCEPTION_HANDLING (vx, special_args, vy_special, n_adjust, vlen);
  VUINT vclass_orig = __riscv_vfclass (x, vlen);
  VBOOL x_pos_inf;
  IDENTIFY (vclass_orig, class_posInf, x_pos_inf, vlen);
  vy_special = __riscv_vfmerge (vy_special, fp_posZero, x_pos_inf, vlen);

  VBOOL non_special_args = __riscv_vmnot (special_args, vlen);
  VBOOL args_lt_1 = __riscv_vmflt (vx, fp_posOne, vlen);
  args_lt_1 = __riscv_vmand (args_lt_1, non_special_args, vlen);
  VBOOL args_ge_1 = __riscv_vmfge (vx, fp_posOne, vlen);
  args_ge_1 = __riscv_vmand (args_ge_1, non_special_args, vlen);
  VFLOAT vy_xlt1 = VFMV_VF (fp_posZero, vlen);
  VFLOAT vy_xge1 = VFMV_VF (fp_posZero, vlen);

  if (__riscv_vcpop (args_lt_1, vlen) > 0)
    {
      VFLOAT vx_lt_1
          = __riscv_vfmerge (vx, RVV_EXPINT1_X_LT_ONE_SUBSTITUTE_D, args_ge_1,
                             vlen);
      VFLOAT neg_logx_hi, neg_logx_lo;
      EXPINT1_NEG_LOG (vx_lt_1, n_adjust, neg_logx_hi, neg_logx_lo, vlen);
      VFLOAT poly_hi, poly_lo;
      EXPINT_POLY_01 (vx_lt_1, poly_hi, poly_lo, vlen);

      VFLOAT AA, aa;
      KNUTH2SUM (neg_logx_hi, poly_hi, AA, aa, vlen);
      aa = __riscv_vfadd (aa, neg_logx_lo, vlen);
      aa = __riscv_vfadd (aa, poly_lo, vlen);
      vy_xlt1 = __riscv_vfadd (AA, aa, vlen);
    }

  if (__riscv_vcpop (args_ge_1, vlen) > 0)
    {
      VFLOAT vx_ge_1
          = __riscv_vfmerge (vx, RVV_EXPINT1_X_GE_ONE_SUBSTITUTE_D, args_lt_1,
                             vlen);
      vx_ge_1 = __riscv_vfmin (vx_ge_1, RVV_EXPINT1_X_GE_ONE_MAX_D, vlen);

      VFLOAT recip_x_hi, recip_x_lo;
      recip_x_hi = __riscv_vfrdiv (vx_ge_1, fp_posOne, vlen);
      recip_x_lo = VFMV_VF (fp_posOne, vlen);
      recip_x_lo = __riscv_vfnmsac (recip_x_lo, vx_ge_1, recip_x_hi, vlen);
      recip_x_lo = __riscv_vfmul (recip_x_hi, recip_x_lo, vlen);

      VFLOAT rat_hi, rat_lo;
      EXPINT1_RAT_GE1 (recip_x_hi, recip_x_lo, RVV_EXPINT1_RAT_SCALE_D, rat_hi,
                       rat_lo, vlen);

      VFLOAT rat_by_x_hi, rat_by_x_lo;
      PROD_X2Y2 (recip_x_hi, recip_x_lo, rat_hi, rat_lo, rat_by_x_hi,
                 rat_by_x_lo, vlen);

      VFLOAT exp_negx_hi, exp_negx_lo;
      VINT n;
      EXP_NEGX (vx_ge_1, n, exp_negx_hi, exp_negx_lo, vlen);

      VFLOAT result_hi, result_lo;
      PROD_X2Y2 (rat_by_x_hi, rat_by_x_lo, exp_negx_hi, exp_negx_lo, result_hi,
                 result_lo, vlen);
      vy_xge1 = __riscv_vfadd (result_hi, result_lo, vlen);
      FAST_LDEXP (vy_xge1, n, vlen);
    }

  vy = __riscv_vmerge (vy_xlt1, vy_xge1, args_ge_1, vlen);
  vy = __riscv_vmerge (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
