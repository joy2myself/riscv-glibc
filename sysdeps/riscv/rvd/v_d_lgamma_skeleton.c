/* Double-precision vector lgamma function.
   lgamma(x) computes the natural logarithm of the absolute value of the gamma function

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
  VFLOAT vx, vx_orig, vy, vy_special;
  VBOOL special_args;
  VFLOAT zero;
  VINT expo_adj;
  const struct v_lgamma_data *data = ptr_barrier (&__v_lgamma_data);
  VBOOL lgamma_mask;

  SET_ROUNDTONEAREST;

  vx = x;
  vx_orig = x;
  zero = VFMV_VF (fp_posZero, vlen);
  vy_special = VFMV_VF (fp_posZero, vlen);

  /* Classify special lanes up front, following the log/logf skeleton style:
     NaN propagates, +/-Inf and poles map to +Inf with the expected flags,
     tiny nonzero inputs take the -log(|x|) path, and very large finite inputs
     overflow directly to +Inf.  */
  VUINT vclass = __riscv_vfclass (vx_orig, vlen);
  VUINT expo_x
    = __riscv_vand (__riscv_vsrl (F_AS_U (vx_orig), MAN_LEN, vlen), 0x7FF,
                    vlen);
  VFLOAT ax = __riscv_vfabs (vx_orig, vlen);
  VBOOL x_small = __riscv_vmsltu (expo_x, EXP_BIAS - 60, vlen);
  VBOOL x_inf_nan = __riscv_vmseq (expo_x, 0x7FF, vlen);
  VFLOAT ax_cmp = __riscv_vfmerge (ax, fp_posZero, x_inf_nan, vlen);
  VFLOAT x_sign_cmp = __riscv_vfmerge (vx_orig, fp_posZero, x_inf_nan,
                                       vlen);
  VBOOL x_large = __riscv_vmfge (ax_cmp, LGAMMA_OVERFLOW_X, vlen);
  VBOOL x_neg = __riscv_vmflt (x_sign_cmp, fp_posZero, vlen);
  VBOOL x_neg_finite = __riscv_vmandn (x_neg, x_inf_nan, vlen);
  VBOOL x_neg_huge_int = __riscv_vmand (x_neg_finite,
                                        __riscv_vmfge (ax_cmp, 0x1.0p+52,
                                                       vlen),
                                        vlen);
  VBOOL x_neg_small = __riscv_vmandn (x_neg_finite, x_neg_huge_int, vlen);
  VFLOAT ax_int_input
      = VMERGE_VVM_F (VFMV_VF (fp_posZero, vlen), ax, x_neg_small, vlen);
  VINT ax_int = __riscv_vfcvt_x (ax_int_input, vlen);
  VFLOAT ax_roundtrip = __riscv_vfcvt_f (ax_int, vlen);
  VBOOL x_neg_small_int = __riscv_vmand (x_neg_small,
                                         __riscv_vmfeq (ax_int_input,
                                                        ax_roundtrip, vlen),
                                         vlen);
  VBOOL x_neg_int = __riscv_vmor (x_neg_huge_int, x_neg_small_int, vlen);
  VBOOL x_front_special = __riscv_vmor (x_small, x_inf_nan, vlen);
  special_args = __riscv_vmor (x_front_special, x_large, vlen);
  special_args = __riscv_vmor (special_args, x_neg_int, vlen);
  UINT nb_special_args = __riscv_vcpop (special_args, vlen);

  if (nb_special_args > 0)
    {
      VBOOL x_inf;
      IDENTIFY (vclass, class_Inf, x_inf, vlen);
      vx = __riscv_vfmerge (vx, fp_posInf, x_inf, vlen);

      VBOOL x_zero;
      IDENTIFY (vclass, class_Zero, x_zero, vlen);
      vx = __riscv_vfmerge (vx, fp_posZero, x_zero, vlen);

      VFLOAT y_tmp = vx;
      VFLOAT y_0 = __riscv_vfrec7 (x_zero, vx, vlen);
      y_tmp = __riscv_vmerge (y_tmp, y_0, x_zero, vlen);
      VFLOAT vy_front = VFADD_VV_M (x_front_special, vx, y_tmp, vlen);
      vy_special = __riscv_vmerge (vy_special, vy_front, x_front_special,
                                   vlen);

      VBOOL x_small_nonzero
        = __riscv_vmandn (x_small, __riscv_vmfeq (vx, fp_posZero, vlen), vlen);
      if (__riscv_vcpop (x_small_nonzero, vlen) > 0)
        {
          VFLOAT x_tmp = VFMV_VF (fp_posOne, vlen);
          x_tmp = __riscv_vmerge (x_tmp, vx, x_small_nonzero, vlen);
          x_tmp = __riscv_vfsgnj (x_tmp, fp_posOne, vlen);
          VFLOAT y_hi, y_lo;
          lgamma_mask = x_small_nonzero;
          LGAMMA_LOG (x_tmp, zero, zero, y_hi, y_lo, vlen);
          y_hi = LGAMMA_VFADD_VV_M (x_small_nonzero, y_hi, y_lo, vlen);
          y_hi = __riscv_vfsgnj (y_hi, fp_posOne, vlen);
          vy_special = __riscv_vmerge (vy_special, y_hi, x_small_nonzero,
                                       vlen);
        }

      vy_special = __riscv_vfmerge (vy_special, fp_posInf, x_large, vlen);
      VFLOAT y_neg_int = __riscv_vfrec7 (x_neg_int, zero, vlen);
      vy_special = __riscv_vmerge (vy_special, y_neg_int, x_neg_int, vlen);
      vx = __riscv_vfmerge (vx, fp_posOne, special_args, vlen);
    }

  lgamma_mask = __riscv_vmnot (special_args, vlen);

  /* Work with the absolute value, sign handled later for negative inputs. */
  vx = __riscv_vfabs (vx, vlen);
  VFLOAT vx_cmp = __riscv_vfmerge (vx_orig, fp_posZero, special_args, vlen);
  VBOOL x_lt_0 = __riscv_vmflt (vx_cmp, fp_posZero, vlen);
  x_lt_0 = __riscv_vmandn (x_lt_0, special_args, vlen);

  VBOOL x_lt_225 = __riscv_vmflt (vx, 0x1.2p+1, vlen);
  VFLOAT c = VFMV_VF (fp_posOne, vlen);
  c = __riscv_vfmerge (c, fp_posZero, x_lt_225, vlen);

  VFLOAT logx_hi, logx_lo;
  LGAMMA_LOG (vx, zero, c, logx_hi, logx_lo, vlen);

  VFLOAT y_left_hi = VFMV_VF (fp_posZero, vlen);
  VFLOAT y_left_lo = VFMV_VF (fp_posZero, vlen);
  if (__riscv_vcpop (x_lt_225, vlen) > 0)
    {
      /* 0 < x < 2.25 */
      VFLOAT vx_tmp = VFMV_VF (0x1.0p0, vlen);
      vx_tmp = __riscv_vmerge (vx_tmp, vx, x_lt_225, vlen);
      LGAMMA_LE_225 (vx_tmp, logx_hi, logx_lo, y_left_hi, y_left_lo, vlen);
    }

  VFLOAT stir_hi, stir_lo;
  VFLOAT stir_corr_hi, stir_corr_lo;
  expo_adj = VMVI_VX (0, vlen);
  LGAMMA_LOG_STIRLING (vx, logx_hi, logx_lo, stir_hi, stir_lo, expo_adj, vlen);
  LOG_STIRLING_CORRECTION (vx, stir_corr_hi, stir_corr_lo, vlen);

  VFLOAT loggamma_hi, loggamma_lo;
  KNUTH2SUM (stir_hi, stir_corr_hi, loggamma_hi, loggamma_lo, vlen);
  loggamma_lo = LGAMMA_VFADD (loggamma_lo, stir_corr_lo, vlen);
  loggamma_lo = LGAMMA_VFADD (loggamma_lo, stir_lo, vlen);

  loggamma_hi = __riscv_vmerge (loggamma_hi, y_left_hi, x_lt_225, vlen);
  loggamma_lo = __riscv_vmerge (loggamma_lo, y_left_lo, x_lt_225, vlen);

  if (__riscv_vcpop (x_lt_0, vlen) > 0)
    {
      VFLOAT vx_for_neg = VFMV_VF (0x1.0p-1, vlen);
      vx_for_neg = __riscv_vmerge (vx_for_neg, vx, x_lt_0, vlen);
      VFLOAT logsin_hi, logsin_lo;
      lgamma_mask = x_lt_0;
      LGAMMA_LOGSIN (vx_for_neg, logsin_hi, logsin_lo, vy_special, special_args,
                     vlen);
      lgamma_mask = __riscv_vmnot (special_args, vlen);

      VFLOAT A, a;
      LGAMMA_KNUTH2SUM_M (x_lt_0, loggamma_hi, logsin_hi, A, a, vlen);
      a = LGAMMA_VFADD_VV_M (x_lt_0, a, logsin_lo, vlen);
      a = LGAMMA_VFADD_VV_M (x_lt_0, a, loggamma_lo, vlen);
      VFLOAT sum_hi, sum_lo;
      LGAMMA_KNUTH2SUM_M (x_lt_0, A, a, sum_hi, sum_lo, vlen);
      A = sum_hi;
      a = sum_lo;
      A = __riscv_vfsgnjx (A, fp_negOne, vlen);
      a = __riscv_vfsgnjx (a, fp_negOne, vlen);
      loggamma_hi = __riscv_vmerge (loggamma_hi, A, x_lt_0, vlen);
      loggamma_lo = __riscv_vmerge (loggamma_lo, a, x_lt_0, vlen);
    }

  VFLOAT abs_loggamma = __riscv_vfabs (loggamma_hi, vlen);
  VBOOL lgamma_small_bnd = __riscv_vmflt (abs_loggamma, 0x1.74p-4, vlen);
  VBOOL lgamma_small_3 = __riscv_vmflt (abs_loggamma, 0x1.168p-4, vlen);
  VBOOL x_gt_neg3 = __riscv_vmfgt (vx_cmp, -0x1.8p+1, vlen);
  VBOOL x_lt_bnd
      = __riscv_vmflt (vx_cmp, -0x1.4e147ae147ae1p+1, vlen);
  VBOOL x_gt_bnd
      = __riscv_vmfgt (vx_cmp, -0x1.4e147ae147ae1p+1, vlen);
  VBOOL x_lt_neg2 = __riscv_vmflt (vx_cmp, -0x1.0p1, vlen);
  VBOOL range_bnd_to_2 = __riscv_vmand (x_gt_bnd, x_lt_neg2, vlen);
  range_bnd_to_2 = __riscv_vmand (range_bnd_to_2, x_lt_0, vlen);
  range_bnd_to_2 = __riscv_vmand (range_bnd_to_2, lgamma_small_bnd, vlen);
  VBOOL range_3_to_bnd = __riscv_vmand (x_gt_neg3, x_lt_bnd, vlen);
  range_3_to_bnd = __riscv_vmand (range_3_to_bnd, x_lt_0, vlen);
  range_3_to_bnd = __riscv_vmand (range_3_to_bnd, lgamma_small_3, vlen);

  const int len_bnd_to_2 = (int) (sizeof (data->neg_bnd_to_2_c)
                                  / sizeof (data->neg_bnd_to_2_c[0]));
  const int len_3_to_bnd = (int) (sizeof (data->neg_3_to_bnd_c)
                                  / sizeof (data->neg_3_to_bnd_c[0]));
  RVV_LIBMVEC_IMPL_NAME (lgamma_zero_apply, LMUL)
    (vx_orig, range_bnd_to_2, RVV_LGAMMA_NEG_BND_TO_2_X0_0_D,
     RVV_LGAMMA_NEG_BND_TO_2_X0_1_D,
     RVV_LGAMMA_NEG_BND_TO_2_X0_2_D,
     data->neg_bnd_to_2_c, len_bnd_to_2, LGAMMA_NEG_BND_TO_2_SC,
     &loggamma_hi, &loggamma_lo, vlen);
  RVV_LIBMVEC_IMPL_NAME (lgamma_zero_apply, LMUL)
    (vx_orig, range_3_to_bnd, RVV_LGAMMA_NEG_3_TO_BND_X0_0_D,
     RVV_LGAMMA_NEG_3_TO_BND_X0_1_D,
     RVV_LGAMMA_NEG_3_TO_BND_X0_2_D,
     data->neg_3_to_bnd_c, len_3_to_bnd, LGAMMA_NEG_3_TO_BND_SC,
     &loggamma_hi, &loggamma_lo, vlen);

  loggamma_hi = LGAMMA_VFADD (loggamma_hi, loggamma_lo, vlen);
  expo_adj = __riscv_vadd (expo_adj, EXP_BIAS, vlen);
  vy = LGAMMA_VFMUL
    (loggamma_hi, I_AS_F (__riscv_vsll (expo_adj, MAN_LEN, vlen)), vlen);

  vy = __riscv_vmerge (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
