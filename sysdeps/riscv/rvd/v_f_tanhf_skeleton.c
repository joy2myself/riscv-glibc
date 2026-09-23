/* Single-precision vector tanhf function.
   Based on SLEEF xtanhf (u10 = 1.0 ULP precision).

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
  VFLOAT y;
  VFLOAT y_special = VFMV_VF (fp_posZero_f, vlen);
  VBOOL special_args, small_x;
  UINT nb_special_args;

  SET_ROUNDTONEAREST;

  VFLOAT abs_x = __riscv_vfabs (vx_orig, vlen);
  VUINT vclass = __riscv_vfclass (vx_orig, vlen);
  VBOOL NaN_Inf;
  IDENTIFY_F (vclass, class_NaN | class_Inf, NaN_Inf, vlen);
  VBOOL x_nan;
  IDENTIFY_F (vclass, class_NaN, x_nan, vlen);

  VUINT expo_x = __riscv_vsrl (F_AS_U_F (abs_x), MAN_LEN_F, vlen);
  expo_x = __riscv_vand (expo_x, 0xff, vlen);
  small_x = __riscv_vmsleu
    (expo_x, EXP_BIAS_F - MAN_LEN_F - 5, vlen);

  /* Ordered comparisons with NaNs raise FE_INVALID on RVV.  Preserve Inf
     for its valid saturation comparison, but sanitize NaN lanes first.  */
  VFLOAT abs_x_compare
    = VFMERGE_VFM (abs_x, fp_posZero_f, x_nan, vlen);
  VBOOL x_large = __riscv_vmfgt
    (abs_x_compare, RVV_TANHF_SATURATION_THRESHOLD_F, vlen);
  special_args = __riscv_vmor
    (NaN_Inf, __riscv_vmor (small_x, x_large, vlen), vlen);
  nb_special_args = __riscv_vcpop (special_args, vlen);

  if (nb_special_args > 0)
    {
      /* tanhf(x) rounds to x for tiny finite x.  Preserve the scalar
         inexact/underflow side effect without sending these lanes through the
         exponential and reciprocal core.  */
      VFLOAT tiny_result = __riscv_vfmadd
        (small_x, vx_orig, RVV_LIBMVEC_HYPERBOLIC_TINY_CORRECTION_F,
         vx_orig, vlen);
      y_special = VMERGE_VVM_F
        (y_special, tiny_result, small_x, vlen);

      /* tanh(+-Inf) = +-1.  */
      VFLOAT signed_one = __riscv_vfsgnj
        (VFMV_VF_F (fp_posOne_f, vlen), vx_orig, vlen);
      y_special = VMERGE_VVM_F
        (y_special, signed_one, x_large, vlen);

      /* Quiet and propagate NaNs through the EH path.  */
      VFLOAT nan_result = vx_orig;
      y_special = VMERGE_VVM_F
        (y_special, nan_result, x_nan, vlen);

      /* A quiet NaN must not manufacture FE_INVALID.  Only signaling NaNs
         need an arithmetic quieting operation, which also preserves the
         scalar invalid-exception behavior.  */
      VBOOL x_snan;
      IDENTIFY_F (vclass, class_sNaN, x_snan, vlen);
      VFLOAT snan_result = __riscv_vfadd
        (x_snan, vx_orig, vx_orig, vlen);
      y_special = VMERGE_VVM_F
        (y_special, snan_result, x_snan, vlen);
    }

  VFLOAT vx_core = VFMERGE_VFM
    (vx_orig, fp_posZero_f, special_args, vlen);

  /* y = |x|.  */
  abs_x = __riscv_vfabs (vx_core, vlen);

  /* Compute tanh using expk2f: tanh(x) = (d - e) / (d + e)
     where d = exp(|x|) and e = 1/d.  */
  y = abs_x;

  /* d = expk2f(y).  */
  VFLOAT d_hi, d_lo;
  RVV_LIBMVEC_IMPL_NAME (sinhf_exp_pair, LMUL)
    (y, VFMV_VF (0.0f, vlen), &d_hi, &d_lo, vlen);

  /* e = 1 / d.  */
  VFLOAT e_hi, e_lo;
  DFREC_VF2_VF2 (d_hi, d_lo, e_hi, e_lo, vlen);

  /* neg_e = -e.  */
  VFLOAT neg_e_hi, neg_e_lo;
  RVV_SLEEF_ORIG_DFNEG_VF2 (e_hi, e_lo, neg_e_hi, neg_e_lo, vlen);

  /* numerator = d + (-e) = d - e.  */
  VFLOAT num_hi, num_lo;
  RVV_SLEEF_ORIG_DFADD2_VF2_VF2 (d_hi, d_lo, neg_e_hi, neg_e_lo, num_hi, num_lo, vlen);

  /* denominator = d + e.  */
  VFLOAT den_hi, den_lo;
  RVV_SLEEF_ORIG_DFADD2_VF2_VF2 (d_hi, d_lo, e_hi, e_lo, den_hi, den_lo, vlen);

  /* result = numerator / denominator.  */
  VFLOAT result_hi, result_lo;
  DFDIV_VF2_VF2_VF2
    (num_hi, num_lo, den_hi, den_lo, result_hi, result_lo, vlen);

  /* y = result_hi + result_lo.  */
  y = VFADD_VV (result_hi, result_lo, vlen);

  /* If |x| > 8.664339742 or y is NaN, return 1.  */
  VBOOL y_is_nan = __riscv_vmfne (y, y, vlen);
  VBOOL use_one = __riscv_vmor (x_large, y_is_nan, vlen);
  y = VFMERGE_VFM (y, 1.0f, use_one, vlen);

  /* Apply sign.  */
  y = __riscv_vfsgnjx (y, vx_orig, vlen);
  y = VMERGE_VVM_F (y, y_special, special_args, vlen);

  RESTORE_FRM;
  return y;
}

#undef LMUL
#undef SIMDLEN
