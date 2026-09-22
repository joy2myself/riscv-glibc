/* Single-precision vector hypotf function.
   hypotf(x, y) computes sqrt(x^2 + y^2) without overflow

   Based on SLEEF xhypotf_u05 (u05 = 0.5 ULP precision).

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

/* Shared float-float helpers live in v_math.h.  */

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT vx_orig, VFLOAT vy_orig,
                                           size_t vlen)
{
  VFLOAT vx, vy, vret, vret_special;
  VBOOL special_args;
  SET_ROUNDTONEAREST;

  vx = vx_orig;
  vy = vy_orig;

  /* Handle NaN/Inf lanes up front so they do not enter the ratio core.
     Match the double routine: preserve NaN propagation via x + y, but let
     Inf dominate unless a signaling NaN is present.  */
  VUINT vclass_x = __riscv_vfclass (vx_orig, vlen);
  VUINT vclass_y = __riscv_vfclass (vy_orig, vlen);
  VBOOL special_x, special_y;
  IDENTIFY (vclass_x, class_NaN | class_Inf, special_x, vlen);
  IDENTIFY (vclass_y, class_NaN | class_Inf, special_y, vlen);
  special_args = __riscv_vmor (special_x, special_y, vlen);

  vret_special = VFMV_VF (fp_posZero_f, vlen);
  if (__riscv_vcpop (special_args, vlen) > 0)
    {
      vret_special = VFADD_VV_M (special_args, vx_orig, vy_orig, vlen);

      VBOOL x_snan, y_snan, any_snan;
      IDENTIFY (vclass_x, class_sNaN, x_snan, vlen);
      IDENTIFY (vclass_y, class_sNaN, y_snan, vlen);
      any_snan = __riscv_vmor (x_snan, y_snan, vlen);

      VBOOL x_inf, y_inf, any_inf;
      IDENTIFY (vclass_x, class_Inf, x_inf, vlen);
      IDENTIFY (vclass_y, class_Inf, y_inf, vlen);
      any_inf = __riscv_vmor (x_inf, y_inf, vlen);

      VBOOL inf_no_snan
          = __riscv_vmand (any_inf, __riscv_vmnot (any_snan, vlen), vlen);
      vret_special = VFMERGE_VFM (vret_special, fp_posInf_f, inf_no_snan, vlen);

      /* Keep later core operations finite and non-zero on replaced lanes.  */
      vx = VFMERGE_VFM (vx, fp_posOne_f, special_args, vlen);
      vy = VFMERGE_VFM (vy, fp_posZero_f, special_args, vlen);
    }

  /* Take absolute values.  */
  vx = __riscv_vfabs (vx, vlen);
  vy = __riscv_vfabs (vy, vlen);

  /* Find min and max.  */
  VFLOAT vmin = __riscv_vfmin (vx, vy, vlen);
  VFLOAT vmax = __riscv_vfmax (vx, vy, vlen);

  /* Scale tiny and overflow-adjacent lanes by a common power of two.  */
  VFLOAT d = vmax;
  VBOOL all_zero = __riscv_vmfeq (vmax, fp_posZero_f, vlen);
  d = VFMERGE_VFM (d, fp_posOne_f, all_zero, vlen);
  VBOOL is_tiny = __riscv_vmflt (vmax, RVV_HYPOTF_CLEAN_SCALE_LIMIT_F, vlen);
  VBOOL is_large = __riscv_vmfge (vmax, RVV_HYPOTF_CLEAN_LARGE_LIMIT_F, vlen);
  VBOOL scale_mask = __riscv_vmor (is_tiny, is_large, vlen);
  UINT nb_scaled = __riscv_vcpop (scale_mask, vlen);
  VFLOAT core_scale = VFMV_VF_F (1.0f, vlen);
  core_scale = VFMERGE_VFM (core_scale, RVV_HYPOTF_CLEAN_SCALE_UP_F, is_tiny,
                            vlen);
  core_scale = VFMERGE_VFM (core_scale, RVV_HYPOTF_CLEAN_SCALE_DOWN_F, is_large,
                            vlen);

  /* If the exponent gap makes the smaller nonzero operand negligible, its
     squared contribution is far below one float ulp of the result.  Quarantine
     it before the DF division so the correction tree cannot raise an
     intermediate underflow.  Preserve the scalar hypotf inexact side effect
     with a masked half-ulp tie operation.  When max is itself tiny, the
     common scaling above keeps the ratio faithful instead.  */
  VUINT min_exp = __riscv_vsrl (F_AS_U_F (vmin), MAN_LEN_F, vlen);
  min_exp = __riscv_vand (min_exp, 0xff, vlen);
  VUINT max_exp = __riscv_vsrl (F_AS_U_F (vmax), MAN_LEN_F, vlen);
  max_exp = __riscv_vand (max_exp, 0xff, vlen);
  VUINT exponent_gap = __riscv_vsub (max_exp, min_exp, vlen);
  VBOOL exponent_negligible
      = __riscv_vmsgeu (exponent_gap, RVV_HYPOTF_CLEAN_EXPONENT_GAP, vlen);
  VBOOL min_nonzero = __riscv_vmfne (vmin, fp_posZero_f, vlen);
  VBOOL ratio_negligible
      = __riscv_vmand (exponent_negligible, min_nonzero, vlen);
  ratio_negligible
      = __riscv_vmand (ratio_negligible, __riscv_vmnot (is_tiny, vlen), vlen);
  VFLOAT vmin_core = VFMERGE_VFM (vmin, fp_posZero_f, ratio_negligible, vlen);
  if (__riscv_vcpop (ratio_negligible, vlen) > 0)
    {
      volatile VFLOAT force_inexact
          = VFADD_VF_M (ratio_negligible,
                        VFMV_VF_F (1.0f, vlen), RVV_HYPOTF_CLEAN_HALF_ULP_F,
                                   vlen);
      (void) force_inexact;
    }

  VFLOAT n = __riscv_vfmul (vmin_core, core_scale, vlen);
  d = __riscv_vfmul (d, core_scale, vlen);
  VFLOAT vmax_core = __riscv_vfmul (vmax, core_scale, vlen);

  /* Compute t = n / d as vfloat2.  */
  VFLOAT t_hi, t_lo;
  RVV_SLEEF_ORIG_DFDIV_VF2_VF2
    (n, VFMV_VF_F (fp_posZero_f, vlen), d,
     VFMV_VF_F (fp_posZero_f, vlen), t_hi, t_lo, vlen);

  /* Compute t^2 as vfloat2.  */
  VFLOAT t2_hi, t2_lo;
  RVV_SLEEF_ORIG_DFSQU_VF2 (t_hi, t_lo, t2_hi, t2_lo, vlen);

  /* Compute t^2 + 1 as vfloat2.  */
  VFLOAT sum_hi, sum_lo;
  RVV_SLEEF_ORIG_DFADD2_VF2_VF_S
    (t2_hi, t2_lo, fp_posOne_f, sum_hi, sum_lo, vlen);

  /* Compute sqrt(t^2 + 1) as vfloat2.  */
  VFLOAT sqrt_hi, sqrt_lo;
  RVV_SLEEF_ORIG_DFSQRT_VF2_VF2
    (sum_hi, sum_lo, sqrt_hi, sqrt_lo, vlen);

  /* Compute result = sqrt * max as vfloat2.  */
  VFLOAT result_hi, result_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (sqrt_hi, sqrt_lo, vmax_core, result_hi, result_lo, vlen);
  if (nb_scaled > 0)
    {
      VFLOAT result_scale = VFMV_VF_F (1.0f, vlen);
      result_scale = VFMERGE_VFM (result_scale, RVV_HYPOTF_CLEAN_SCALE_DOWN_F,
                                  is_tiny, vlen);
      result_scale = VFMERGE_VFM (result_scale, RVV_HYPOTF_CLEAN_SCALE_UP_F,
                                  is_large, vlen);
      VFLOAT result_hi_scaled
          = __riscv_vfmul (result_hi, result_scale, vlen);
      /* Fuse the low correction with the unscaled high result.  A separate
         low*scale would report an intermediate underflow; on overflow, an
         opposite-signed low infinity would instead create an invalid inf-inf.
         Masking the low term in high-infinite lanes preserves scalar
         overflow/inexact semantics.  */
      VUINT result_hi_class = __riscv_vfclass (result_hi_scaled, vlen);
      VBOOL result_hi_inf;
      IDENTIFY_F (result_hi_class, class_Inf, result_hi_inf, vlen);
      VFLOAT low_for_scale
          = VFMERGE_VFM (result_lo, fp_posZero_f, result_hi_inf, vlen);
      VFLOAT result_scaled
          = VFMADD_VV (low_for_scale, result_scale, result_hi_scaled, vlen);
      result_hi = __riscv_vmerge (result_hi, result_scaled, scale_mask, vlen);
      result_lo = VFMERGE_VFM (result_lo, fp_posZero_f, scale_mask, vlen);
    }

  /* Extract result from vfloat2.  */
  vret = __riscv_vfadd (result_hi, result_lo, vlen);

  /* Handle NaN result: if result is NaN, return Inf.  */
  VBOOL is_nan_result = __riscv_vmfne (vret, vret, vlen);
  vret = __riscv_vfmerge (vret, fp_posInf_f, is_nan_result, vlen);

  /* Special case: if min == 0, return max.  */
  VBOOL min_is_zero = __riscv_vmfeq (vmin, fp_posZero_f, vlen);
  vret = __riscv_vmerge (vret, vmax, min_is_zero, vlen);

  /* Merge special lanes back in after the finite core.  */
  vret = VMERGE_VVM_F (vret, vret_special, special_args, vlen);

  RESTORE_FRM;
  return vret;
}

#undef LMUL
#undef SIMDLEN
