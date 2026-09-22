/* Single-precision vector atan2 function.

   Based on SLEEF xatan2f_u1 (u1 precision).

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
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT vy, VFLOAT vx, size_t vlen)
{
  VFLOAT vx_orig = vx;
  VFLOAT vy_orig = vy;
  SET_ROUNDTONEAREST;

  VUINT ux = F_AS_U (vx_orig);
  VUINT uy = F_AS_U (vy_orig);
  VUINT ax_bits = __riscv_vand (ux, 0x7fffffff, vlen);
  VUINT ay_bits = __riscv_vand (uy, 0x7fffffff, vlen);
  VBOOL x_is_inf = __riscv_vmseq (ax_bits, 0x7f800000, vlen);
  VBOOL y_is_inf = __riscv_vmseq (ay_bits, 0x7f800000, vlen);
  VBOOL x_is_zero = __riscv_vmseq (ax_bits, 0, vlen);
  VBOOL y_is_zero = __riscv_vmseq (ay_bits, 0, vlen);
  VBOOL x_is_nan = __riscv_vmsgtu (ax_bits, 0x7f800000, vlen);
  VBOOL y_is_nan = __riscv_vmsgtu (ay_bits, 0x7f800000, vlen);
  VBOOL x_sign = __riscv_vmsne (__riscv_vand (ux, 0x80000000, vlen), 0, vlen);
  VBOOL any_nan = __riscv_vmor (x_is_nan, y_is_nan, vlen);
  VBOOL special_args = __riscv_vmor (x_is_inf, y_is_inf, vlen);
  special_args = __riscv_vmor (special_args, x_is_zero, vlen);
  special_args = __riscv_vmor (special_args, y_is_zero, vlen);
  special_args = __riscv_vmor (special_args, any_nan, vlen);
#if defined(COMPILE_FOR_ATAN2PI)
  /* atan2pif has exact normalized results for finite equal magnitudes:
     +/-1/4 when x is positive and +/-3/4 when x is negative.  Keep those
     lanes out of the approximation core so they do not raise an inexact
     flag merely while reconstructing an exact quadrant result.  */
  VBOOL equal_abs = __riscv_vmseq (ax_bits, ay_bits, vlen);
  VBOOL finite_equal = __riscv_vmandn (equal_abs, special_args, vlen);
  special_args = __riscv_vmor (special_args, finite_equal, vlen);
#endif

  /* Prepare special results before sanitizing the finite-core inputs.  This
     keeps the clean path in the VECLIBM-shaped order:
     classify -> prepare -> sanitize -> finite core -> final merge.  */
  VFLOAT r_special = VFMV_VF (0.0f, vlen);

#if defined(COMPILE_FOR_ATAN2PI)
  VFLOAT m = __riscv_vfsgnj (VFMV_VF (0.5f, vlen), vx_orig, vlen);
  VFLOAT x_inf_part = __riscv_vfsub (
    VFMV_VF (0.5f, vlen),
    RVV_LIBMVEC_IMPL_NAME (visinf2_f, LMUL) (vx_orig, m, vlen), vlen);
#else
  VFLOAT m = __riscv_vfsgnj (VFMV_VF (RVV_ATANF_PI2_F, vlen), vx_orig, vlen);
  VFLOAT x_inf_part = __riscv_vfsub (
    VFMV_VF (RVV_ATANF_PI2_F, vlen),
    RVV_LIBMVEC_IMPL_NAME (visinf2_f, LMUL) (vx_orig, m, vlen), vlen);
#endif
  VBOOL x_inf_or_zero = __riscv_vmor (x_is_inf, x_is_zero, vlen);
  r_special = __riscv_vmerge (r_special, x_inf_part, x_inf_or_zero, vlen);

#if defined(COMPILE_FOR_ATAN2PI)
  VFLOAT m2 = __riscv_vfsgnj (VFMV_VF (0.25f, vlen), vx_orig, vlen);
  VFLOAT y_inf_part = __riscv_vfsub (
    VFMV_VF (0.5f, vlen),
    RVV_LIBMVEC_IMPL_NAME (visinf2_f, LMUL) (vx_orig, m2, vlen), vlen);
#else
  VFLOAT m2 = __riscv_vfsgnj (VFMV_VF (RVV_ATANF_PI4_F, vlen), vx_orig, vlen);
  VFLOAT y_inf_part = __riscv_vfsub (
    VFMV_VF (RVV_ATANF_PI2_F, vlen),
    RVV_LIBMVEC_IMPL_NAME (visinf2_f, LMUL) (vx_orig, m2, vlen), vlen);
#endif
  r_special = __riscv_vmerge (r_special, y_inf_part, y_is_inf, vlen);

#if defined(COMPILE_FOR_ATAN2PI)
  VFLOAT r_xzero = __riscv_vfsgnj (VFMV_VF (0.5f, vlen), vy_orig, vlen);
#else
  VFLOAT r_xzero = __riscv_vfsgnj (VFMV_VF (RVV_ATANF_PI2_F, vlen), vy_orig,
                                   vlen);
#endif
  VBOOL x_only_zero = __riscv_vmandn (x_is_zero, y_is_zero, vlen);
  r_special = __riscv_vmerge (r_special, r_xzero, x_only_zero, vlen);

#if defined(COMPILE_FOR_ATAN2PI)
  VFLOAT r_zero = __riscv_vmerge (VFMV_VF (0.0f, vlen),
                                  VFMV_VF (1.0f, vlen), x_sign, vlen);
#else
  VFLOAT r_zero = __riscv_vmerge (VFMV_VF (0.0f, vlen),
                                  VFMV_VF (RVV_ATANF_PI_F, vlen), x_sign,
                                  vlen);
#endif
  r_special = __riscv_vmerge (r_special, r_zero, y_is_zero, vlen);
  r_special = __riscv_vfsgnj (r_special, vy_orig, vlen);

#if defined(COMPILE_FOR_ATAN2PI)
  VFLOAT equal_result = __riscv_vmerge (VFMV_VF (0.25f, vlen),
                                        VFMV_VF (0.75f, vlen), x_sign,
                                        vlen);
  equal_result = __riscv_vfsgnj (equal_result, vy_orig, vlen);
  r_special = __riscv_vmerge (r_special, equal_result, finite_equal, vlen);
#endif

#if !defined(COMPILE_FOR_ATAN2PI)
  /* The scalar atan2f path raises inexact for both signed-zero cases with
     x = -0, because the result is the pi branch.  Preserve that EH contract
     without changing the result bits: under the RNE envelope, adding a
     sub-half-ULP probe to +/-pi rounds back to the same value.  The
     normalized atan2pif endpoints are exact, so this probe must not be
     applied to that composition.  */
  VBOOL x_is_neg_zero = __riscv_vmand (x_is_zero, x_sign, vlen);
  VBOOL zero_pi = __riscv_vmand (x_is_neg_zero, y_is_zero, vlen);
  VFLOAT zero_pi_probe
    = VFADD_VF_M (zero_pi, r_special, 0x1p-25f, vlen);
  r_special = __riscv_vmerge (r_special, zero_pi_probe, zero_pi, vlen);
#endif

  /* Propagate the actual NaN-bearing operand with a masked x+x operation.
     This preserves qNaN behavior while allowing sNaN lanes to raise the
     required invalid flag without exposing them to the finite core.  */
  VFLOAT nan_input = __riscv_vmerge (vy_orig, vx_orig, x_is_nan, vlen);
  VFLOAT nan_result = VFADD_VV_M (any_nan, nan_input, nan_input, vlen);
  r_special = __riscv_vmerge (r_special, nan_result, any_nan, vlen);

  VBOOL main_args = __riscv_vmnot (special_args, vlen);
  if (__riscv_vcpop (main_args, vlen) == 0)
    {
      RESTORE_FRM;
      return r_special;
    }

  vx = __riscv_vfmerge (vx, fp_posOne, special_args, vlen);
  vy = __riscv_vfmerge (vy, fp_posZero, special_args, vlen);

  VFLOAT ax = __riscv_vfabs (vx, vlen);
  VBOOL tiny = __riscv_vmflt (ax, RVV_ATANF_TINY_F, vlen);

  if (__riscv_vcpop (tiny, vlen) > 0)
    {
      /* Do not compute discarded products in non-tiny lanes: a large finite
         operand times 2^24 can otherwise raise overflow/inexact before the
         result is selected away.  */
      vx = VFMUL_VF_M (tiny, vx, RVV_ATANF_SCALE_F, vlen);
      vy = VFMUL_VF_M (tiny, vy, RVV_ATANF_SCALE_F, vlen);
      ax = __riscv_vfabs (vx, vlen);
    }

  VFLOAT ay = __riscv_vfabs (vy, vlen);
  VFLOAT r_hi, r_lo;
  RVV_LIBMVEC_IMPL_NAME (atan2kf_u1, LMUL)
    (ay, VFMV_VF (0.0f, vlen), vx, VFMV_VF (0.0f, vlen), &r_hi, &r_lo,
     vlen);
#if defined(COMPILE_FOR_ATAN2PI)
  VFLOAT scaled_hi, scaled_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (r_hi, r_lo, VFMV_VF (RVV_ATANF_INV_PI_HI_F, vlen),
     VFMV_VF (RVV_ATANF_INV_PI_LO_F, vlen), scaled_hi, scaled_lo, vlen);
  VFLOAT r = VFADD_VV (scaled_hi, scaled_lo, vlen);
#else
  VFLOAT r = VFADD_VV (r_hi, r_lo, vlen);
#endif
  r = __riscv_vfsgnjx (r, vx, vlen);
  r = __riscv_vfsgnj (r, vy_orig, vlen);
  r = __riscv_vmerge (r, r_special, special_args, vlen);

  RESTORE_FRM;
  return r;
}

#undef LMUL
#undef SIMDLEN
