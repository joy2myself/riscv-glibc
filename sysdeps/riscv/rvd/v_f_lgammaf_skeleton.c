/* Single-precision vector lgamma function.
   lgammaf(x) computes the natural logarithm of the absolute value of the
   gamma function.

   The finite kernel is the pinned SLEEF gammafk/xlgammaf_u1 graph.  The
   surrounding classification and exception handling follows the clean
   VECLIBM-style path used by the other RVV routines.

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

#include "v_f_lgammaf_clean_repair.c"

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT x_orig = x;
  VFLOAT y, y_special;

  SET_ROUNDTONEAREST;

  /* Classify first.  All comparisons and the integer conversion below use
     sanitized operands for unordered lanes, so sNaN/Inf lanes do not leak
     signaling operations into the finite preparation path.  */
  VUINT xclass = __riscv_vfclass (x_orig, vlen);
  VBOOL x_is_nan, x_is_inf, x_is_subnormal;
  IDENTIFY_F (xclass, class_NaN, x_is_nan, vlen);
  IDENTIFY_F (xclass, class_Inf, x_is_inf, vlen);
  IDENTIFY_F (xclass, class_Subnormal, x_is_subnormal, vlen);
  VBOOL x_is_nan_inf = __riscv_vmor (x_is_nan, x_is_inf, vlen);
  VBOOL x_is_number = __riscv_vmnot (x_is_nan_inf, vlen);

  VFLOAT ax = __riscv_vfabs (x_orig, vlen);
  VFLOAT ax_cmp = VFMERGE_VFM (ax, fp_posOne_f, x_is_nan_inf, vlen);
  VFLOAT x_cmp = VFMERGE_VFM (x_orig, fp_posOne_f, x_is_nan_inf, vlen);

  VBOOL x_is_zero = __riscv_vmfeq (ax_cmp, fp_posZero_f, vlen);
  VBOOL ax_ge_2p24 = __riscv_vmfge (ax_cmp, 0x1.0p24f, vlen);
  VBOOL ax_lt_2p24 = __riscv_vmflt (ax_cmp, 0x1.0p24f, vlen);
  VBOOL safe_int = __riscv_vmand (ax_lt_2p24, x_is_number, vlen);
  VFLOAT x_for_int
    = VFMERGE_VFM (x_orig, fp_posZero_f, __riscv_vmnot (safe_int, vlen),
                   vlen);
  VINT x_trunc_i = RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_tgammaf_truncate_m, LMUL) (x_for_int, vlen);
  VFLOAT x_trunc = __riscv_vfcvt_f (x_trunc_i, vlen);
  VBOOL x_is_int_small = __riscv_vmfeq (x_trunc, x_for_int, vlen);
  x_is_int_small = __riscv_vmand (x_is_int_small, safe_int, vlen);
  VBOOL x_is_int_large = __riscv_vmand (ax_ge_2p24, x_is_number, vlen);
  VBOOL x_is_int = __riscv_vmor (x_is_int_small, x_is_int_large, vlen);

  VUINT x_sign = __riscv_vand (F_AS_U_F (x_orig), 0x80000000u, vlen);
  VBOOL x_is_neg
    = __riscv_vmseq (U_AS_I_F (x_sign), VMVI_VX_F (0x80000000, vlen), vlen);
  x_is_neg = __riscv_vmand (x_is_neg, x_is_number, vlen);
  VBOOL x_is_pole = __riscv_vmor
    (x_is_zero, __riscv_vmand (x_is_neg, x_is_int, vlen), vlen);

  VBOOL x_eq_1 = __riscv_vmfeq (x_cmp, 1.0f, vlen);
  x_eq_1 = __riscv_vmand (x_eq_1, x_is_number, vlen);
  VBOOL x_eq_2 = __riscv_vmfeq (x_cmp, 2.0f, vlen);
  x_eq_2 = __riscv_vmand (x_eq_2, x_is_number, vlen);
  VBOOL x_is_1_or_2 = __riscv_vmor (x_eq_1, x_eq_2, vlen);

  VBOOL special_args = __riscv_vmor (x_is_nan_inf, x_is_pole, vlen);
  special_args = __riscv_vmor (special_args, x_is_1_or_2, vlen);
  VBOOL quarantine_args = __riscv_vmor (special_args, x_is_subnormal, vlen);

  UINT n_special = __riscv_vcpop (quarantine_args, vlen);
  y_special = VFMV_VF_F (0.0f, vlen);
  if (n_special > 0)
    {
      /* NaN provider: operate only on NaN lanes so a discarded finite lane
         cannot contribute an unrelated overflow or inexact flag.  */
      VFLOAT nan_input = VFMERGE_VFM
        (x_orig, fp_posZero_f, __riscv_vmnot (x_is_nan, vlen), vlen);
      VFLOAT nan_value = __riscv_vfadd (nan_input, nan_input, vlen);
      y_special = VMERGE_VVM_F (y_special, nan_value, x_is_nan, vlen);

      /* Poles are +Inf and must raise FE_DIVBYZERO.  The inactive lanes use
         one/one, which is exact and keeps the provider lane-local.  */
      VFLOAT pole_denom = VFMERGE_VFM
        (VFMV_VF_F (1.0f, vlen), fp_posZero_f, x_is_pole, vlen);
      VFLOAT pole_value = __riscv_vfdiv
        (VFMV_VF_F (1.0f, vlen), pole_denom, vlen);
      y_special = VMERGE_VVM_F (y_special, pole_value, x_is_pole, vlen);

      /* x == 1 and x == 2 are exact zeroes.  */
      y_special = VMERGE_VVM_F
        (y_special, VFMV_VF_F (0.0f, vlen), x_is_1_or_2, vlen);

      /* lgammaf(+/-Inf) = +Inf and does not require an invalid provider.  */
      y_special = VFMERGE_VFM
        (y_special, fp_posInf_f, x_is_inf, vlen);

      /* Subnormal finite inputs have lgammaf(x) = -log(|x|) to binary32
         precision.  The provider normalizes before logk2f, keeping the
         clean path free of source reciprocal/scale exceptions.  */
      if (__riscv_vcpop (x_is_subnormal, vlen) > 0)
        {
          VFLOAT subnormal_input = VFMERGE_VFM
            (__riscv_vfabs (x_orig, vlen), 0x1.0p-64f,
             __riscv_vmnot (x_is_subnormal, vlen), vlen);
          VFLOAT subnormal_value = RVV_LGAMMAF_REPAIR_CALL
            (rvv_lgammaf_clean_subnormal, LMUL) (subnormal_input, vlen);
          y_special = VMERGE_VVM_F
            (y_special, subnormal_value, x_is_subnormal, vlen);
        }

      /* The finite kernel is source-shaped, but special lanes are never
         allowed to enter it with their original exceptional values.  */
      x = VFMERGE_VFM (x_orig, fp_posOne_f, quarantine_args, vlen);
    }

  VBOOL main_args = __riscv_vmnot (quarantine_args, vlen);
  if (__riscv_vcpop (main_args, vlen) == 0)
    {
      RESTORE_FRM;
      return y_special;
    }

  /* This is the pinned SLEEF xlgammaf_u1 finite graph:

       gammafk(x)
       -> dfdiv(clln, clld)
       -> dfabs
       -> logk2f
       -> dfadd(clc, log(abs(clln/clld)))
       -> collapse

     No Lanczos replacement, local log/sin approximation, or LMUL-specific
     precision split belongs in this clean routine.  Only the outer EH path
     above and the final source-compatible selection below are clean-side
     adaptations.  */
  /* Keep the source graph away from the large-x intermediate overflow in
     (x - 1/2) * log(x).  The clean-only stable rearrangement below is used
     for finite values beginning before that source intermediate overflows;
     expected-overflow lanes are quarantined and provided separately.  */
  VBOOL large_positive = __riscv_vmand
    (main_args, __riscv_vmfge (x_cmp, 0x1.8p+121f, vlen), vlen);
  VBOOL large_overflow = __riscv_vmand
    (main_args, __riscv_vmfge (x_cmp, 0x1.8964cp+121f, vlen), vlen);
  VBOOL large_finite = __riscv_vmand
    (large_positive, __riscv_vmnot (large_overflow, vlen), vlen);
  VFLOAT kernel_x = VFMERGE_VFM
    (x, 4.0f, large_positive, vlen);

  VFLOAT clc_hi, clc_lo, clln_hi, clln_lo, clld_hi, clld_lo;
  RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_tgammaf_gammafk_m, LMUL)
    (kernel_x, &clc_hi, &clc_lo, &clln_hi, &clln_lo, &clld_hi,
     &clld_lo, vlen);

  VFLOAT quotient_hi, quotient_lo;
  RVV_SLEEF_ORIG_DFDIV_VF2_VF2
    (clln_hi, clln_lo, clld_hi, clld_lo, quotient_hi, quotient_lo, vlen);

  VFLOAT quotient_abs_hi = __riscv_vfabs (quotient_hi, vlen);
  VUINT quotient_sign = __riscv_vand
    (F_AS_U_F (quotient_hi), 0x80000000u, vlen);
  VFLOAT quotient_abs_lo = U_AS_F_F
    (__riscv_vxor (F_AS_U_F (quotient_lo), quotient_sign, vlen));

  VFLOAT log_abs_hi, log_abs_lo;
  RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_tgammaf_logk2f_m, LMUL)
    (quotient_abs_hi, quotient_abs_lo, &log_abs_hi, &log_abs_lo, vlen);

  VFLOAT sum_hi, sum_lo;
  RVV_SLEEF_ORIG_DFADD2_VF2_VF2
    (clc_hi, clc_lo, log_abs_hi, log_abs_lo, sum_hi, sum_lo, vlen);
  y = VFADD_VV (sum_hi, sum_lo, vlen);

  /* Preserve xlgammaf_u1's late policy: +/-Inf input, non-positive integer
     input, or a finite input whose reduction produced NaN maps to +Inf.  */
  VBOOL negative_integer = __riscv_vmand
    (__riscv_vmfle (x_cmp, fp_posZero_f, vlen), x_is_int, vlen);
  negative_integer = __riscv_vmand
    (negative_integer, x_is_number, vlen);
  VBOOL reduction_nan = __riscv_vmfne (y, y, vlen);
  VBOOL finite_nan = __riscv_vmand (x_is_number, reduction_nan, vlen);
  VBOOL inf_result = __riscv_vmor
    (x_is_inf, __riscv_vmor (negative_integer, finite_nan, vlen), vlen);
  y = VFMERGE_VFM (y, fp_posInf_f, inf_result, vlen);

  if (__riscv_vcpop (large_positive, vlen) > 0)
    {
      VFLOAT large_x = VFMERGE_VFM (x, 4.0f,
                                    __riscv_vmnot (large_finite, vlen),
                                    vlen);
      VFLOAT stable_y = RVV_LGAMMAF_REPAIR_CALL
        (rvv_lgammaf_clean_large, LMUL) (large_x, vlen);
      y = VMERGE_VVM_F (y, stable_y, large_finite, vlen);

      /* Generate +Inf only on lanes whose correctly rounded result is
         outside the finite range.  This also preserves the required
         overflow/inexact provider without evaluating the stable core on
         those lanes.  */
      VFLOAT overflow_x = VFMERGE_VFM (x, 4.0f,
                                       __riscv_vmnot (large_overflow, vlen),
                                       vlen);
      VFLOAT overflow_value = VFMUL_VF_M
        (large_overflow, overflow_x, 0x1.fffffep+127f, vlen);
      y = VMERGE_VVM_F (y, overflow_value, large_overflow, vlen);
    }

  /* The source reflection expression is accurate in its individual terms,
     but loses more than one ulp when the final value is close to a zero.
     Re-evaluate only those negative lanes with the widened-double clean
     fallback; all other lanes are sanitized before the helper.  The widened
     fallback is intentionally unavailable at LMUL 8.  */
#if LMUL != 8
  VBOOL negative_repair = __riscv_vmand
    (x_is_neg, __riscv_vmand
       (main_args, __riscv_vmflt (__riscv_vfabs (y, vlen), 0x1.0p-5f, vlen),
        vlen), vlen);
  if (__riscv_vcpop (negative_repair, vlen) > 0)
    {
      VFLOAT reflection_x = VFMERGE_VFM
        (x, -2.5f, __riscv_vmnot (negative_repair, vlen), vlen);
      VFLOAT reflection_y = RVV_LGAMMAF_REPAIR_CALL
        (rvv_lgammaf_clean_reflection, LMUL) (reflection_x, vlen);
      y = VMERGE_VVM_F (y, reflection_y, negative_repair, vlen);
    }
#endif

  y = VMERGE_VVM_F (y, y_special, quarantine_args, vlen);

  RESTORE_FRM;
  return y;
}

#undef LMUL
#undef SIMDLEN
