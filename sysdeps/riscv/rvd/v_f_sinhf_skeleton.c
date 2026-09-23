/* Single-precision vector sinhf function.
   Based on SLEEF xsinhf (u10 = 1.0 ULP precision).

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

/* These two quarantine boundaries belong only to this clean entry, not
   the shared sinhf/tanhf exponential component.  */
#define RVV_SINHF_CLEAN_TINY_INPUT_F 0x1.0p-12f
#define RVV_SINHF_CLEAN_RECIPROCAL_QUARANTINE_F 70.0f

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT vx = x;
  VFLOAT vx_orig, vy, vy_special, abs_x;
  VBOOL special_args, tiny_args;
  UINT nb_special_args;

  SET_ROUNDTONEAREST;

  vx_orig = vx;
  abs_x = __riscv_vfabs (vx_orig, vlen);

  /* Special cases: NaN/Inf and overflow.  */
  VUINT vclass = __riscv_vfclass (vx_orig, vlen);
  VBOOL NaN_Inf;
  IDENTIFY_F (vclass, class_NaN | class_Inf, NaN_Inf, vlen);

  /* Overflow threshold for sinhf: log(2*FLT_MAX) ~= 89.415985.  */
  VBOOL overflow = __riscv_vmfgt
    (abs_x, RVV_LIBMVEC_HYPERBOLIC_OVERFLOW, vlen);
  /* Keep the near-zero result in the same prepared/sanitized lane class as
     the other exceptional paths.  This prevents the post-core repair from
     touching a NaN/Inf lane through vx_orig after its core value was
     sanitized.  */
  tiny_args = __riscv_vmflt
    (abs_x, RVV_SINHF_CLEAN_TINY_INPUT_F, vlen);

  special_args = __riscv_vmor
    (NaN_Inf, __riscv_vmor (overflow, tiny_args, vlen), vlen);
  nb_special_args = __riscv_vcpop (special_args, vlen);

  if (nb_special_args > 0)
    {
      /* Propagate NaNs/Infs via x*1, and use signed Inf for overflow.  */
      vy_special = __riscv_vfmul (vx_orig, fp_posOne_f, vlen);

      /* Materialize the scalar overflow/inexact side effect on finite
         overflow lanes before replacing the value with signed infinity.  */
      VFLOAT overflow_probe
          = __riscv_vfmul (overflow, abs_x, RVV_FLT_MAX_F,
                            vlen);
      vy_special = __riscv_vmerge
        (vy_special, overflow_probe, overflow, vlen);

      VFLOAT inf_sign
          = __riscv_vfsgnj (VFMV_VF_F (fp_posInf_f, vlen), vx_orig, vlen);
      vy_special = __riscv_vmerge (vy_special, inf_sign, overflow, vlen);

      /* For tiny finite x, sinhf(x) rounds to x; use x+x*2^-60 to retain
         the scalar inexact/underflow side effect.  */
      VFLOAT tiny_result
          = __riscv_vfmadd (tiny_args, vx_orig,
                            RVV_LIBMVEC_HYPERBOLIC_TINY_CORRECTION_F, vx_orig,
                            vlen);
      vy_special = __riscv_vmerge (vy_special, tiny_result, tiny_args, vlen);

      /* Avoid spurious exceptions from the core path (e.g. NaN->int).  */
      abs_x = __riscv_vfmerge (abs_x, fp_posZero_f, special_args, vlen);
    }

  /* Use a scaled exp path for large |x| so we don't overflow exp(|x|) when
     sinhf(|x|) is still finite.  */
  VBOOL use_scaled = __riscv_vmfge
    (abs_x, RVV_LIBMVEC_HYPERBOLIC_LARGE_INPUT, vlen);

  /* Small path: sinh(x) = 0.5 * (exp(|x|) - 1/exp(|x|)).  */
  VFLOAT ax_small = __riscv_vfmerge (abs_x, fp_posZero_f, use_scaled, vlen);
  VFLOAT d_hi, d_lo;
  RVV_LIBMVEC_IMPL_NAME (sinhf_exp_pair, LMUL)
    (ax_small, VFMV_VF_F (0.0f, vlen), &d_hi, &d_lo, vlen);
  /* The source DF reciprocal computes a low correction as t*u.  Once
     |x| > 70, exp(|x|) is large enough that this correction is discarded by
     the final float result but can itself underflow in the target FPU.  Keep
     the source graph for the useful lanes and quarantine that discarded
     correction through an exact neutral reciprocal input.  The resulting
     reciprocal pair is then zeroed for the quarantined lanes, which is
     mathematically invisible at sinhf's output precision and keeps the
     clean EH path from manufacturing FE_UNDERFLOW.  */
  VBOOL reciprocal_quarantine = __riscv_vmfgt
    (ax_small, RVV_SINHF_CLEAN_RECIPROCAL_QUARANTINE_F, vlen);
  VFLOAT reciprocal_d_hi
      = __riscv_vfmerge (d_hi, fp_posOne_f, reciprocal_quarantine, vlen);
  VFLOAT reciprocal_d_lo
      = __riscv_vfmerge (d_lo, fp_posZero_f, reciprocal_quarantine, vlen);
  VFLOAT e_hi, e_lo;
  RVV_SLEEF_ORIG_DFREC_VF2_VF2
    (reciprocal_d_hi, reciprocal_d_lo, e_hi, e_lo, vlen);
  e_hi = __riscv_vfmerge (e_hi, fp_posZero_f, reciprocal_quarantine, vlen);
  e_lo = __riscv_vfmerge (e_lo, fp_posZero_f, reciprocal_quarantine, vlen);
  VFLOAT diff_hi, diff_lo;
  RVV_SLEEF_ORIG_DFSUB_VF2_VF2
    (d_hi, d_lo, e_hi, e_lo, diff_hi, diff_lo, vlen);
  VFLOAT y_small
      = __riscv_vfmul (__riscv_vfadd (diff_hi, diff_lo, vlen),
                       RVV_HALF_F, vlen);

  /* Large path: sinh(x) ~ exp(|x|) / 2 = exp(|x| - ln2).  */
  VFLOAT y_large = VFMV_VF_F (fp_posZero_f, vlen);
  if (__riscv_vcpop (use_scaled, vlen) > 0)
    {
      VFLOAT ax_large
          = __riscv_vfmerge (abs_x, fp_posZero_f,
                             __riscv_vmnot (use_scaled, vlen), vlen);
      VFLOAT scaled_hi, scaled_lo;
      DFADD2_VF2_VF_S (ax_large, VFMV_VF_F (0.0f, vlen),
                       -RVV_LN2_SPLIT_HI_F, scaled_hi, scaled_lo, vlen);
      DFADD2_VF2_VF_S (scaled_hi, scaled_lo, -RVV_LN2_SPLIT_LO_F,
                       scaled_hi, scaled_lo, vlen);

      VFLOAT d_hi_s, d_lo_s;
      RVV_LIBMVEC_IMPL_NAME (sinhf_exp_pair, LMUL)
        (scaled_hi, scaled_lo, &d_hi_s, &d_lo_s, vlen);
      y_large = __riscv_vfadd (d_hi_s, d_lo_s, vlen);
    }

  vy = __riscv_vmerge (y_small, y_large, use_scaled, vlen);

  /* sinh is odd.  */
  vy = __riscv_vfsgnj (vy, vx_orig, vlen);

  if (nb_special_args > 0)
    vy = __riscv_vmerge (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef RVV_SINHF_CLEAN_TINY_INPUT_F
#undef RVV_SINHF_CLEAN_RECIPROCAL_QUARANTINE_F

#undef LMUL
#undef SIMDLEN
