/* Single-precision vector asinh function.
   asinh(x) = log(x + sqrt(x*x + 1)) for all finite x.

   Based on SLEEF xasinhf implementation.

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

/* The importer includes this skeleton once for every LMUL.  Keep the
   routine-local helpers in distinct namespaces so their vector parameter
   types cannot collide across those inclusions.  */
static inline void
RVV_ASINHF_CLEAN_LOGK2F
  (VFLOAT d_hi, VFLOAT d_lo, VFLOAT *r_hi, VFLOAT *r_lo,
   VBOOL active, size_t vlen)
{
  /* Extract exponent: e = ilogb(d_hi * (1/0.75)) */
  /* Inactive lanes are sanitized to d=1.  Keep the first potentially
     inexact normalization multiply masked so discarded special lanes do not
     alter the caller's fenv state.  The masked-off value d_hi=1 then follows
     an exact x=0/log(1) neutral path through the helper.  */
  VFLOAT scaled = VFMUL_VF_M (active, d_hi, 1.0f / 0.75f, vlen);
  VINT e = __riscv_vsra (F_AS_I (scaled), MAN_LEN, vlen);
  e = __riscv_vsub (e, EXP_BIAS, vlen);

  /* Scale (d_hi, d_lo) by 2^(-e) to get m in [0.75, 1.5) */
  VINT neg_e = __riscv_vneg (e, vlen);
  VINT scale_bits = __riscv_vsll (__riscv_vadd (neg_e, EXP_BIAS, vlen), MAN_LEN,
                                  vlen);
  VFLOAT scale = I_AS_F (scale_bits);

  VFLOAT m_hi, m_lo;
  DFSCALE_VF2_VF2_VF (d_hi, d_lo, scale, m_hi, m_lo, vlen);

  /* Compute x = (m - 1) / (m + 1) using double-float division */
  VFLOAT num_hi, num_lo;
  DFADD2_VF2_VF_S (m_hi, m_lo, -1.0f, num_hi, num_lo, vlen);

  VFLOAT den_hi, den_lo;
  DFADD2_VF2_VF_S (m_hi, m_lo, 1.0f, den_hi, den_lo, vlen);

  VFLOAT x_hi, x_lo;
  RVV_SLEEF_ORIG_DFDIV_VF2_VF2
    (num_hi, num_lo, den_hi, den_lo, x_hi, x_lo, vlen);

  /* x2 = x * x (double-float square) */
  VFLOAT x2_hi, x2_lo;
  RVV_SLEEF_ORIG_DFSQU_VF2
    (x_hi, x_lo, x2_hi, x2_lo, vlen);

  /* Polynomial: t = P0 + x2*(P1 + x2*(P2 + x2*P3)) */
  VFLOAT poly = VFMV_VF (RVV_ASINHF_LOG_P0, vlen);
  poly = __riscv_vfmadd
    (poly, x2_hi, VFMV_VF (RVV_ASINHF_LOG_P1, vlen), vlen);
  poly = __riscv_vfmadd
    (poly, x2_hi, VFMV_VF (RVV_ASINHF_LOG_P2, vlen), vlen);
  poly = __riscv_vfmadd
    (poly, x2_hi, VFMV_VF (RVV_ASINHF_LOG_P3, vlen), vlen);

  /* s = e * log(2) as double-float */
  VFLOAT e_flt = __riscv_vfcvt_f (e, vlen);
  VFLOAT s_hi, s_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (VFMV_VF (RVV_ASINHF_L2U, vlen),
     VFMV_VF (RVV_ASINHF_L2L, vlen), e_flt,
     s_hi, s_lo, vlen);

  /* s = s + 2*x (double-float addition) */
  VFLOAT two_x_hi, two_x_lo;
  DFSCALE_VF2_VF2_S (x_hi, x_lo, 2.0f, two_x_hi, two_x_lo, vlen);

  VFLOAT sum_hi, sum_lo;
  RVV_SLEEF_ORIG_DFADD_VF2_VF2
    (s_hi, s_lo, two_x_hi, two_x_lo, sum_hi, sum_lo, vlen);

  /* s = s + x2 * x * poly */
  VFLOAT x2x_hi, x2x_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (x2_hi, x2_lo, x_hi, x_lo, x2x_hi, x2x_lo, vlen);

  VFLOAT term_hi, term_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (x2x_hi, x2x_lo, poly, term_hi, term_lo, vlen);

  RVV_SLEEF_ORIG_DFADD_VF2_VF2
    (sum_hi, sum_lo, term_hi, term_lo, *r_hi, *r_lo, vlen);
}

static inline void
RVV_ASINHF_CLEAN_LOGF_LARGE
  (VFLOAT x, VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  /* Extract exponent and mantissa without overflowing the scale step. */
  VINT x_bits = F_AS_I (x);
  VINT e = __riscv_vsra (x_bits, MAN_LEN, vlen);
  e = __riscv_vsub (e, EXP_BIAS, vlen);

  VINT mant_mask = VMVI_VX ((1 << MAN_LEN) - 1, vlen);
  VINT mant_bits = __riscv_vand (x_bits, mant_mask, vlen);
  VINT exp_one = VMVI_VX (EXP_BIAS << MAN_LEN, vlen);
  VFLOAT m_hi = I_AS_F (__riscv_vor (mant_bits, exp_one, vlen));
  VFLOAT m_lo = VFMV_VF (0.0f, vlen);

  /* Keep m in [0.75, 1.5) for polynomial accuracy. */
  VINT half_bit = VMVI_VX (1 << (MAN_LEN - 1), vlen);
  VBOOL need_adjust = __riscv_vmsne (__riscv_vand (mant_bits, half_bit, vlen),
                                     0, vlen);
  VFLOAT m_half = __riscv_vfmul (m_hi, 0.5f, vlen);
  m_hi = __riscv_vmerge (m_hi, m_half, need_adjust, vlen);

  VINT e_plus1 = __riscv_vadd (e, 1, vlen);
  e = __riscv_vmerge (e, e_plus1, need_adjust, vlen);
  /* Return log(2*x) for the large-x asinh path. */
  e = __riscv_vadd (e, 1, vlen);

  /* Compute x = (m - 1) / (m + 1) using double-float division. */
  VFLOAT num_hi, num_lo;
  DFADD2_VF2_VF_S (m_hi, m_lo, -1.0f, num_hi, num_lo, vlen);

  VFLOAT den_hi, den_lo;
  DFADD2_VF2_VF_S (m_hi, m_lo, 1.0f, den_hi, den_lo, vlen);

  VFLOAT x_hi, x_lo;
  DFDIV_VF2_VF2_VF2 (num_hi, num_lo, den_hi, den_lo, x_hi, x_lo, vlen);

  VFLOAT x2_hi, x2_lo;
  DFSQU_VF2_VF2 (x_hi, x_lo, x2_hi, x2_lo, vlen);

  VFLOAT poly = VFMV_VF (RVV_ASINHF_LOG_P0, vlen);
  poly = __riscv_vfmadd
    (poly, x2_hi, VFMV_VF (RVV_ASINHF_LOG_P1, vlen), vlen);
  poly = __riscv_vfmadd
    (poly, x2_hi, VFMV_VF (RVV_ASINHF_LOG_P2, vlen), vlen);
  poly = __riscv_vfmadd
    (poly, x2_hi, VFMV_VF (RVV_ASINHF_LOG_P3, vlen), vlen);

  VFLOAT e_flt = __riscv_vfcvt_f (e, vlen);
  VFLOAT s_hi = __riscv_vfmul (e_flt, RVV_ASINHF_L2U, vlen);
  VFLOAT s_lo = __riscv_vfmul (e_flt, RVV_ASINHF_L2L, vlen);

  VFLOAT two_x_hi, two_x_lo;
  DFSCALE_VF2_VF2_S (x_hi, x_lo, 2.0f, two_x_hi, two_x_lo, vlen);

  VFLOAT sum_hi, sum_lo;
  DFADD_VF2_VF2_VF2 (s_hi, s_lo, two_x_hi, two_x_lo, sum_hi, sum_lo, vlen);

  VFLOAT x2x_hi, x2x_lo;
  DFMUL_VF2_VF2_VF2 (x2_hi, x2_lo, x_hi, x_lo, x2x_hi, x2x_lo, vlen);

  VFLOAT term_hi, term_lo;
  DFMUL_VF2_VF2_VF (x2x_hi, x2x_lo, poly, term_hi, term_lo, vlen);

  DFADD_VF2_VF2_VF2 (sum_hi, sum_lo, term_hi, term_lo, *r_hi, *r_lo, vlen);
}

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT vx, vy, vy_special;
  VFLOAT ax = __riscv_vfabs (x, vlen);
  VBOOL special_args;

  SET_ROUNDTONEAREST;
  vx = x;

  VUINT vclass = __riscv_vfclass (vx, vlen);
  VBOOL x_inf, x_nan, x_snan, x_subnormal, x_zero;
  IDENTIFY_F (vclass, class_Inf, x_inf, vlen);
  IDENTIFY_F (vclass, class_NaN, x_nan, vlen);
  IDENTIFY_F (vclass, class_sNaN, x_snan, vlen);
  IDENTIFY_F (vclass, class_Subnormal, x_subnormal, vlen);
  IDENTIFY_F (vclass, class_Zero, x_zero, vlen);

  special_args = __riscv_vmor (x_zero, x_subnormal, vlen);
  special_args = __riscv_vmor (special_args, x_inf, vlen);
  special_args = __riscv_vmor (special_args, x_nan, vlen);

  vy_special = vx;
  if (__riscv_vcpop (special_args, vlen) > 0)
    {
      /* Propagate NaNs via x + x so sNaN quieting stays lane-local.  */
      VFLOAT vy_nan = __riscv_vfadd (x_nan, vx, vx, vlen);
      /* Some RVV masked implementations do not expose the sNaN trap from
         the original payload consistently.  Re-issue it with an explicit
         signaling operand in an otherwise exact-safe vector, without
         touching qNaN lanes.  */
      VFLOAT s_nan = VFMV_VF (fp_sNaN, vlen);
      VFLOAT s_nan_input = __riscv_vmerge
        (VFMV_VF (0.0f, vlen), s_nan, x_snan, vlen);
      VFLOAT vy_snan = __riscv_vfadd (s_nan_input, s_nan_input, vlen);
      vy_nan = __riscv_vmerge (vy_nan, vy_snan, x_snan, vlen);
      vy_special = __riscv_vmerge (vy_special, vy_nan, x_nan, vlen);

      /* Match scalar asinhf for infinities with x + x.  */
      VFLOAT vy_inf = __riscv_vfadd (x_inf, vx, vx, vlen);
      vy_special = __riscv_vmerge (vy_special, vy_inf, x_inf, vlen);

    }

  /* No finite core lane remains for an all-zero/Inf/NaN vector.  Return the
     already prepared results before normalization, reciprocal, or log work
     can manufacture flags on discarded lanes.  Subnormal-only vectors are
     deliberately not included: their scalar contract requires the core to
     raise underflow/inexact.  */
  VBOOL no_core_args = __riscv_vmor (x_zero, x_inf, vlen);
  no_core_args = __riscv_vmor (no_core_args, x_nan, vlen);
  if (__riscv_vcpop (no_core_args, vlen) == vlen)
    {
      RESTORE_FRM;
      return vy_special;
    }

  /* Materialize the large finite result before the ordinary expression can
     overflow in sqrt(d*d + 1) or d + x.  Infinite and NaN lanes are excluded
     from this class and receive their source-compatible special result above.
     Use a benign finite operand for all other lanes while evaluating the
     helper, so discarded lanes cannot raise an exception.  */
  VBOOL nonfinite = __riscv_vmor (x_inf, x_nan, vlen);
  VBOOL finite = __riscv_vmnot (nonfinite, vlen);
  /* Classify NaN/Inf lanes before the signaling magnitude comparison.  A
     later mask cannot undo an exception already raised by comparing qNaN. */
  VFLOAT ax_compare = __riscv_vfmerge
    (ax, fp_posZero, nonfinite, vlen);
  VBOOL x_large = __riscv_vmfgt
    (ax_compare, RVV_ASINHF_SQRT_FLT_MAX, vlen);
  x_large = __riscv_vmand (x_large, finite, vlen);

  /* 0.5 is an exact neutral input for the normalized log helper: it maps to
     m=1 and e=0, so discarded lanes do not manufacture FE_INEXACT while the
     selected large lanes retain their original argument.  */
  VFLOAT large_arg = __riscv_vfmerge
    (ax, 0.5f, __riscv_vmnot (x_large, vlen), vlen);
  VFLOAT large_hi, large_lo;
  RVV_ASINHF_CLEAN_LOGF_LARGE (large_arg, &large_hi, &large_lo, vlen);
  VFLOAT large_val = __riscv_vfadd (large_hi, large_lo, vlen);
  large_val = __riscv_vfsgnj (large_val, x, vlen);
  vy_special = __riscv_vmerge (vy_special, large_val, x_large, vlen);

  /* Feed benign zeros to the finite core for every lane whose result is
     already prepared.  This is the same classify -> prepare -> sanitize ->
     core -> merge ordering used by the VECLIBM EH path.  */
  VBOOL core_quarantine = __riscv_vmor (special_args, x_large, vlen);
  vx = __riscv_vfmerge (vx, fp_posZero, core_quarantine, vlen);
  /* GCC may otherwise move an exceptional finite-core operation across this
     classify/prepare/sanitize boundary under -O2.  Keep discarded special
     lanes sanitized before the core, matching the VECLIBM EH contract.  */
  asm volatile ("" ::: "memory");

  VFLOAT ay = __riscv_vfabs (vx, vlen);
  VBOOL use_recip = __riscv_vmfgt (ay, 1.0f, vlen);

  /* d = |x| (lo = 0) or 1/x for |x| > 1 */
  VFLOAT d_hi = ay;
  VFLOAT d_lo = VFMV_VF (0.0f, vlen);

  VFLOAT rec_hi, rec_lo;
  VFLOAT rec_input = __riscv_vfmerge (vx, fp_posOne,
                                      __riscv_vmnot (use_recip, vlen), vlen);
  RVV_SLEEF_ORIG_DFREC_VF2_VF (rec_input, rec_hi, rec_lo, vlen);

  d_hi = __riscv_vmerge (d_hi, rec_hi, use_recip, vlen);
  d_lo = __riscv_vmerge (d_lo, rec_lo, use_recip, vlen);

  /* d = sqrt(d*d + 1) */
  VFLOAT d2_hi, d2_lo;
  RVV_SLEEF_ORIG_DFSQU_VF2
    (d_hi, d_lo, d2_hi, d2_lo, vlen);

  VFLOAT sum_hi, sum_lo;
  RVV_SLEEF_ORIG_DFADD2_VF2_VF_S
    (d2_hi, d2_lo, 1.0f, sum_hi, sum_lo, vlen);

  VFLOAT sqrt_hi, sqrt_lo;
  RVV_SLEEF_ORIG_DFSQRT_VF2_VF2
    (sum_hi, sum_lo, sqrt_hi, sqrt_lo, vlen);

  /* If |x| > 1, d = sqrt(d*d + 1) * |x| */
  VFLOAT mul_hi, mul_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (sqrt_hi, sqrt_lo, ay, mul_hi, mul_lo, vlen);

  d_hi = __riscv_vmerge (sqrt_hi, mul_hi, use_recip, vlen);
  d_lo = __riscv_vmerge (sqrt_lo, mul_lo, use_recip, vlen);

  /* logk2f(normalize(d + x)) */
  VFLOAT arg_hi, arg_lo;
  RVV_SLEEF_ORIG_DFADD2_VF2_VF
    (d_hi, d_lo, vx, arg_hi, arg_lo, vlen);

  VFLOAT norm_hi, norm_lo;
  RVV_SLEEF_ORIG_DFNORMALIZE_VF2
    (arg_hi, arg_lo, norm_hi, norm_lo, vlen);

  VFLOAT log_hi, log_lo;
  VBOOL core_args = __riscv_vmnot (core_quarantine, vlen);
  RVV_ASINHF_CLEAN_LOGK2F
    (norm_hi, norm_lo, &log_hi, &log_lo, core_args, vlen);

  vy = __riscv_vfadd (log_hi, log_lo, vlen);

  /* If the main path produced NaN, return signed Inf.  */
  VBOOL y_nan = __riscv_vmfne (vy, vy, vlen);
  VFLOAT inf_signed = __riscv_vfsgnj (VFMV_VF (fp_posInf, vlen), vx, vlen);
  vy = __riscv_vmerge (vy, inf_signed, y_nan, vlen);

  VBOOL merge_args = __riscv_vmor (special_args, x_large, vlen);
  vy = __riscv_vmerge (vy, vy_special, merge_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
