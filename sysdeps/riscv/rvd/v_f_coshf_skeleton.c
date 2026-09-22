/* Single-precision vector coshf function.
   Based on SLEEF xcoshf (u10 = 1.0 ULP precision).

   Copyright (C) 2026 Free Software Foundation, Inc.
   This file is part of the GNU C Library. */

#include <vecmath_config.h>
#include <v_math.h>

#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

/* ============================================================================
   Double-float helper macros for RVV (vfloat2 operations)
   ============================================================================ */

/* Shared float-float helpers live in v_math.h.  */

/* dfsqu: (a_hi, a_lo)^2 -> (r_hi, r_lo) */

/* Keep the source-shaped power-of-two construction local to this clean
   translation unit.  The two floating-point multiplies are required for a
   double-float value: exponent-bit injection is only valid for a normalized
   high part and corrupts zero, subnormal, or negative low parts.  */
static inline VFLOAT
RVV_COSHF_SOURCE_MANGLE (coshf_source_vpow2i) (VINT e, size_t vlen)
{
  VINT biased = __riscv_vadd (e, EXP_BIAS_F, vlen);
  return I_AS_F (__riscv_vsll (biased, MAN_LEN_F, vlen));
}

static inline VFLOAT
RVV_COSHF_SOURCE_MANGLE (coshf_source_vldexp2) (VFLOAT d, VINT e, size_t vlen)
{
  VINT e_half = __riscv_vsra (e, 1, vlen);
  VINT e_rest = __riscv_vsub (e, e_half, vlen);
  VFLOAT p_half = RVV_COSHF_SOURCE_MANGLE (coshf_source_vpow2i)
    (e_half, vlen);
  VFLOAT p_rest = RVV_COSHF_SOURCE_MANGLE (coshf_source_vpow2i)
    (e_rest, vlen);
  return VFMUL_VV (VFMUL_VV (d, p_half, vlen), p_rest, vlen);
}

/* expk2f: exp((d_hi, d_lo)) -> (r_hi, r_lo)
   SLEEF expk2f helper for extended precision exponential */
static inline void
RVV_COSHF_EXPK2F_VF2 (VFLOAT d_hi, VFLOAT d_lo, VFLOAT *r_hi, VFLOAT *r_lo,
                      size_t vlen)
{
  /* Argument reduction: exp(d) = 2^q * exp(s) where |s| <= ln(2)/2 */
  VFLOAT u = VFMUL_VF (VFADD_VV (d_hi, d_lo, vlen),
                       RVV_SLEEF_EXPM1F_R_LN2,
                       vlen);
  VINT q = __riscv_vfcvt_x (u, vlen);
  VFLOAT q_flt = __riscv_vfcvt_f (q, vlen);

  /* s = d - q * ln(2) using double-float arithmetic */
  VFLOAT s_hi, s_lo;
  VFLOAT tmp = VFMUL_VF (q_flt, -RVV_SLEEF_EXPM1F_L2U, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF
    (d_hi, d_lo, tmp, s_hi, s_lo, vlen);

  tmp = VFMUL_VF (q_flt, -RVV_SLEEF_EXPM1F_L2L, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF
    (s_hi, s_lo, tmp, s_hi, s_lo, vlen);

  /* Polynomial approximation for exp(s) - 1 */
  u = VFMV_VF (RVV_SLEEF_EXPM1F_POLY0, vlen);
  u = VFMADD_VV (u, s_hi, VFMV_VF (RVV_SLEEF_EXPM1F_POLY1, vlen), vlen);
  u = VFMADD_VV (u, s_hi, VFMV_VF (RVV_SLEEF_EXPM1F_POLY2, vlen), vlen);
  u = VFMADD_VV (u, s_hi, VFMV_VF (RVV_SLEEF_EXPM1F_POLY3, vlen), vlen);

  /* t = s * u + 0.166666659414234244790680580464 */
  VFLOAT t_hi, t_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (s_hi, s_lo, u, t_hi, t_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF_S
    (t_hi, t_lo, RVV_SLEEF_EXPM1F_TAYLOR_C,
     t_hi, t_lo, vlen);

  /* t = s * t + 0.5 */
  VFLOAT tmp_hi, tmp_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (s_hi, s_lo, t_hi, t_lo, tmp_hi, tmp_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF_S
    (tmp_hi, tmp_lo, RVV_SLEEF_EXPM1F_HALF, t_hi, t_lo, vlen);

  /* t = s + s^2 * t */
  VFLOAT s2_hi, s2_lo;
  RVV_SLEEF_ORIG_DFSQU_VF2
    (s_hi, s_lo, s2_hi, s2_lo, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (s2_hi, s2_lo, t_hi, t_lo, tmp_hi, tmp_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF2
    (s_hi, s_lo, tmp_hi, tmp_lo, t_hi, t_lo, vlen);

  /* t = 1 + t */
  RVV_SLEEF_ORIG_DFADD_VF_VF2
    (VFMV_VF (1.0f, vlen), t_hi, t_lo, t_hi, t_lo, vlen);

  /* Scale by 2^q using the source's two-multiply vldexp2 operation. */
  t_hi = RVV_COSHF_SOURCE_MANGLE (coshf_source_vldexp2) (t_hi, q, vlen);
  t_lo = RVV_COSHF_SOURCE_MANGLE (coshf_source_vldexp2) (t_lo, q, vlen);

  /* Handle underflow: if d.x < -104, return 0 */
  VBOOL underflow = __riscv_vmflt (d_hi, RVV_SLEEF_EXPM1F_UNDERFLOW, vlen);
  t_hi = __riscv_vfmerge (t_hi, fp_posZero_f, underflow, vlen);
  t_lo = __riscv_vfmerge (t_lo, fp_posZero_f, underflow, vlen);

  *r_hi = t_hi;
  *r_lo = t_lo;
}

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT vx = x;
  VFLOAT y, vy_special;
  VBOOL special_args;
  UINT nb_special_args;

  SET_ROUNDTONEAREST;

  /* y = |x| */
  y = __riscv_vfabs (vx, vlen);

  /* Special cases: NaN/Inf and overflow.  */
  VUINT vclass = __riscv_vfclass (vx, vlen);
  VBOOL NaN_Inf, inf_args;
  IDENTIFY_F (vclass, class_NaN | class_Inf, NaN_Inf, vlen);
  IDENTIFY_F (vclass, class_Inf, inf_args, vlen);

  /* Overflow threshold for coshf: log(2*FLT_MAX) ~= 89.415985.  Keep raw
     NaN/Inf lanes out of the ordered comparison so the compare itself does
     not become an EH source.  */
  VFLOAT y_for_compare
      = __riscv_vfmerge (y, fp_posZero_f, NaN_Inf, vlen);
  VBOOL overflow = __riscv_vmfgt
    (y_for_compare, RVV_LIBMVEC_COSHF_OVERFLOW, vlen);

  special_args = __riscv_vmor (NaN_Inf, overflow, vlen);
  nb_special_args = __riscv_vcpop (special_args, vlen);
  if (nb_special_args > 0)
    {
      /* Propagate NaNs via x*1; both input infinities and finite overflow
         produce the positive infinity required by the even cosh function.  */
      vy_special = __riscv_vfmul (vx, fp_posOne_f, vlen);
      VBOOL inf_or_overflow = __riscv_vmor (inf_args, overflow, vlen);
      vy_special = __riscv_vfmerge
        (vy_special, fp_posInf_f, inf_or_overflow, vlen);

      /* Keep special lanes out of exp/log/int-convert core work.  */
      y = __riscv_vfmerge (y, fp_posZero_f, special_args, vlen);
    }

  /* Use exp(|x| - ln2) for large |x| to avoid overflow. */
  VBOOL large_x = __riscv_vmfge
    (y, RVV_LIBMVEC_COSHF_LARGE_INPUT, vlen);

  /* The scalar coshf contract treats sufficiently small nonzero inputs as
     an inexact result of one, without exposing the underflow generated by
     squaring the input in the vector expk2f core.  Quarantine those lanes
     from the core and raise only the corresponding inexact below.  */
  VBOOL tiny_x = __riscv_vmand
    (__riscv_vmfne (y, fp_posZero_f, vlen),
     __riscv_vmflt (y, RVV_LIBMVEC_COSHF_TINY_INPUT, vlen), vlen);

  /* d = expk2f(|x|) for lanes that will use the small path.  Large,
     exceptional, and tiny lanes are replaced before entering the finite
     core. */
  VBOOL core_quarantine = __riscv_vmor (large_x, tiny_x, vlen);
  VFLOAT small_input
    = __riscv_vfmerge (y, fp_posZero_f, core_quarantine, vlen);
  VFLOAT d_hi, d_lo;
  RVV_COSHF_EXPK2F_VF2
    (small_input, VFMV_VF (fp_posZero_f, vlen), &d_hi, &d_lo, vlen);

  /* e = 1/d */
  VFLOAT e_hi, e_lo;
  DFREC_VF2_VF2 (d_hi, d_lo, e_hi, e_lo, vlen);

  /* d = d + e, y = 0.5 * (d_hi + d_lo) */
  VFLOAT sum_hi, sum_lo;
  RVV_SLEEF_ORIG_DFADD_VF2_VF2
    (d_hi, d_lo, e_hi, e_lo, sum_hi, sum_lo, vlen);
  VFLOAT y_small = __riscv_vfmul
    (__riscv_vfadd (sum_hi, sum_lo, vlen), RVV_SLEEF_EXPM1F_HALF, vlen);

  /* 1 + min-subnormal is an inexact-only operation under RNE.  Use it for
     quarantined tiny lanes so the public EH result matches scalar coshf
     without allowing the discarded x^2 path to set FE_UNDERFLOW.  */
  VFLOAT min_subnormal = I_AS_F (VMVI_VX (1, vlen));
  VFLOAT tiny_operand
    = VMERGE_VVM_F (VFMV_VF (fp_posZero_f, vlen), min_subnormal,
                     tiny_x, vlen);
  VFLOAT tiny_result = __riscv_vfadd
    (VFMV_VF (fp_posOne_f, vlen), tiny_operand, vlen);
  y_small = VMERGE_VVM_F (y_small, tiny_result, tiny_x, vlen);

  /* y_large = 2*expk2f(|x| - 2*ln2) for large inputs.  This is the same
     0.5*exp(|x|) value as exp(|x|-ln2), but its argument stays below the
     overflow boundary so a discarded candidate cannot set sticky overflow. */
  VBOOL non_large = __riscv_vmnot (large_x, vlen);
  VFLOAT y_large_input_hi = VFMERGE_VFM
    (y, RVV_LIBMVEC_COSHF_LARGE_NEUTRAL_HI, non_large, vlen);
  VFLOAT y_large_input_lo = VFMERGE_VFM
    (VFMV_VF (fp_posZero_f, vlen), RVV_LIBMVEC_COSHF_LARGE_NEUTRAL_LO,
     non_large, vlen);
  VFLOAT y_scaled_hi, y_scaled_lo;
  RVV_SLEEF_ORIG_DFADD2_VF2_VF
    (y_large_input_hi, y_large_input_lo,
     VFMV_VF (-RVV_SLEEF_COSHF_TWO_L2U, vlen),
     y_scaled_hi, y_scaled_lo, vlen);
  y_scaled_hi = VFMERGE_VFM
    (y_scaled_hi, RVV_LIBMVEC_COSHF_LARGE_NEUTRAL_LO, non_large, vlen);
  y_scaled_lo = VFMERGE_VFM
    (y_scaled_lo, fp_posZero_f, non_large, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF
    (y_scaled_hi, y_scaled_lo, VFMV_VF (-RVV_SLEEF_COSHF_TWO_L2L, vlen),
     y_scaled_hi, y_scaled_lo, vlen);
  VFLOAT d_hi_scaled, d_lo_scaled;
  RVV_COSHF_EXPK2F_VF2
    (y_scaled_hi, y_scaled_lo, &d_hi_scaled, &d_lo_scaled, vlen);
  VFLOAT y_large = __riscv_vfmul
    (__riscv_vfadd (d_hi_scaled, d_lo_scaled, vlen),
     RVV_SLEEF_EXPM1F_TWO, vlen);

  y = __riscv_vmerge (y_small, y_large, large_x, vlen);

  /* Keep finite overflow flags visible even though the unsafe input lanes
     were quarantined from the finite core.  The probe is discarded by the
     final special merge below. */
  if (nb_special_args > 0)
    {
      VFLOAT overflow_probe = VFMV_VF (fp_posZero_f, vlen);
      overflow_probe = VMERGE_VVM_F
        (overflow_probe, VFMV_VF (RVV_SLEEF_EXPM1F_MAX_FINITE, vlen),
         overflow, vlen);
      overflow_probe = VFMUL_VF
        (overflow_probe, RVV_SLEEF_EXPM1F_TWO, vlen);
      y = VMERGE_VVM_F (y, overflow_probe, overflow, vlen);
    }

  if (nb_special_args > 0)
    y = __riscv_vmerge (y, vy_special, special_args, vlen);

  RESTORE_FRM;
  return y;
}

#undef LMUL
#undef SIMDLEN
