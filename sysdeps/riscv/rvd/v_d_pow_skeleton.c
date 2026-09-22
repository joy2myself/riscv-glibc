/* Double-precision vector pow function.
   pow(x, y) computes x raised to the power y.

   This implementation follows the core algorithm of glibc scalar pow:
     pow(x, y) = exp (y * log (x)),
   using a log kernel that returns (hi, tail) and an exp kernel that
   consumes (x, xtail) to reduce error amplification for large |y|.

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

/* Pow lookup tables are owned by __v_pow_data.  */
static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (pow_log_inline, LMUL)
  (const struct v_pow_data *data, VFLOAT x, VFLOAT *tail, size_t vlen)
{
  /* Normalize subnormals: ix = asuint64 (x * 0x1p52) - 52<<52.  */
  VUINT ix = F_AS_U (x);
  VUINT top12 = __riscv_vsrl (ix, 52, vlen);
  VBOOL top12_is_zero = __riscv_vmseq (top12, 0, vlen);
  if (__riscv_vcpop (top12_is_zero, vlen) > 0)
    {
      VFLOAT x_scaled = VFMUL_VF_M (top12_is_zero, x, 0x1.0p52, vlen);
      VUINT ix_scaled = F_AS_U (x_scaled);
      ix_scaled = __riscv_vand (ix_scaled, 0x7fffffffffffffffULL, vlen);
      ix_scaled = __riscv_vsub (ix_scaled, 52ULL << 52, vlen);
      ix = __riscv_vmerge (ix, ix_scaled, top12_is_zero, vlen);
    }

  /* x = 2^k z; where z is in range [OFF,2*OFF) and exact.  */
  VUINT tmp = __riscv_vsub (ix, RVV_POW_LOG_OFF, vlen);
  VUINT i = __riscv_vand (
      __riscv_vsrl (tmp, 52 - V_POW_LOG_TABLE_BITS, vlen),
      V_POW_LOG_N - 1, vlen);
  VINT k = __riscv_vsra (U_AS_I (tmp), 52, vlen);
  VUINT iz = __riscv_vsub (ix, __riscv_vand (tmp, 0xfffULL << 52, vlen), vlen);
  VFLOAT z = U_AS_F (iz);
  VFLOAT kd = __riscv_vfcvt_f (k, vlen);

  VUINT i_off = __riscv_vsll (i, 3, vlen);
  VFLOAT invc = __riscv_vluxei64 (data->log_invc, i_off, vlen);
  VFLOAT logc = __riscv_vluxei64 (data->log_logc, i_off, vlen);
  VFLOAT logctail = __riscv_vluxei64 (data->log_logctail, i_off, vlen);

  /* r = fma (z, invc, -1.0).  */
  VFLOAT r = VFMV_VF (-1.0, vlen);
  r = VFMACC_VV (r, z, invc, vlen);

  /* k*ln2 + log(c) + r.  */
  VFLOAT t1 = __riscv_vfmadd (kd, RVV_POW_LOG_LN2_HI_D, logc, vlen);
  VFLOAT t2 = __riscv_vfadd (t1, r, vlen);
  VFLOAT lo1 = __riscv_vfmadd (kd, RVV_POW_LOG_LN2_LO_D, logctail, vlen);
  VFLOAT lo2 = __riscv_vfadd (__riscv_vfsub (t1, t2, vlen), r, vlen);

  /* Polynomial evaluation.  */
  VFLOAT ar = __riscv_vfmul (r, RVV_POW_LOG_POLY0_D, vlen);
  VFLOAT ar2 = __riscv_vfmul (r, ar, vlen);
  VFLOAT ar3 = __riscv_vfmul (r, ar2, vlen);

  VFLOAT hi = __riscv_vfadd (t2, ar2, vlen);
  VFLOAT lo3 = __riscv_vfmsub (ar, r, ar2, vlen);
  VFLOAT lo4 = __riscv_vfadd (__riscv_vfsub (t2, hi, vlen), ar2, vlen);

  VFLOAT t = __riscv_vfmadd (r, RVV_POW_LOG_POLY6_D,
                             VFMV_VF (RVV_POW_LOG_POLY5_D, vlen),
                             vlen);
  t = __riscv_vfmadd (ar2, t,
                      __riscv_vfmadd (r, RVV_POW_LOG_POLY4_D,
                                      VFMV_VF (RVV_POW_LOG_POLY3_D, vlen), vlen),
                      vlen);
  t = __riscv_vfmadd (ar2, t,
                      __riscv_vfmadd (r, RVV_POW_LOG_POLY2_D,
                                      VFMV_VF (RVV_POW_LOG_POLY1_D, vlen), vlen),
                      vlen);
  VFLOAT p = __riscv_vfmul (ar3, t, vlen);

  VFLOAT lo = __riscv_vfadd (lo1, lo2, vlen);
  lo = __riscv_vfadd (lo, lo3, vlen);
  lo = __riscv_vfadd (lo, lo4, vlen);
  lo = __riscv_vfadd (lo, p, vlen);

  VFLOAT y = __riscv_vfadd (hi, lo, vlen);
  *tail = __riscv_vfadd (__riscv_vfsub (hi, y, vlen), lo, vlen);
  return y;
}

static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (pow_exp_specialcase, LMUL)
  (VFLOAT tmp, VUINT sbits, VINT ki, VBOOL special, size_t vlen)
{
  /* k > 0 may overflow; k < 0 may underflow in scale*(1+tmp).  */
  VBOOL k_pos = __riscv_vmsge (ki, 0, vlen);

  VUINT sbits_pos = __riscv_vsub (sbits, 1009ULL << 52, vlen);
  VFLOAT scale_pos = U_AS_F (sbits_pos);
  VFLOAT y_pos
      = __riscv_vfmacc (special, scale_pos, scale_pos, tmp, vlen);
  y_pos = __riscv_vfmul (special, y_pos, 0x1.0p1009, vlen);

  VUINT sbits_neg = __riscv_vadd (sbits, 1022ULL << 52, vlen);
  VFLOAT scale_neg = U_AS_F (sbits_neg);
  VFLOAT y_neg
      = __riscv_vfmacc (special, scale_neg, scale_neg, tmp, vlen);
  y_neg = __riscv_vfmul (special, y_neg, 0x1.0p-1022, vlen);

  return VMERGE_VVM_F (y_neg, y_pos, k_pos, vlen);
}

static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (pow_exp_inline, LMUL)
  (const struct v_pow_data *data, VFLOAT x, VFLOAT xtail, VUINT sign_bias,
   size_t vlen)
{
  VUINT top12 = __riscv_vsrl (F_AS_U (x), 52, vlen);
  VUINT abstop = __riscv_vand (top12, 0x7ff, vlen);

  /* |x| is tiny (<= 2^-54) or large (>= 512).  */
  VUINT diff = __riscv_vsub (abstop, RVV_POW_EXP_SMALL_TOP, vlen);
  VBOOL uoflow = __riscv_vmsgeu (diff, RVV_POW_EXP_THRES_TOP, vlen);
  UINT nb_uoflow = __riscv_vcpop (uoflow, vlen);

  VBOOL tiny = __riscv_vmsltu (abstop, RVV_POW_EXP_SMALL_TOP, vlen);
  tiny = __riscv_vmand (uoflow, tiny, vlen);

  VBOOL huge = __riscv_vmsgeu (abstop, RVV_POW_EXP_HUGE_TOP, vlen);
  huge = __riscv_vmand (uoflow, huge, vlen);

  VBOOL special = __riscv_vmand (
      uoflow, __riscv_vmnot (__riscv_vmor (tiny, huge, vlen), vlen), vlen);

  /* Mask out lanes we will override, to keep conversions in-range.  */
  VBOOL not_core = __riscv_vmor (tiny, huge, vlen);
  VFLOAT x_core = VFMERGE_VFM (x, fp_posZero, not_core, vlen);
  VFLOAT xtail_core = VFMERGE_VFM (xtail, fp_posZero, not_core, vlen);

  /* Argument reduction: exp(x) = 2^(k/N) * exp(r), r in [-ln2/2N, ln2/2N].  */
  VFLOAT z = __riscv_vfmul (x_core, RVV_POW_EXP_INV_LN2_N_D, vlen);
  VINT ki = __riscv_vfcvt_x (z, vlen);
  VFLOAT kd = __riscv_vfcvt_f (ki, vlen);

  VFLOAT r = x_core;
  r = VFMACC_VF (r, RVV_POW_EXP_NEG_LN2_HI_N_D, kd, vlen);
  r = VFMACC_VF (r, RVV_POW_EXP_NEG_LN2_LO_N_D, kd, vlen);
  r = __riscv_vfadd (r, xtail_core, vlen);

  /* Table lookup for scale and tail.  */
  VUINT kiu = I_AS_U (ki);
  VUINT idx = __riscv_vand (kiu, V_POW_EXP_N - 1, vlen);
  idx = __riscv_vsll (idx, 1, vlen);
  VUINT off0 = __riscv_vsll (idx, 3, vlen);
  VUINT off1 = __riscv_vadd (off0, 8, vlen);
  VUINT tail_bits = __riscv_vluxei64 (data->exp_tab, off0, vlen);
  VUINT sbits_base = __riscv_vluxei64 (data->exp_tab, off1, vlen);

  VUINT top = __riscv_vsll (__riscv_vadd (kiu, sign_bias, vlen),
                            52 - V_POW_EXP_TABLE_BITS, vlen);
  VUINT sbits = __riscv_vadd (sbits_base, top, vlen);
  VFLOAT scale = U_AS_F (sbits);
  VFLOAT tail = U_AS_F (tail_bits);

  /* Polynomial approximation for exp(r) - 1.  */
  VFLOAT r2 = __riscv_vfmul (r, r, vlen);
  VFLOAT t = VFMACC_VF (VFMV_VF (RVV_POW_EXP_C2_D, vlen),
                         RVV_POW_EXP_C3_D, r, vlen);
  VFLOAT tmp = VFMACC_VV (r, r2, t, vlen);
  VFLOAT r4 = __riscv_vfmul (r2, r2, vlen);
  VFLOAT t2 = VFMACC_VF (VFMV_VF (RVV_POW_EXP_C4_D, vlen),
                         RVV_POW_EXP_C5_D, r, vlen);
  tmp = VFMACC_VV (tmp, r4, t2, vlen);
  tmp = __riscv_vfadd (tmp, tail, vlen);

  /* exp(x) ~= scale + scale * tmp.  */
  VFLOAT y = VFMACC_VV (scale, scale, tmp, vlen);

  if (__riscv_vcpop (special, vlen) > 0)
    {
      VFLOAT y_special
          = RVV_LIBMVEC_IMPL_NAME (pow_exp_specialcase, LMUL)
            (tmp, sbits, ki, special, vlen);
      y = VMERGE_VVM_F (y, y_special, special, vlen);
    }

  if (nb_uoflow > 0)
    {
      VUINT sign_mask = __riscv_vsll (sign_bias, 52 - V_POW_EXP_TABLE_BITS, vlen);

      /* Tiny: exp(x+xtail) ~= 1 + (x+xtail).  */
      VFLOAT tsmall = __riscv_vfadd (x, xtail, vlen);
      VFLOAT y_tiny = __riscv_vfadd (tsmall, fp_posOne, vlen);
      y_tiny = U_AS_F (VXOR_VV (F_AS_U (y_tiny), sign_mask, vlen));
      y = VMERGE_VVM_F (y, y_tiny, tiny, vlen);

      /* Huge: underflow to 0 for x<0, overflow to +Inf for x>0.  */
      VBOOL x_is_neg = __riscv_vmflt (x, fp_posZero, vlen);
      VFLOAT y_huge = VFMV_VF (fp_posInf, vlen);
      y_huge = VFMERGE_VFM (y_huge, fp_posZero, x_is_neg, vlen);
      y_huge = U_AS_F (VXOR_VV (F_AS_U (y_huge), sign_mask, vlen));
      y = VMERGE_VVM_F (y, y_huge, huge, vlen);
    }

  return y;
}

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT vx, VFLOAT vy, size_t vlen)
{
  VFLOAT vz, vz_special;
  const struct v_pow_data *data = ptr_barrier (&__v_pow_data);

  SET_ROUNDTONEAREST;

  VUINT xclass = __riscv_vfclass (vx, vlen);
  VUINT yclass = __riscv_vfclass (vy, vlen);

  VBOOL x_is_zero, y_is_zero, x_is_inf, y_is_inf, x_is_nan, y_is_nan;
  VBOOL x_is_snan, y_is_snan;
  IDENTIFY (xclass, class_Zero, x_is_zero, vlen);
  IDENTIFY (yclass, class_Zero, y_is_zero, vlen);
  IDENTIFY (xclass, class_Inf, x_is_inf, vlen);
  IDENTIFY (yclass, class_Inf, y_is_inf, vlen);
  IDENTIFY (xclass, class_NaN, x_is_nan, vlen);
  IDENTIFY (yclass, class_NaN, y_is_nan, vlen);
  IDENTIFY (xclass, class_sNaN, x_is_snan, vlen);
  IDENTIFY (yclass, class_sNaN, y_is_snan, vlen);

  /* Extract sign bit of x for handling negative bases, including -0/-Inf.  */
  const uint64_t sign_bit = 0x8000000000000000ULL;
  VUINT x_bits = F_AS_U (vx);
  VBOOL x_sign =
      __riscv_vmsgtu (__riscv_vand (x_bits, sign_bit, vlen), 0, vlen);

  /* Determine if y is an integer, and if it is odd, using the same bit
     decomposition as scalar checkint().  In particular, do not convert
     non-integral lanes to integers: vfcvt.x.f would raise NX for those lanes
     even though scalar pow's classification path is exception-free.  */
  VFLOAT abs_y = __riscv_vfabs (vy, vlen);
  VUINT abs_y_bits = F_AS_U (abs_y);
  VUINT y_exp = __riscv_vand (__riscv_vsrl (abs_y_bits, 52, vlen), 0x7ff,
                              vlen);
  VUINT y_frac = __riscv_vand (abs_y_bits, 0x000fffffffffffffULL, vlen);
  VBOOL y_exp_ge_bias = __riscv_vmsgeu (y_exp, 0x3ff, vlen);
  VBOOL y_exp_lt_2p52 = __riscv_vmsltu (y_exp, 0x433, vlen);
  VBOOL y_exp_eq_2p52 = __riscv_vmseq (y_exp, 0x433, vlen);
  VBOOL y_exp_ge_2p53 = __riscv_vmsgeu (y_exp, 0x434, vlen);

  VUINT shift = __riscv_vrsub (y_exp, 0x433, vlen);
  VBOOL shift_zero = __riscv_vmor (__riscv_vmsltu (y_exp, 0x3ff, vlen),
                                   __riscv_vmnot (y_exp_lt_2p52, vlen),
                                   vlen);
  shift = __riscv_vmerge (shift, 0, shift_zero, vlen);
  VUINT low_mask = __riscv_vsub (
      __riscv_vsll (VMVU_VX (1, vlen), shift, vlen), 1, vlen);
  VBOOL low_fraction_zero
      = __riscv_vmseq (__riscv_vand (y_frac, low_mask, vlen), 0, vlen);
  VBOOL y_is_int_bounded = __riscv_vmand (
      __riscv_vmand (y_exp_ge_bias, y_exp_lt_2p52, vlen),
      low_fraction_zero, vlen);
  VBOOL y_is_int_large = __riscv_vmand (
      __riscv_vmor (y_exp_eq_2p52, y_exp_ge_2p53, vlen),
      __riscv_vmnot (y_is_nan, vlen), vlen);
  VBOOL y_is_int = __riscv_vmor (y_is_zero,
                                 __riscv_vmor (y_is_int_bounded,
                                               y_is_int_large, vlen), vlen);

  VUINT significand = __riscv_vor (y_frac, 0x0010000000000000ULL, vlen);
  VUINT odd_bit = __riscv_vsll (VMVU_VX (1, vlen), shift, vlen);
  VBOOL y_int_is_odd = __riscv_vmsgtu (
      __riscv_vand (significand, odd_bit, vlen), 0, vlen);
  VBOOL y_odd_range = __riscv_vmand (
      y_exp_ge_bias, __riscv_vmsltu (y_exp, 0x434, vlen), vlen);
  VBOOL y_is_odd = __riscv_vmand (
      __riscv_vmand (y_is_int, y_odd_range, vlen), y_int_is_odd, vlen);

  /* sign_bias sets the sign of the final result in the exp core.  */
  VBOOL neg_and_odd = __riscv_vmand (x_sign, y_is_odd, vlen);
  VUINT sign_bias =
      __riscv_vmerge (VMVU_VX (0, vlen), RVV_POW_EXP_SIGN_BIAS, neg_and_odd,
                      vlen);

  VBOOL x_eq_one = __riscv_vmseq (F_AS_U (vx), 0x3ff0000000000000ULL, vlen);
  VBOOL x_is_special =
      __riscv_vmor (__riscv_vmor (x_is_zero, x_is_inf, vlen), x_is_nan, vlen);
  VBOOL y_is_special =
      __riscv_vmor (__riscv_vmor (y_is_zero, y_is_inf, vlen), y_is_nan, vlen);
  VBOOL x_is_finite = __riscv_vmnot (x_is_special, vlen);
  VBOOL x_neg_finite_nonzero
      = __riscv_vmand (x_sign,
                       __riscv_vmand (x_is_finite,
                                      __riscv_vmnot (x_is_zero, vlen), vlen),
                       vlen);
  VBOOL x_neg_nonint
      = __riscv_vmand (x_neg_finite_nonzero, __riscv_vmnot (y_is_int, vlen),
                       vlen);
  VBOOL domain_err = __riscv_vmand (x_neg_nonint,
                                    __riscv_vmnot (y_is_special, vlen), vlen);
  VBOOL core_special = __riscv_vmor (x_is_special, y_is_special, vlen);
  core_special = __riscv_vmor (core_special, domain_err, vlen);
  core_special = __riscv_vmor (core_special, x_eq_one, vlen);

  /* Prepare every special result before sanitizing the finite core.  */
  vz_special = VFMV_VF (fp_posZero, vlen);

  if (__riscv_vcpop (domain_err, vlen) > 0)
    {
      VFLOAT invalid_input
          = VFMERGE_VFM (vx, fp_sNaN, domain_err, vlen);
      VFLOAT invalid
          = VFADD_VV_M (domain_err, invalid_input, invalid_input, vlen);
      vz_special = VMERGE_VVM_F (vz_special, invalid, domain_err, vlen);
    }

  VBOOL xy_nan = __riscv_vmor (x_is_nan, y_is_nan, vlen);
  if (__riscv_vcpop (xy_nan, vlen) > 0)
    {
      VFLOAT xy_sum = VFADD_VV_M (xy_nan, vx, vy, vlen);
      vz_special = VMERGE_VVM_F (vz_special, xy_sum, xy_nan, vlen);
    }

  /* y is +/-Inf: result depends on |x| relative to 1.  */
  if (__riscv_vcpop (y_is_inf, vlen) > 0)
    {
      VFLOAT abs_x = __riscv_vfabs (vx, vlen);
      VFLOAT abs_x_cmp = VFMERGE_VFM (abs_x, fp_posZero, x_is_nan, vlen);
      VBOOL absx_eq1 = __riscv_vmfeq (abs_x_cmp, fp_posOne, vlen);
      VBOOL absx_lt1 = __riscv_vmflt (abs_x_cmp, fp_posOne, vlen);
      VBOOL absx_gt1 = __riscv_vmfgt (abs_x_cmp, fp_posOne, vlen);
      VBOOL y_neg_inf;
      IDENTIFY (yclass, class_negInf, y_neg_inf, vlen);
      VBOOL no_x_nan = __riscv_vmnot (x_is_nan, vlen);
      VBOOL y_pos_inf = __riscv_vmand (y_is_inf,
                                       __riscv_vmnot (y_neg_inf, vlen), vlen);

      VBOOL want_zero =
          __riscv_vmor (__riscv_vmand (absx_lt1, y_pos_inf, vlen),
                        __riscv_vmand (absx_gt1, y_neg_inf, vlen), vlen);
      VBOOL y_inf_valid = __riscv_vmand (y_is_inf, no_x_nan, vlen);
      VBOOL y_inf_to_inf = y_inf_valid;
      y_inf_to_inf = __riscv_vmandn (y_inf_to_inf, want_zero, vlen);
      y_inf_to_inf = __riscv_vmandn (y_inf_to_inf, absx_eq1, vlen);
      VBOOL y_inf_to_zero = __riscv_vmand (y_inf_valid, want_zero, vlen);
      vz_special = VFMERGE_VFM (vz_special, fp_posInf, y_inf_to_inf, vlen);
      vz_special = VFMERGE_VFM (vz_special, fp_posZero, y_inf_to_zero, vlen);
      VBOOL y_inf_to_one = __riscv_vmand (y_inf_valid, absx_eq1, vlen);
      vz_special = VFMERGE_VFM (vz_special, fp_posOne, y_inf_to_one, vlen);
    }

  /* x is +/-0; keep the divide-by-zero operation restricted to its lanes.  */
  if (__riscv_vcpop (x_is_zero, vlen) > 0)
    {
      VFLOAT vy_cmp = VFMERGE_VFM (vy, fp_posZero, y_is_nan, vlen);
      VBOOL y_lt0 = __riscv_vmflt (vy_cmp, fp_posZero, vlen);
      VBOOL y_gt0 = __riscv_vmfgt (vy_cmp, fp_posZero, vlen);
      VBOOL zero_sign = __riscv_vmand (x_sign, y_is_odd, vlen);
      VUINT zero_sign_mask =
          __riscv_vmerge (VMVU_VX (0, vlen), sign_bit, zero_sign, vlen);
      VFLOAT pos_zero = VFMV_VF (fp_posZero, vlen);
      VFLOAT zero_signed
          = U_AS_F (VXOR_VV (F_AS_U (pos_zero), zero_sign_mask, vlen));
      VBOOL not_y_inf = __riscv_vmnot (y_is_inf, vlen);
      VBOOL mask_inf = __riscv_vmand (x_is_zero, y_lt0, vlen);
      mask_inf = __riscv_vmand (mask_inf, not_y_inf, vlen);
      VBOOL mask_zero = __riscv_vmand (x_is_zero, y_gt0, vlen);
      mask_zero = __riscv_vmand (mask_zero, not_y_inf, vlen);
      if (__riscv_vcpop (mask_inf, vlen) > 0)
        {
          VFLOAT zero_for_div
              = VFMERGE_VFM (zero_signed, fp_posOne,
                             __riscv_vmnot (mask_inf, vlen), vlen);
          VFLOAT divzero
              = __riscv_vfdiv (VFMV_VF (fp_posOne, vlen), zero_for_div, vlen);
          vz_special = VMERGE_VVM_F (vz_special, divzero, mask_inf, vlen);
        }
      vz_special = VMERGE_VVM_F (vz_special, zero_signed, mask_zero, vlen);
    }

  /* x is +/-Inf.  */
  if (__riscv_vcpop (x_is_inf, vlen) > 0)
    {
      VFLOAT vy_cmp = VFMERGE_VFM (vy, fp_posZero, y_is_nan, vlen);
      VBOOL y_lt0 = __riscv_vmflt (vy_cmp, fp_posZero, vlen);
      VBOOL y_gt0 = __riscv_vmfgt (vy_cmp, fp_posZero, vlen);
      VBOOL inf_sign = __riscv_vmand (x_sign, y_is_odd, vlen);
      VUINT inf_sign_mask =
          __riscv_vmerge (VMVU_VX (0, vlen), sign_bit, inf_sign, vlen);
      VFLOAT pos_inf = VFMV_VF (fp_posInf, vlen);
      VFLOAT pos_zero = VFMV_VF (fp_posZero, vlen);
      VFLOAT inf_signed
          = U_AS_F (VXOR_VV (F_AS_U (pos_inf), inf_sign_mask, vlen));
      VFLOAT zero_signed
          = U_AS_F (VXOR_VV (F_AS_U (pos_zero), inf_sign_mask, vlen));
      VBOOL mask_inf = __riscv_vmand (x_is_inf, y_gt0, vlen);
      VBOOL mask_zero = __riscv_vmand (x_is_inf, y_lt0, vlen);
      mask_inf = __riscv_vmand (mask_inf, __riscv_vmnot (y_is_nan, vlen), vlen);
      mask_zero = __riscv_vmand (mask_zero, __riscv_vmnot (y_is_nan, vlen), vlen);
      vz_special = VMERGE_VVM_F (vz_special, inf_signed, mask_inf, vlen);
      vz_special = VMERGE_VVM_F (vz_special, zero_signed, mask_zero, vlen);
    }

  /* x == 1 and y == 0 are identity overrides unless the other operand is
     signaling NaN; scalar pow must then preserve the NaN and invalid event.  */
  VBOOL x_one_identity = __riscv_vmandn (x_eq_one, y_is_snan, vlen);
  VBOOL y_zero_identity = __riscv_vmandn (y_is_zero, x_is_snan, vlen);
  vz_special = VFMERGE_VFM (vz_special, fp_posOne, x_one_identity, vlen);
  vz_special = VFMERGE_VFM (vz_special, fp_posOne, y_zero_identity, vlen);

  /* Main computation uses log/exp only on sanitized finite lanes.  */
  VFLOAT ax = __riscv_vfabs (vx, vlen);
  ax = VFMERGE_VFM (ax, fp_posOne, core_special, vlen);
  VFLOAT vy_core = VFMERGE_VFM (vy, fp_posZero, core_special, vlen);

  /* Compute log(ax) = hi + tail.  */
  VFLOAT log_tail;
  VFLOAT log_hi = RVV_LIBMVEC_IMPL_NAME (pow_log_inline, LMUL)
    (data, ax, &log_tail, vlen);

  /* ehi = y * log_hi; elo = y * log_tail + fma (y, log_hi, -ehi).  */
  VFLOAT ehi = __riscv_vfmul (vy_core, log_hi, vlen);
  VFLOAT emid = __riscv_vfmsub (vy_core, log_hi, ehi, vlen);
  VFLOAT elo = VFMACC_VV (emid, vy_core, log_tail, vlen);

  vz = RVV_LIBMVEC_IMPL_NAME (pow_exp_inline, LMUL)
    (data, ehi, elo, sign_bias, vlen);
  /* Merge the prepared exceptional results only after the finite core.  */
  vz = VMERGE_VVM_F (vz, vz_special, core_special, vlen);

  RESTORE_FRM;
  return vz;
}
