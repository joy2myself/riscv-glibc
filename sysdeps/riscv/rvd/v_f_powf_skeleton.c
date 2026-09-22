/* Single-precision vector powf function.
   powf(x, y) computes x raised to the power y.

   Based on glibc's scalar powf (sysdeps/ieee754/flt-32/e_powf.c) and the
   associated log2/exp2 tables (sysdeps/ieee754/flt-32/e_powf_log2_data.c,
   sysdeps/ieee754/flt-32/e_exp2f_data.c), adapted to RVV.

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
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT vx, VFLOAT vy, size_t vlen)
{
  const struct v_powf_data *data = ptr_barrier (&__v_powf_data);
  VFLOAT_F vz;
  VUINT_F ix, iy, ax;
  VUINT_F vclass_x, vclass_y;
  VBOOL_F x_zero, x_inf, x_nan, x_snan, x_special;
  VBOOL_F y_zero, y_inf, y_nan, y_snan, y_special;
  VBOOL_F x_neg, y_neg;
  VBOOL_F y_conv_ok, y_large;
  VINT_F y_i;
  VFLOAT_F y_i_f;
  VBOOL_F y_is_int_small, y_is_int, y_is_odd;
  VBOOL_F sign_flip;
  VBOOL_F domain_err;

  SET_ROUNDTONEAREST;

  /* Classify inputs up front for special-case handling and to avoid spurious
     exceptions in the main path (align with v_f_logf_skeleton style).  */
  vclass_x = __riscv_vfclass (vx, vlen);
  vclass_y = __riscv_vfclass (vy, vlen);

  IDENTIFY_F (vclass_x, class_sNaN, x_snan, vlen);
  IDENTIFY_F (vclass_x, class_NaN, x_nan, vlen);
  IDENTIFY_F (vclass_x, class_Inf, x_inf, vlen);
  IDENTIFY_F (vclass_x, class_Zero, x_zero, vlen);

  IDENTIFY_F (vclass_y, class_sNaN, y_snan, vlen);
  IDENTIFY_F (vclass_y, class_NaN, y_nan, vlen);
  IDENTIFY_F (vclass_y, class_Inf, y_inf, vlen);
  IDENTIFY_F (vclass_y, class_Zero, y_zero, vlen);

  x_special = __riscv_vmor (x_zero, __riscv_vmor (x_inf, x_nan, vlen), vlen);
  y_special = __riscv_vmor (y_zero, __riscv_vmor (y_inf, y_nan, vlen), vlen);

  ix = F_AS_U_F (vx);
  iy = F_AS_U_F (vy);
  ax = __riscv_vand (ix, RVV_POWF_ABS_MASK, vlen);

  x_neg = __riscv_vmsne (__riscv_vand (ix, RVV_POWF_SIGN_MASK, vlen), 0,
                         vlen);
  y_neg = __riscv_vmsne (__riscv_vand (iy, RVV_POWF_SIGN_MASK, vlen), 0,
                         vlen);

  /* Integer classification for y (only needed for finite y; avoid conversion
     on huge or non-finite lanes).  Sanitize exceptional lanes before the
     floating-point comparison: comparing a signaling/quiet NaN would
     otherwise raise an exception even though that lane is handled below.  */
  VFLOAT_F y_class = VFMERGE_VFM (vy, fp_posZero_f, y_special, vlen);
  y_large
      = __riscv_vmfge (__riscv_vfabs (y_class, vlen), RVV_POWF_INT_BOUND_F,
                       vlen);
  y_conv_ok = __riscv_vmandn (__riscv_vmnot (y_special, vlen), y_large, vlen);

  VFLOAT_F y_safe
      = VFMERGE_VFM (vy, fp_posZero_f, __riscv_vmnot (y_conv_ok, vlen), vlen);
  y_i = __riscv_vfcvt_x (y_safe, vlen);
  y_i_f = __riscv_vfcvt_f (y_i, vlen);
  y_is_int_small = __riscv_vmfeq (y_i_f, y_safe, vlen);
  y_is_int = __riscv_vmor (y_large, y_is_int_small, vlen);

  VINT_F y_lsb = __riscv_vand (y_i, 1, vlen);
  VBOOL_F y_lsb_one = __riscv_vmseq (y_lsb, 1, vlen);
  y_is_odd = __riscv_vmand (y_is_int_small, y_lsb_one, vlen);

  sign_flip = __riscv_vmand (x_neg, y_is_odd, vlen);

  /* Domain error: x < 0 with finite non-integer y.  (x < 0 with y = +/-Inf is
     handled by the y-special path per IEEE rules.)  */
  domain_err
      = __riscv_vmand (x_neg,
                       __riscv_vmandn (__riscv_vmnot (y_inf, vlen), y_is_int,
                                       vlen),
                       vlen);
  domain_err = __riscv_vmand (domain_err, __riscv_vmnot (x_special, vlen),
                              vlen);
  domain_err = __riscv_vmand (domain_err, __riscv_vmnot (y_special, vlen), vlen);

  /* Main path: powf(x, y) = exp2(y * log2(|x|)), using glibc tables and
     double-precision intermediate for accuracy.  Mask special lanes out to
     avoid spurious exceptions.  */
  VBOOL_F special_args
      = __riscv_vmor (x_special, __riscv_vmor (y_special, domain_err, vlen),
                      vlen);
  VFLOAT_F x_abs = U_AS_F_F (ax);
  VFLOAT_F x_main = VFMERGE_VFM (x_abs, fp_posOne_f, special_args, vlen);
#if LMUL != 8
  VFLOAT_F y_main = VFMERGE_VFM (vy, fp_posZero_f, special_args, vlen);
#endif

  /* Normalize subnormal x for the log2 core.  */
  VUINT_F ix_main = F_AS_U_F (x_main);
  VBOOL_F x_sub = __riscv_vmsltu (ix_main, RVV_POWF_MIN_NORM, vlen);
  /* Only subnormal lanes need the scaling.  Keep all other lanes benign so a
     large finite x cannot raise overflow/inexact in a result that is later
     discarded.  */
  VFLOAT_F x_scale_input
      = VFMERGE_VFM (x_main, fp_posOne_f, __riscv_vmnot (x_sub, vlen), vlen);
  VFLOAT_F x_scaled
      = __riscv_vfmul (x_scale_input, RVV_POWF_SUBNORM_SCALE_F, vlen);
  VUINT_F ix_scaled = __riscv_vand (F_AS_U_F (x_scaled), RVV_POWF_ABS_MASK,
                                    vlen);
  VUINT_F ix_norm = __riscv_vsub (ix_scaled,
                                  RVV_POWF_SUBNORM_EXPONENT_BIAS, vlen);
  ix_main = __riscv_vmerge (ix_main, ix_norm, x_sub, vlen);

#if LMUL == 8
  /* LMUL=8 cannot widen f32 -> f64 with a matching element count.  Split the
     register group into two LMUL=4 subvectors and reuse the same
     double-intermediate glibc powf port used by the narrower variants.  */
  size_t vlmax_m4 = __riscv_vsetvlmax_e32m4 ();
  size_t lo_vlen = vlen < vlmax_m4 ? vlen : vlmax_m4;
  size_t hi_vlen = vlen > vlmax_m4 ? vlen - vlmax_m4 : 0;

  vfloat32m4_t vx_lo = __riscv_vget_v_f32m8_f32m4 (vx, 0);
  vfloat32m4_t vy_lo = __riscv_vget_v_f32m8_f32m4 (vy, 0);
  vfloat32m4_t vz_lo
      = RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, 4) (vx_lo, vy_lo, lo_vlen);

  vfloat32m4_t vz_hi = __riscv_vfmv_v_f_f32m4 (0.0f, vlmax_m4);
  if (hi_vlen > 0)
    {
      vfloat32m4_t vx_hi = __riscv_vget_v_f32m8_f32m4 (vx, 1);
      vfloat32m4_t vy_hi = __riscv_vget_v_f32m8_f32m4 (vy, 1);
      vz_hi = RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, 4) (vx_hi, vy_hi, hi_vlen);
    }

  VFLOAT_F vz_m8 = __riscv_vcreate_v_f32m4_f32m8 (vz_lo, vz_hi);
  RESTORE_FRM;
  return vz_m8;

#else
  /* Double-precision path for LMUL <= 4.  */
  /* log2_inline(ix_main) in double.  */
  RVV_POWF_VU32_IDX tmp
      = __riscv_vsub (ix_main, RVV_POWF_LOG2_OFF, vlen);
  RVV_POWF_VU32_IDX idx
      = __riscv_vsrl (tmp, RVV_POWF_FLOAT_MANTISSA_BITS
                            - RVV_POWF_LOG2_TABLE_BITS,
                      vlen);
  idx = __riscv_vand (idx, RVV_POWF_LOG2_TABLE_MASK, vlen);
  RVV_POWF_VU32_IDX top
      = __riscv_vand (tmp, RVV_POWF_LOG2_TOP_MASK, vlen);
  RVV_POWF_VU32_IDX iz = __riscv_vsub (ix_main, top, vlen);

  VINT_F k_i = __riscv_vsra (U_AS_I_F (top), 23, vlen);
  VFLOAT_F k_f = __riscv_vfcvt_f (k_i, vlen);
  RVV_POWF_VDOUBLE k_d = RVV_POWF_VFWCVT (k_f, vlen);

  VFLOAT_F z_f = U_AS_F_F (iz);
  RVV_POWF_VDOUBLE z_d = RVV_POWF_VFWCVT (z_f, vlen);

  RVV_POWF_VU32_IDX idx_bytes = __riscv_vsll (idx, 3, vlen);
  RVV_POWF_VDOUBLE invc
      = RVV_POWF_VGATHER_F64 (data->log2_invc, idx_bytes, vlen);
  RVV_POWF_VDOUBLE logc
      = RVV_POWF_VGATHER_F64 (data->log2_logc, idx_bytes, vlen);

  RVV_POWF_VDOUBLE r
      = __riscv_vfsub (__riscv_vfmul (z_d, invc, vlen), 1.0, vlen);
  RVV_POWF_VDOUBLE y0 = __riscv_vfadd (logc, k_d, vlen);

  RVV_POWF_VDOUBLE r2 = __riscv_vfmul (r, r, vlen);
  RVV_POWF_VDOUBLE r4 = __riscv_vfmul (r2, r2, vlen);

  RVV_POWF_VDOUBLE ypoly
      = __riscv_vfmadd (RVV_POWF_VDFMV_VF (RVV_POWF_LOG2_POLY0_D, vlen), r,
                        RVV_POWF_VDFMV_VF (RVV_POWF_LOG2_POLY1_D, vlen),
                        vlen);
  RVV_POWF_VDOUBLE ppoly
      = __riscv_vfmadd (RVV_POWF_VDFMV_VF (RVV_POWF_LOG2_POLY2_D, vlen), r,
                        RVV_POWF_VDFMV_VF (RVV_POWF_LOG2_POLY3_D, vlen),
                        vlen);
  RVV_POWF_VDOUBLE qpoly = __riscv_vfmadd (
      RVV_POWF_VDFMV_VF (RVV_POWF_LOG2_POLY4_D, vlen), r, y0, vlen);
  qpoly = __riscv_vfmadd (ppoly, r2, qpoly, vlen);
  RVV_POWF_VDOUBLE log2x = __riscv_vfmadd (ypoly, r4, qpoly, vlen);

  RVV_POWF_VDOUBLE y_d = RVV_POWF_VFWCVT (y_main, vlen);
  RVV_POWF_VDOUBLE ylogx = __riscv_vfmul (y_d, log2x, vlen);

  /* Clamp huge |ylogx| before exp2 bit-hack.  Without this, exp2_inline's
     shift/reinterpret sequence can generate incorrect NaN/Inf for inputs far
     outside float range (e.g. ylogx ~ -1e26 should underflow to 0).  */
  VBOOL_F ylogx_under = __riscv_vmfle (ylogx, RVV_POWF_YLOGX_UNDER_D, vlen);
  VBOOL_F ylogx_over = __riscv_vmfgt (ylogx, RVV_POWF_YLOGX_OVER_D, vlen);
  VBOOL_F ylogx_clamp = __riscv_vmor (ylogx_under, ylogx_over, vlen);
  RVV_POWF_VDOUBLE ylogx_safe = __riscv_vmerge (
      ylogx, RVV_POWF_VDFMV_VF (0.0, vlen), ylogx_clamp, vlen);

  if (__riscv_vcpop (ylogx_under, vlen) > 0)
    {
      volatile VFLOAT_F uf_probe
          = __riscv_vfmul (VFMERGE_VFM (VFMV_VF_F (0.0f, vlen),
                                        RVV_POWF_MIN_SUBNORMAL_F,
                                        ylogx_under, vlen),
                           VFMERGE_VFM (VFMV_VF_F (0.0f, vlen),
                                        RVV_POWF_MIN_SUBNORMAL_F,
                                        ylogx_under, vlen),
                           vlen);
      (void) uf_probe;
    }

  if (__riscv_vcpop (ylogx_over, vlen) > 0)
    {
      volatile VFLOAT_F of_probe
          = __riscv_vfmul (VFMERGE_VFM (VFMV_VF_F (1.0f, vlen),
                                        RVV_POWF_MAX_FINITE_F,
                                        ylogx_over, vlen),
                           2.0f, vlen);
      (void) of_probe;
    }

  /* exp2_inline(ylogx) in double.  */
  RVV_POWF_VDOUBLE kd
      = __riscv_vfadd (ylogx_safe, RVV_POWF_EXP2_SHIFT_D, vlen);
  RVV_POWF_VU64 ki = RVV_POWF_F64_AS_U64 (kd);
  kd = __riscv_vfsub (kd, RVV_POWF_EXP2_SHIFT_D, vlen);
  RVV_POWF_VDOUBLE rr = __riscv_vfsub (ylogx_safe, kd, vlen);

  RVV_POWF_VU64 tab_i
      = __riscv_vand (ki, RVV_POWF_EXP2_TABLE_MASK, vlen);
  tab_i = __riscv_vsll (tab_i, 3, vlen);
  RVV_POWF_VU64 t
      = RVV_POWF_VGATHER_U64 (data->exp2_tab, tab_i, vlen);
  RVV_POWF_VU64 ki_shift
      = __riscv_vsll (ki, RVV_POWF_DOUBLE_MANTISSA_BITS
                            - RVV_POWF_EXP2_TABLE_BITS,
                      vlen);
  t = __riscv_vadd (t, ki_shift, vlen);
  RVV_POWF_VDOUBLE s = RVV_POWF_U64_AS_F64 (t);

  RVV_POWF_VDOUBLE z
      = __riscv_vfmadd (RVV_POWF_VDFMV_VF (RVV_POWF_EXP2_C0_D, vlen), rr,
                        RVV_POWF_VDFMV_VF (RVV_POWF_EXP2_C1_D, vlen), vlen);
  RVV_POWF_VDOUBLE rr2 = __riscv_vfmul (rr, rr, vlen);
  RVV_POWF_VDOUBLE p
      = __riscv_vfmadd (RVV_POWF_VDFMV_VF (RVV_POWF_EXP2_C2_D, vlen), rr,
                        RVV_POWF_VDFMV_VF (1.0, vlen), vlen);
  p = __riscv_vfmadd (z, rr2, p, vlen);
  RVV_POWF_VDOUBLE pow_d = __riscv_vfmul (p, s, vlen);

  vz = RVV_POWF_VFNCVT (pow_d, vlen);
  vz = VMERGE_VVM_F (vz, VFMV_VF_F (fp_posZero_f, vlen), ylogx_under, vlen);
  vz = VMERGE_VVM_F (vz, VFMV_VF_F (fp_posInf_f, vlen), ylogx_over, vlen);
#endif /* LMUL == 8 */

  /* Apply sign for negative base and odd integer exponent.  */
  VINT_F sign_mask_i = VMVI_VX_F (0, vlen);
  sign_mask_i
      = VMERGE_VXM_I (sign_mask_i, (int32_t) RVV_POWF_SIGN_MASK,
                       sign_flip, vlen);
  vz = U_AS_F_F (VXOR_VV (F_AS_U_F (vz), I_AS_U_F (sign_mask_i), vlen));

  /* Special-case: negative base with non-integer exponent -> NaN.  */
  if (__riscv_vcpop (domain_err, vlen) > 0)
    {
      VFLOAT_F nan = VFMV_VF_F (fp_sNaN_f, vlen);
      nan = __riscv_vfadd (domain_err, nan, nan, vlen);
      vz = VMERGE_VVM_F (vz, nan, domain_err, vlen);
    }

  /* y in {0, +/-Inf, NaN}: handle per scalar powf rules.  */
  if (__riscv_vcpop (y_special, vlen) > 0)
    {
      VFLOAT_F vy_special = VFMV_VF_F (fp_posOne_f, vlen);

      /* y == 0: return 1, except propagate sNaN x via x+y.  */
      VBOOL_F mask = __riscv_vmand (y_zero, x_snan, vlen);
      if (__riscv_vcpop (mask, vlen) > 0)
        {
          VFLOAT_F xy = __riscv_vfadd (mask, vx, vy, vlen);
          vy_special = VMERGE_VVM_F (vy_special, xy, mask, vlen);
        }

      /* x == 1: return 1, except propagate sNaN y via x+y.  */
      VBOOL_F x_eq_one = __riscv_vmseq (ix, RVV_POWF_ONE_BITS, vlen);
      mask
          = __riscv_vmand (y_special, __riscv_vmand (x_eq_one, y_snan, vlen),
                           vlen);
      if (__riscv_vcpop (mask, vlen) > 0)
        {
          VFLOAT_F xy = __riscv_vfadd (mask, vx, vy, vlen);
          vy_special = VMERGE_VVM_F (vy_special, xy, mask, vlen);
        }

      /* Propagate NaNs (excluding y==0 and x==1 handled above).  */
      VBOOL_F nan_any = __riscv_vmor (x_nan, y_nan, vlen);
      mask = __riscv_vmand (y_special, nan_any, vlen);
      mask = __riscv_vmandn (mask, __riscv_vmor (y_zero, x_eq_one, vlen), vlen);
      if (__riscv_vcpop (mask, vlen) > 0)
        {
          VFLOAT_F xy = __riscv_vfadd (mask, vx, vy, vlen);
          vy_special = VMERGE_VVM_F (vy_special, xy, mask, vlen);
        }

      /* y is +/-Inf and |x| != 1: result is 0 or +Inf depending on |x| < 1.  */
      VBOOL_F abs_eq_one = __riscv_vmseq (ax, RVV_POWF_ONE_BITS, vlen);
      VBOOL_F abs_lt_one = __riscv_vmsltu (ax, RVV_POWF_ONE_BITS, vlen);
      VBOOL_F abs_ne_one = __riscv_vmnot (abs_eq_one, vlen);
      VBOOL_F y_inf_ne_one = __riscv_vmand (
          __riscv_vmand (y_inf, abs_ne_one, vlen),
          __riscv_vmnot (nan_any, vlen), vlen);

      /* cond_zero = abs_lt_one XOR y_neg.  */
      VBOOL_F abs_lt_one_n = __riscv_vmnot (abs_lt_one, vlen);
      VBOOL_F y_neg_n = __riscv_vmnot (y_neg, vlen);
      VBOOL_F cond_zero
          = __riscv_vmor (__riscv_vmand (abs_lt_one, y_neg_n, vlen),
                          __riscv_vmand (abs_lt_one_n, y_neg, vlen), vlen);
      VBOOL_F to_zero = __riscv_vmand (y_inf_ne_one, cond_zero, vlen);
      VBOOL_F to_inf = __riscv_vmand (y_inf_ne_one, __riscv_vmnot (cond_zero, vlen),
                                      vlen);
      vy_special = VMERGE_VVM_F (vy_special, VFMV_VF_F (fp_posZero_f, vlen),
                                 to_zero, vlen);
      vy_special = VMERGE_VVM_F (vy_special, VFMV_VF_F (fp_posInf_f, vlen),
                                 to_inf, vlen);

      vz = VMERGE_VVM_F (vz, vy_special, y_special, vlen);
    }

  /* x in {0, +/-Inf, NaN} with finite nonzero y: compute via x2 and reciprocal,
     with sign for odd integer y.  */
  if (__riscv_vcpop (x_special, vlen) > 0)
    {
      VBOOL_F x_special_finite_y = __riscv_vmandn (x_special, y_special, vlen);
      if (__riscv_vcpop (x_special_finite_y, vlen) > 0)
        {
          VFLOAT_F x_fix = VFMERGE_VFM (
              vx, fp_posOne_f, __riscv_vmnot (x_special_finite_y, vlen), vlen);
          VFLOAT_F x2 = __riscv_vfmul (x_fix, x_fix, vlen);

          VBOOL_F x_sign_flip = __riscv_vmand (x_neg, y_is_odd, vlen);
          VINT_F sgn = VMVI_VX_F (0, vlen);
          sgn = VMERGE_VXM_I (sgn, (int32_t) RVV_POWF_SIGN_MASK,
                              x_sign_flip, vlen);
          x2 = U_AS_F_F (VXOR_VV (F_AS_U_F (x2), I_AS_U_F (sgn), vlen));

          VBOOL_F use_inv = __riscv_vmand (x_special_finite_y, y_neg, vlen);
          VFLOAT_F x2_inv = VFMERGE_VFM (
              x2, fp_posOne_f, __riscv_vmnot (use_inv, vlen), vlen);
          VFLOAT_F inv_x2 = __riscv_vfdiv (VFMV_VF_F (fp_posOne_f, vlen),
                                           x2_inv, vlen);
          VFLOAT_F x_res = __riscv_vmerge (x2, inv_x2, y_neg, vlen);
          vz = __riscv_vmerge (vz, x_res, x_special_finite_y, vlen);
        }
    }

  RESTORE_FRM;
  return vz;
}

#undef LMUL
#undef SIMDLEN
