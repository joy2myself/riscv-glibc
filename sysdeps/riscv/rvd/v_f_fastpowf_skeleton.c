/* Single-precision vector fastpowf function.
   fastpowf(x, y) computes x^y with reduced precision for speed

   Migrated from SLEEF (src/libm/sleefsimdsp.c:xfastpowf_u3500)

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

#ifndef V_F_FASTPOWF_SKELETON_ONCE
#define V_F_FASTPOWF_SKELETON_ONCE

/* SLEEF constants for exp/log.  */
#define R_LN2f 0x1.715476p+0f
#define L2Uf 0x1.62e400p-1f
#define L2Lf 0x1.7f7d1cp-20f
#define LN2f 0x1.62e430p-1f
#define LOGK3F_INV_075 0x1.555556p+0f
#define LOGK3F_SCALE_SMALL 0x1.0p64f
#define LOGK3F_MIN 0x1.0p-126f
#define POWF_INT_BOUND 16777216.0f
#ifndef FE_TOWARDZERO
# define FE_TOWARDZERO 0x001
#endif

#endif

#undef MAKE_VSTORE
#undef VFSTORE
#define MAKE_VSTORE(A)                                                        \
  __PASTE3 (__PASTE3 (__riscv_vse, BIT_WIDTH, _v_), A,                        \
            __PASTE3 (BIT_WIDTH, m, LMUL))
#define VFSTORE MAKE_VSTORE (f)

#define FASTPOWF_VILOGB2K __PASTE4 (fastpowf_vilogb2k_, LMUL, _, SIMDLEN)
#define FASTPOWF_VLDEXP3 __PASTE4 (fastpowf_vldexp3_, LMUL, _, SIMDLEN)
#define FASTPOWF_VPOW2I __PASTE4 (fastpowf_vpow2i_, LMUL, _, SIMDLEN)
#define FASTPOWF_VLDEXP2 __PASTE4 (fastpowf_vldexp2_, LMUL, _, SIMDLEN)
#define FASTPOWF_VTRUNCATE __PASTE4 (fastpowf_vtruncate_, LMUL, _, SIMDLEN)
#define FASTPOWF_LOGK3F __PASTE4 (fastpowf_logk3f_, LMUL, _, SIMDLEN)
#define FASTPOWF_EXPK3F __PASTE4 (fastpowf_expk3f_, LMUL, _, SIMDLEN)

static inline VINT_F
FASTPOWF_VILOGB2K (VFLOAT_F d, size_t vlen)
{
  VUINT_F u = F_AS_U (d);
  VUINT_F e = __riscv_vsrl (u, MAN_LEN_F, vlen);
  e = __riscv_vand (e, 0xff, vlen);
  return __riscv_vsub (U_AS_I (e), VMVI_VX (EXP_BIAS_F, vlen), vlen);
}

static inline VFLOAT_F
FASTPOWF_VLDEXP3 (VFLOAT_F d, VINT_F q, size_t vlen)
{
  VINT_F q_shift = __riscv_vsll (q, MAN_LEN_F, vlen);
  return I_AS_F (__riscv_vadd (F_AS_I (d), q_shift, vlen));
}

static inline VFLOAT_F
FASTPOWF_VPOW2I (VINT_F q, size_t vlen)
{
  VINT_F q_bias = __riscv_vadd (q, VMVI_VX (EXP_BIAS_F, vlen), vlen);
  q_bias = __riscv_vsll (q_bias, MAN_LEN_F, vlen);
  return I_AS_F (q_bias);
}

static inline VFLOAT_F
FASTPOWF_VLDEXP2 (VFLOAT_F d, VINT_F e, size_t vlen)
{
  VINT_F e1 = __riscv_vsra (e, 1, vlen);
  VINT_F e2 = __riscv_vsub (e, e1, vlen);
  VFLOAT_F s1 = FASTPOWF_VPOW2I (e1, vlen);
  VFLOAT_F s2 = FASTPOWF_VPOW2I (e2, vlen);
  return __riscv_vfmul (__riscv_vfmul (d, s1, vlen), s2, vlen);
}

static inline VFLOAT_F
FASTPOWF_VTRUNCATE (VFLOAT_F x, size_t vlen)
{
  unsigned long frm = read_frm ();
  if (frm != FE_TOWARDZERO)
    write_frm (FE_TOWARDZERO);
  VINT_F xi = __riscv_vfcvt_x (x, vlen);
  VFLOAT_F xf = __riscv_vfcvt_f (xi, vlen);
  if (frm != FE_TOWARDZERO)
    write_frm (frm);
  return xf;
}

static inline VFLOAT_F
FASTPOWF_LOGK3F (VFLOAT_F d, size_t vlen)
{
  VFLOAT_F x, x2, t, m;
  VBOOL_F o = __riscv_vmflt (d, VFMV_VF (LOGK3F_MIN, vlen), vlen);
  VFLOAT_F d_scaled
      = __riscv_vfmul (d, VFMV_VF (LOGK3F_SCALE_SMALL, vlen), vlen);
  d = VMERGE_VVM_F (d, d_scaled, o, vlen);

  VFLOAT_F scaled = __riscv_vfmul (d, VFMV_VF (LOGK3F_INV_075, vlen), vlen);
  VINT_F e = FASTPOWF_VILOGB2K (scaled, vlen);
  m = FASTPOWF_VLDEXP3 (d, __riscv_vneg (e, vlen), vlen);
  e = VMERGE_VVM_I (e, __riscv_vsub (e, VMVI_VX (64, vlen), vlen), o, vlen);

  x = __riscv_vfdiv (__riscv_vfsub (m, VFMV_VF (1.0f, vlen), vlen),
                     __riscv_vfadd (m, VFMV_VF (1.0f, vlen), vlen), vlen);
  x2 = __riscv_vfmul (x, x, vlen);

  t = VFMV_VF (0.2392828464508056640625f, vlen);
  t = VFMADD_VV (t, x2, VFMV_VF (0.28518211841583251953125f, vlen), vlen);
  t = VFMADD_VV (t, x2, VFMV_VF (0.400005877017974853515625f, vlen), vlen);
  t = VFMADD_VV (t, x2, VFMV_VF (0.666666686534881591796875f, vlen), vlen);
  t = VFMADD_VV (t, x2, VFMV_VF (2.0f, vlen), vlen);

  VFLOAT_F e_f = __riscv_vfcvt_f (e, vlen);
  VFLOAT_F ln2e = __riscv_vfmul (e_f, VFMV_VF (LN2f, vlen), vlen);
  x = VFMADD_VV (x, t, ln2e, vlen);

  return x;
}

static inline VFLOAT_F
FASTPOWF_EXPK3F (VFLOAT_F d, size_t vlen)
{
  VINT_F q
      = __riscv_vfcvt_x (__riscv_vfmul (d, VFMV_VF (R_LN2f, vlen), vlen), vlen);
  VFLOAT_F q_f = __riscv_vfcvt_f (q, vlen);
  VFLOAT_F s = VFMADD_VV (q_f, VFMV_VF (-L2Uf, vlen), d, vlen);
  s = VFMADD_VV (q_f, VFMV_VF (-L2Lf, vlen), s, vlen);

  VFLOAT_F u = VFMV_VF (0.000198527617612853646278381f, vlen);
  u = VFMADD_VV (u, s, VFMV_VF (0.00139304355252534151077271f, vlen), vlen);
  u = VFMADD_VV (u, s, VFMV_VF (0.00833336077630519866943359f, vlen), vlen);
  u = VFMADD_VV (u, s, VFMV_VF (0.0416664853692054748535156f, vlen), vlen);
  u = VFMADD_VV (u, s, VFMV_VF (0.166666671633720397949219f, vlen), vlen);
  u = VFMADD_VV (u, s, VFMV_VF (0.5f, vlen), vlen);

  VFLOAT_F s2 = __riscv_vfmul (s, s, vlen);
  u = VFMADD_VV (s2, u, __riscv_vfadd (s, VFMV_VF (1.0f, vlen), vlen), vlen);
  u = FASTPOWF_VLDEXP2 (u, q, vlen);

  VBOOL_F too_small = __riscv_vmflt (d, VFMV_VF (-104.0f, vlen), vlen);
  u = VMERGE_VVM_F (u, VFMV_VF (0.0f, vlen), too_small, vlen);

  return u;
}

/* SLEEF xfastpowf_u3500 implementation (ULP = 3500).  */

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, VFLOAT y, size_t vlen)
{
  VFLOAT_F x_work, y_work;
  VFLOAT_F result;
  VFLOAT_F vy_special;
  VFLOAT_F abs_y;
  VBOOL_F y_is_int;
  VBOOL_F y_is_odd;
  VBOOL_F x_eq_zero;
  VBOOL_F x_eq_one;
  VBOOL_F x_abs_eq_one;
  VBOOL_F y_eq_zero;
  VBOOL_F y_eq_one;
  VBOOL_F core_skip;
  VBOOL_F special_args;

  SET_ROUNDTONEAREST;

  VUINT_F vclass_x = __riscv_vfclass (x, vlen);
  VUINT_F vclass_y = __riscv_vfclass (y, vlen);
  VBOOL_F x_nan, y_nan, x_inf, y_inf;
  IDENTIFY_F (vclass_x, class_NaN, x_nan, vlen);
  IDENTIFY_F (vclass_y, class_NaN, y_nan, vlen);
  IDENTIFY_F (vclass_x, class_Inf, x_inf, vlen);
  IDENTIFY_F (vclass_y, class_Inf, y_inf, vlen);
  VBOOL_F x_nan_inf = __riscv_vmor (x_nan, x_inf, vlen);
  VBOOL_F y_nan_inf = __riscv_vmor (y_nan, y_inf, vlen);
  VBOOL_F x_finite = __riscv_vmnot (x_nan_inf, vlen);
  VBOOL_F y_finite = __riscv_vmnot (y_nan_inf, vlen);

  /* Use class/bitwise identities instead of signaling FP compares.  */
  IDENTIFY_F (vclass_x, class_Zero, x_eq_zero, vlen);
  IDENTIFY_F (vclass_y, class_Zero, y_eq_zero, vlen);
  VUINT_F x_bits = F_AS_U (x);
  VUINT_F abs_x_bits = __riscv_vand (x_bits, 0x7fffffffu, vlen);
  x_eq_one = __riscv_vmseq (x_bits, 0x3f800000u, vlen);
  x_abs_eq_one = __riscv_vmseq (U_AS_I (abs_x_bits),
                                VMVI_VX (0x3f800000, vlen), vlen);
  y_eq_one = __riscv_vmseq (F_AS_U (y), 0x3f800000u, vlen);

  VFLOAT_F y_cmp = VFMERGE_VFM (y, 0.0f, y_nan_inf, vlen);
  abs_y = __riscv_vfabs (y_cmp, vlen);
  VBOOL_F y_large = __riscv_vmand (y_finite,
                                   __riscv_vmfge (abs_y, POWF_INT_BOUND, vlen),
                                   vlen);
  VBOOL_F y_conv_special
      = __riscv_vmor (__riscv_vmnot (y_finite, vlen), y_large, vlen);
  VFLOAT_F y_int_input = VFMERGE_VFM (y, 0.0f, y_conv_special, vlen);
  VFLOAT_F y_trunc = FASTPOWF_VTRUNCATE (y_int_input, vlen);
  y_is_int = __riscv_vmor (__riscv_vmfeq (y_trunc, y_int_input, vlen),
                           y_large, vlen);
  y_is_int = __riscv_vmand (y_is_int, y_finite, vlen);

  VINT_F y_int = __riscv_vfcvt_x (y_trunc, vlen);
  VBOOL_F y_lt_bound = __riscv_vmflt (abs_y, POWF_INT_BOUND, vlen);
  VINT_F y_int_odd = VAND_VV_I (y_int, VMVI_VX (1, vlen), vlen);
  VBOOL_F y_int_is_odd = __riscv_vmseq (y_int_odd, VMVI_VX (1, vlen), vlen);
  y_is_odd = __riscv_vmand (y_is_int, y_int_is_odd, vlen);
  y_is_odd = __riscv_vmand (y_is_odd, y_lt_bound, vlen);

  VUINT_F sign_bits = __riscv_vand (x_bits, 0x80000000u, vlen);
  VBOOL_F x_neg
      = __riscv_vmseq (U_AS_I (sign_bits), VMVI_VX (0x80000000, vlen), vlen);
  VBOOL_F x_neg_nonzero = __riscv_vmandn (x_neg, x_eq_zero, vlen);
  VBOOL_F x_neg_nonint = __riscv_vmand (x_neg_nonzero,
                                        __riscv_vmnot (y_is_int, vlen), vlen);
  x_neg_nonint = __riscv_vmand (x_neg_nonint, x_finite, vlen);
  x_neg_nonint = __riscv_vmand (x_neg_nonint, y_finite, vlen);

  special_args = __riscv_vmor (x_nan_inf, y_nan_inf, vlen);
  special_args = __riscv_vmor (special_args, x_neg_nonint, vlen);

  core_skip = __riscv_vmor (x_eq_zero, x_eq_one, vlen);
  core_skip = __riscv_vmor (core_skip, y_eq_zero, vlen);
  core_skip = __riscv_vmor (core_skip, y_eq_one, vlen);
  core_skip = __riscv_vmor (core_skip, special_args, vlen);

  /* Prepare special and identity results before entering the finite core.  */
  vy_special = VFMV_VF (0.0f, vlen);
  if (__riscv_vcpop (special_args, vlen) > 0)
    {
      VBOOL_F nan_args = __riscv_vmor (x_nan, y_nan, vlen);
      if (__riscv_vcpop (nan_args, vlen) > 0)
        {
          VFLOAT_F nan_result = VFADD_VV_M (nan_args, x, y, vlen);
          vy_special = VMERGE_VVM_F (vy_special, nan_result, nan_args, vlen);
        }

      if (__riscv_vcpop (x_neg_nonint, vlen) > 0)
        {
          VFLOAT_F invalid_input
              = VFMERGE_VFM (x, fp_sNaN_f, x_neg_nonint, vlen);
          VFLOAT_F invalid
              = VFADD_VV_M (x_neg_nonint, invalid_input, invalid_input, vlen);
          vy_special = VMERGE_VVM_F (vy_special, invalid, x_neg_nonint, vlen);
        }

      VBOOL_F y_pos_inf, y_neg_inf;
      IDENTIFY_F (vclass_y, class_posInf, y_pos_inf, vlen);
      IDENTIFY_F (vclass_y, class_negInf, y_neg_inf, vlen);
      VFLOAT_F ax = __riscv_vfabs (VFMERGE_VFM (x, 0.0f, x_nan, vlen), vlen);
      VBOOL_F ax_gt_1 = __riscv_vmfgt (ax, 1.0f, vlen);
      VBOOL_F ax_lt_1 = __riscv_vmflt (ax, 1.0f, vlen);
      VBOOL_F no_x_nan = __riscv_vmnot (x_nan, vlen);
      VBOOL_F y_inf_to_inf
          = __riscv_vmor (__riscv_vmand (y_pos_inf, ax_gt_1, vlen),
                          __riscv_vmand (y_neg_inf, ax_lt_1, vlen), vlen);
      VBOOL_F y_inf_to_zero
          = __riscv_vmor (__riscv_vmand (y_pos_inf, ax_lt_1, vlen),
                          __riscv_vmand (y_neg_inf, ax_gt_1, vlen), vlen);
      y_inf_to_inf = __riscv_vmand (y_inf_to_inf, no_x_nan, vlen);
      y_inf_to_zero = __riscv_vmand (y_inf_to_zero, no_x_nan, vlen);
      vy_special = VMERGE_VVM_F (vy_special, VFMV_VF (fp_posInf_f, vlen),
                                 y_inf_to_inf, vlen);
      vy_special = VMERGE_VVM_F (vy_special, VFMV_VF (fp_posZero_f, vlen),
                                 y_inf_to_zero, vlen);

      VBOOL_F y_gt_0 = __riscv_vmfgt (y_cmp, 0.0f, vlen);
      VBOOL_F y_lt_0 = __riscv_vmflt (y_cmp, 0.0f, vlen);
      VFLOAT_F xinf_pos = VFMV_VF (fp_posInf_f, vlen);
      VFLOAT_F xinf_zero = VFMV_VF (fp_posZero_f, vlen);
      VFLOAT_F signed_inf = __riscv_vfsgnj (xinf_pos, x, vlen);
      VFLOAT_F signed_zero = __riscv_vfsgnj (xinf_zero, x, vlen);
      signed_inf = VMERGE_VVM_F (xinf_pos, signed_inf, y_is_odd, vlen);
      signed_zero = VMERGE_VVM_F (xinf_zero, signed_zero, y_is_odd, vlen);
      VBOOL_F x_inf_to_inf = __riscv_vmand (x_inf, y_gt_0, vlen);
      VBOOL_F x_inf_to_zero = __riscv_vmand (x_inf, y_lt_0, vlen);
      x_inf_to_inf = __riscv_vmand (x_inf_to_inf, no_x_nan, vlen);
      x_inf_to_inf = __riscv_vmand (x_inf_to_inf,
                                    __riscv_vmnot (y_nan, vlen), vlen);
      x_inf_to_zero = __riscv_vmand (x_inf_to_zero, no_x_nan, vlen);
      x_inf_to_zero = __riscv_vmand (x_inf_to_zero,
                                     __riscv_vmnot (y_nan, vlen), vlen);
      vy_special = VMERGE_VVM_F (vy_special, signed_inf, x_inf_to_inf, vlen);
      vy_special = VMERGE_VVM_F (vy_special, signed_zero, x_inf_to_zero, vlen);
      vy_special = VMERGE_VVM_F (vy_special, VFMV_VF (1.0f, vlen),
                                 __riscv_vmand (y_inf, x_abs_eq_one, vlen), vlen);
    }

  /* Preserve the existing identity precedence in one prepared result.  */
  if (__riscv_vcpop (x_eq_zero, vlen) > 0)
    {
      /* Zero bases need the scalar powf pole/sign rules.  Keep NaN and
         infinite exponents on their already-prepared paths, and trigger
         divide-by-zero only for finite negative exponents.  */
      VBOOL_F zero_finite_y = __riscv_vmand (x_eq_zero, y_finite, vlen);
      VBOOL_F zero_to_inf
          = __riscv_vmand (zero_finite_y, __riscv_vmflt (y_cmp, 0.0f, vlen),
                           vlen);
      VBOOL_F zero_to_zero
          = __riscv_vmand (zero_finite_y, __riscv_vmfgt (y_cmp, 0.0f, vlen),
                           vlen);

      VFLOAT_F signed_zero = __riscv_vfsgnj (VFMV_VF (fp_posZero_f, vlen), x,
                                             vlen);
      signed_zero = VMERGE_VVM_F (VFMV_VF (fp_posZero_f, vlen), signed_zero,
                                  y_is_odd, vlen);
      VFLOAT_F signed_inf = __riscv_vfsgnj (VFMV_VF (fp_posInf_f, vlen), x,
                                            vlen);
      signed_inf = VMERGE_VVM_F (VFMV_VF (fp_posInf_f, vlen), signed_inf,
                                 y_is_odd, vlen);

      if (__riscv_vcpop (zero_to_inf, vlen) > 0)
        {
          VFLOAT_F zero_for_div
              = VFMERGE_VFM (signed_zero, fp_posOne_f,
                             __riscv_vmnot (zero_to_inf, vlen), vlen);
          VFLOAT_F divzero
              = __riscv_vfdiv (VFMV_VF (fp_posOne_f, vlen), zero_for_div,
                               vlen);
          vy_special = VMERGE_VVM_F (vy_special, divzero, zero_to_inf, vlen);
        }
      vy_special = VMERGE_VVM_F (vy_special, signed_zero, zero_to_zero, vlen);
    }

  vy_special = VMERGE_VVM_F (vy_special, VFMV_VF (1.0f, vlen), x_eq_one,
                             vlen);
  /* Do not re-expose a signaling NaN after the special path has quieted it.  */
  VBOOL_F y_eq_one_non_nan
      = __riscv_vmand (y_eq_one, __riscv_vmnot (x_nan, vlen), vlen);
  vy_special = VMERGE_VVM_F (vy_special, x, y_eq_one_non_nan, vlen);
  vy_special = VMERGE_VVM_F (vy_special, VFMV_VF (1.0f, vlen), y_eq_zero,
                             vlen);

  x_work = VFMERGE_VFM (x, 1.0f, core_skip, vlen);
  y_work = VFMERGE_VFM (y, 0.0f, core_skip, vlen);

  result = FASTPOWF_EXPK3F (
      __riscv_vfmul (FASTPOWF_LOGK3F (__riscv_vfabs (x_work, vlen), vlen), y_work,
                     vlen),
      vlen);

  VBOOL_F sign_flip = __riscv_vmand (x_neg, y_is_odd, vlen);
  VUINT_F sign_mask
      = __riscv_vmerge (VMVU_VX (0, vlen), 0x80000000, sign_flip, vlen);
  result = U_AS_F (VXOR_VV (F_AS_U (result), sign_mask, vlen));

  result = VMERGE_VVM_F (result, vy_special, core_skip, vlen);

  RESTORE_FRM;
  return result;
}

#undef LMUL
#undef SIMDLEN
#undef FASTPOWF_EXPK3F
#undef FASTPOWF_LOGK3F
#undef FASTPOWF_VTRUNCATE
#undef FASTPOWF_VLDEXP2
#undef FASTPOWF_VPOW2I
#undef FASTPOWF_VLDEXP3
#undef FASTPOWF_VILOGB2K
#undef MAKE_VSTORE
#undef VFSTORE
