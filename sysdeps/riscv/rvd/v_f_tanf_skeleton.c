/* Single-precision vector tanf function.
   The finite path follows SLEEF xtanf_u1 (Sleef_tanf_u10).

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
  const struct v_sleef_trig_data *data =
    ptr_barrier (&__v_sleef_trig_data);

  VFLOAT vx_orig = x;
  VFLOAT d, vy, vy_special;
  VBOOL special_args;

  SET_ROUNDTONEAREST;

  VUINT vclass = __riscv_vfclass (vx_orig, vlen);
  VBOOL NaN_Inf;
  IDENTIFY_F (vclass, class_NaN | class_Inf, NaN_Inf, vlen);

  VFLOAT abs_x = __riscv_vfabs (vx_orig, vlen);
  VUINT expo_x = __riscv_vsrl (F_AS_U_F (abs_x), MAN_LEN_F, vlen);
  expo_x = __riscv_vand (expo_x, 0xff, vlen);
  VBOOL small_x = __riscv_vmsleu
    (expo_x, EXP_BIAS_F - MAN_LEN_F - 5, vlen);
  special_args = __riscv_vmor (NaN_Inf, small_x, vlen);

  /* tanf(x) rounds to x for tiny finite x.  The scaled product preserves the
     scalar inexact/underflow side effect without sending these lanes through
     integer conversion and argument reduction.  */
  vy_special = __riscv_vfmadd
    (small_x, vx_orig, 0x1.0p-60f, vx_orig, vlen);

  d = vx_orig;
  if (__riscv_vcpop (special_args, vlen) > 0)
    {
      VBOOL id_mask;
      IDENTIFY_F (vclass, class_Inf, id_mask, vlen);
      VFLOAT vy_nan = __riscv_vfmerge (vx_orig, fp_sNaN_f, id_mask, vlen);
      vy_nan = __riscv_vfadd (NaN_Inf, vy_nan, vy_nan, vlen);
      vy_special = VMERGE_VVM_F (vy_special, vy_nan, NaN_Inf, vlen);
      /* Keep NaN/Inf and tiny lanes out of the finite reduction graph.  */
      d = VFMERGE_VFM (d, fp_posZero_f, special_args, vlen);
    }

  VFLOAT abs_d = __riscv_vfabs (d, vlen);
  VBOOL g_small = __riscv_vmflt
    (abs_d, RVV_SLEEF_TRIG_RANGE_MAX2, vlen);
  VBOOL g_large = __riscv_vmnot (g_small, vlen);

  /* rempif handles the large lanes.  Quarantine them before the initial
     reduction so no discarded lane can enter an out-of-range conversion.  */
  VFLOAT d_small = VFMERGE_VFM (d, fp_posZero_f, g_large, vlen);

  VFLOAT rounded_t = __riscv_vfcvt_f
    (RVV_VFCVT_X_F_RNE
       (VFMUL_VF (d_small, RVV_SLEEF_TRIG_TWO_OVER_PI, vlen), vlen), vlen);
  VINT q = RVV_VFCVT_X_F_RNE (rounded_t, vlen);
  VFLOAT q_flt = rounded_t;

  VFLOAT v = __riscv_vfmadd
    (q_flt, VFMV_VF_F (RVV_SLEEF_TRIG_NEG_PIBY2_A2_F, vlen), d_small, vlen);
  VFLOAT s_hi;
  VFLOAT s_lo;
  VFLOAT tmp = __riscv_vfmul
    (q_flt, VFMV_VF_F (RVV_SLEEF_TRIG_NEG_PIBY2_B2_F, vlen), vlen);
  DFADD2_VF2_VF_VF (v, tmp, s_hi, s_lo, vlen);
  tmp = __riscv_vfmul
    (q_flt, VFMV_VF_F (RVV_SLEEF_TRIG_NEG_PIBY2_C2_F, vlen), vlen);
  DFADD_VF2_VF2_VF (s_hi, s_lo, tmp, s_hi, s_lo, vlen);

  if (__riscv_vcpop (g_large, vlen) > 0)
    {
      /* rempif implementation (vectorized) to handle large arguments.  */
      VFLOAT a = __riscv_vfmerge (d, fp_posZero_f, g_small, vlen);

      /* ex = ilogb2k(a) - 25 */
      VUINT ia = __riscv_vsrl (F_AS_U_F (a), 23, vlen);
      ia = __riscv_vand (ia, 0xff, vlen);
      VINT ex = __riscv_vsub (U_AS_I_F (ia), 0x7f + 25, vlen);

      /* Scale extremely large values down by 2^-64 before multiplication.
         Table selection still uses the original exponent below.  */
      VBOOL ex_gt = __riscv_vmsgt (ex, 90 - 25, vlen);
      VINT scal = VMVI_VX_F (0, vlen);
      scal = VMERGE_VXM_I (scal, -64, ex_gt, vlen);
      VINT a_bits = F_AS_I_F (a);
      a = I_AS_F_F (
          __riscv_vadd (a_bits, __riscv_vsll (scal, 23, vlen), vlen));

      /* Clamp ex at 0 and convert to byte offsets for vluxei32.
         SLEEF indexes Sleef_rempitabsp in groups of 4 floats; each group
         is 16 bytes.  */
      VBOOL ex_neg = __riscv_vmslt (ex, 0, vlen);
      ex = VMERGE_VXM_I (ex, 0, ex_neg, vlen);
      ex = __riscv_vsll (ex, 4, vlen);
      VUINT ex_u = I_AS_U_F (ex);

      VFLOAT tab0 = RVV_VLUXEI32_F (data->rempitabsp, ex_u, vlen);
      VFLOAT tab1 = RVV_VLUXEI32_F (data->rempitabsp + 1, ex_u, vlen);
      VFLOAT tab2 = RVV_VLUXEI32_F (data->rempitabsp + 2, ex_u, vlen);
      VFLOAT tab3 = RVV_VLUXEI32_F (data->rempitabsp + 3, ex_u, vlen);

      VFLOAT r_hi, r_lo;
      DFMUL_VF2_VF_VF (a, tab0, r_hi, r_lo, vlen);

      /* Source rempisubf: FMA reconstruction and RTZ quotient conversion.  */
      VINT q_rem;
      r_hi = RVV_LIBMVEC_IMPL_NAME (sinf_reduce_quarter, LMUL)
        (r_hi, &q_rem, vlen);
      DFNORMALIZE_VF2_VF2 (r_hi, r_lo, r_hi, r_lo, vlen);

      VFLOAT y_hi, y_lo;
      DFMUL_VF2_VF_VF (a, tab1, y_hi, y_lo, vlen);
      DFADD2_VF2_VF2_VF2 (r_hi, r_lo, y_hi, y_lo, r_hi, r_lo, vlen);

      /* The second source rempisubf is applied after the double-float add.  */
      VINT q_part;
      r_hi = RVV_LIBMVEC_IMPL_NAME (sinf_reduce_quarter, LMUL)
        (r_hi, &q_part, vlen);
      q_rem = __riscv_vadd (q_rem, q_part, vlen);
      DFNORMALIZE_VF2_VF2 (r_hi, r_lo, r_hi, r_lo, vlen);

      /* Add the last table term (as double-float) and scale by 2*pi.  */
      y_hi = tab2;
      y_lo = tab3;
      DFMUL_VF2_VF2_VF (y_hi, y_lo, a, y_hi, y_lo, vlen);
      DFADD2_VF2_VF2_VF2 (r_hi, r_lo, y_hi, y_lo, r_hi, r_lo, vlen);
      DFNORMALIZE_VF2_VF2 (r_hi, r_lo, r_hi, r_lo, vlen);
      DFMUL_VF2_VF2_VF2 (r_hi, r_lo,
                         VFMV_VF_F (3.1415927410125732422f * 2.0f, vlen),
                         VFMV_VF_F (-8.7422776573475857731e-08f * 2.0f, vlen),
                         r_hi, r_lo, vlen);

      /* For small |a|, rempif returns {a, 0}.  */
      VBOOL a_small = __riscv_vmflt (__riscv_vfabs (a, vlen), 0.7f, vlen);
      r_hi = __riscv_vmerge (r_hi, a, a_small, vlen);
      r_lo = __riscv_vfmerge (r_lo, fp_posZero_f, a_small, vlen);

      q = __riscv_vmerge (q, q_rem, g_large, vlen);
      s_hi = __riscv_vmerge (s_hi, r_hi, g_large, vlen);
      s_lo = __riscv_vmerge (s_lo, r_lo, g_large, vlen);
    }

  VUINT q_and_1 = __riscv_vand (I_AS_U_F (q), 0x1, vlen);
  VBOOL odd = __riscv_vmsne (q_and_1, 0, vlen);
  VUINT sign_mask =
      __riscv_vmerge (VMVU_VX_F (0, vlen), 0x80000000, odd, vlen);
  s_hi = U_AS_F_F (__riscv_vxor (F_AS_U_F (s_hi), sign_mask, vlen));
  s_lo = U_AS_F_F (__riscv_vxor (F_AS_U_F (s_lo), sign_mask, vlen));

  VFLOAT t_hi = s_hi;
  VFLOAT t_lo = s_lo;

  DFSQU_VF2_VF2 (s_hi, s_lo, s_hi, s_lo, vlen);
  DFNORMALIZE_VF2_VF2 (s_hi, s_lo, s_hi, s_lo, vlen);

  VFLOAT u = VFMV_VF_F (RVV_TANF_POLY0, vlen);
  u = __riscv_vfmadd (u, s_hi,
                      VFMV_VF_F (RVV_TANF_POLY1, vlen),
                      vlen);
  u = __riscv_vfmadd (u, s_hi,
                      VFMV_VF_F (RVV_TANF_POLY2, vlen), vlen);
  u = __riscv_vfmadd (u, s_hi,
                      VFMV_VF_F (RVV_TANF_POLY3, vlen), vlen);
  u = __riscv_vfmadd (u, s_hi,
                      VFMV_VF_F (RVV_TANF_POLY4, vlen), vlen);

  VFLOAT x_hi;
  VFLOAT x_lo;
  VFLOAT ux = __riscv_vfmul (u, s_hi, vlen);
  DFADD2_VF2_VF_VF (VFMV_VF_F (RVV_TANF_DF_C0, vlen), ux,
                    x_hi, x_lo, vlen);

  VFLOAT t1_hi;
  VFLOAT t1_lo;
  VFLOAT t2_hi;
  VFLOAT t2_lo;
  DFMUL_VF2_VF2_VF2 (s_hi, s_lo, x_hi, x_lo, t1_hi, t1_lo, vlen);
  DFADD_VF2_VF_VF2 (VFMV_VF_F (RVV_TANF_DF_C1, vlen), t1_hi,
                    t1_lo, t2_hi, t2_lo, vlen);
  DFMUL_VF2_VF2_VF2 (t2_hi, t2_lo, s_hi, s_lo, t1_hi, t1_lo, vlen);
  DFADD_VF2_VF_VF2 (VFMV_VF_F (1.0f, vlen), t1_hi, t1_lo, x_hi, x_lo,
                    vlen);
  DFMUL_VF2_VF2_VF2 (t_hi, t_lo, x_hi, x_lo, x_hi, x_lo, vlen);

  VFLOAT rec_hi;
  VFLOAT rec_lo;
  VFLOAT rec_in_hi = VFMERGE_VFM (x_hi, 1.0f, special_args, vlen);
  VFLOAT rec_in_lo = VFMERGE_VFM (x_lo, 0.0f, special_args, vlen);
  DFREC_VF2_VF2 (rec_in_hi, rec_in_lo, rec_hi, rec_lo, vlen);
  VBOOL recip_lanes
      = __riscv_vmand (odd, __riscv_vmnot (special_args, vlen), vlen);
  x_hi = __riscv_vmerge (x_hi, rec_hi, recip_lanes, vlen);
  x_lo = __riscv_vmerge (x_lo, rec_lo, recip_lanes, vlen);

  VFLOAT y = __riscv_vfadd (x_hi, x_lo, vlen);

  VBOOL is_negzero = __riscv_vmseq (F_AS_U_F (vx_orig),
                                      VMVU_VX_F (0x80000000, vlen), vlen);
  y = __riscv_vmerge (y, vx_orig, is_negzero, vlen);

  vy = VMERGE_VVM_F (y, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
