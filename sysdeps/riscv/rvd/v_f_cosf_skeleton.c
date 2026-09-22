/* Single-precision vector cos function.
   cosf(x) computes the cosine of x

   Based on SLEEF xcosf_u1 implementation (1-ULP).

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

/* Helpers for SLEEF rempif.  */
/* Shared float-float helpers live in v_math.h.  */

/* This version reduces argument to [-pi/4, pi/4] and computes cos(r)
   by merging the appropriate coefficients into a vector register. */

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  const struct v_sleef_trig_data *trig_data =
    ptr_barrier (&__v_sleef_trig_data);

  VFLOAT vx_orig = x;
  VFLOAT vx, vy, vy_special;
  VBOOL special_args;

  SET_ROUNDTONEAREST;

  vx = vx_orig;

  /* Special cases: NaN/Inf and |x| very small (cosf(x) rounds to 1).  */
  VUINT vclass = __riscv_vfclass (vx, vlen);
  VBOOL NaN_Inf;
  IDENTIFY_F (vclass, class_NaN | class_Inf, NaN_Inf, vlen);

  VFLOAT ax = __riscv_vfabs (vx, vlen);
  VUINT expo_x = __riscv_vsrl (F_AS_U_F (ax), MAN_LEN_F, vlen);
  expo_x = __riscv_vand (expo_x, 0xff, vlen);
  VBOOL small_x = __riscv_vmsleu (expo_x, EXP_BIAS_F - MAN_LEN_F - 5, vlen);
  special_args = __riscv_vmor (NaN_Inf, small_x, vlen);

  UINT nb_special_args = __riscv_vcpop (special_args, vlen);
  if (nb_special_args > 0)
    {
      VBOOL id_mask;
      IDENTIFY_F (vclass, class_Inf, id_mask, vlen);
      /* For tiny |x|, use 1 + x so nonzero x raises inexact.  */
      vy_special = __riscv_vfadd (small_x, vx, fp_posOne_f, vlen);
      /* Substitute Inf with sNaN and trigger invalid for Inf.  */
      VFLOAT vy_NaN_Inf = __riscv_vfmerge (vx, fp_sNaN_f, id_mask, vlen);
      vy_NaN_Inf = __riscv_vfadd (NaN_Inf, vy_NaN_Inf, vy_NaN_Inf, vlen);
      vy_special = VMERGE_VVM_F (vy_special, vy_NaN_Inf, NaN_Inf, vlen);
      /* Mask special lanes out from the main path to avoid spurious
         exceptions (e.g. NaN->int conversion).  */
      vx = VFMERGE_VFM (vx, fp_posZero_f, special_args, vlen);
    }

  /* Small/large argument split for reduction.  */
  VFLOAT abs_vx = __riscv_vfabs (vx, vlen);
  VBOOL g_small = __riscv_vmflt
    (abs_vx, RVV_SLEEF_TRIG_RANGE_MAX2, vlen);
  VBOOL g_large = __riscv_vmnot (g_small, vlen);

  /* For lanes that will use rempif, clamp vx to 0 for the initial dq/q
     computation to avoid out-of-range conversions.  */
  VFLOAT vx_small = VFMERGE_VFM (vx, fp_posZero_f, g_large, vlen);

  /* Range reduction (SLEEF xcosf_u1): dq = 2*rint(vx*(1/pi) - 0.5) + 1.  */
  VFLOAT t = VFMADD_VV
    (vx_small, VFMV_VF_F (RVV_SLEEF_TRIG_ONE_OVER_PI, vlen),
     VFMV_VF_F (RVV_SLEEF_COSF_REDUCE_NEG_HALF, vlen), vlen);
  VINT n = __riscv_vfcvt_x (t, vlen);
  VFLOAT rint_t = __riscv_vfcvt_f (n, vlen);
  VFLOAT dq = VFMADD_VV
    (rint_t, VFMV_VF_F (RVV_SLEEF_COSF_REDUCE_TWO, vlen),
     VFMV_VF_F (fp_posOne_f, vlen), vlen);
  VINT q = __riscv_vfcvt_x (dq, vlen);

  VFLOAT s_hi, s_lo;
  VFLOAT tmp = VFMUL_VF
    (dq, RVV_SLEEF_COSF_NEG_PI_A2_HALF, vlen);
  DFADD2_VF2_VF_VF (vx, tmp, s_hi, s_lo, vlen);
  tmp = VFMUL_VF
    (dq, RVV_SLEEF_COSF_NEG_PI_B2_HALF, vlen);
  DFADD2_VF2_VF2_VF (s_hi, s_lo, tmp, s_hi, s_lo, vlen);
  tmp = VFMUL_VF
    (dq, RVV_SLEEF_COSF_NEG_PI_C2_HALF, vlen);
  DFADD2_VF2_VF2_VF (s_hi, s_lo, tmp, s_hi, s_lo, vlen);

  if (__riscv_vcpop (g_large, vlen) > 0)
    {
      /* rempif implementation (vectorized) to handle large arguments.  */
      VFLOAT a = vx;

      /* ex = ilogb2k(a) - 25 */
      VUINT ia = __riscv_vsrl (F_AS_U_F (a), 23, vlen);
      ia = __riscv_vand (ia, 0xff, vlen);
      VINT ex = __riscv_vsub (U_AS_I_F (ia), 0x7f + 25, vlen);

      /* Scale extremely large values down by 2^(-64) to keep indices bounded
         (SLEEF behavior).  */
      VBOOL ex_gt = __riscv_vmsgt (ex, 90 - 25, vlen);
      VINT scal = VMVI_VX_F (0, vlen);
      scal = VMERGE_VXM_I (scal, -64, ex_gt, vlen);
      VINT a_bits = F_AS_I_F (a);
      a = I_AS_F_F (__riscv_vadd (a_bits, __riscv_vsll (scal, 23, vlen), vlen));

      /* Clamp ex at 0 and convert to byte offsets for vluxei32.
         SLEEF indexes Sleef_rempitabsp in groups of 4 floats; each group
         is 16 bytes.  */
      VBOOL ex_neg = __riscv_vmslt (ex, 0, vlen);
      ex = VMERGE_VXM_I (ex, 0, ex_neg, vlen);
      ex = __riscv_vsll (ex, 4, vlen);
      VUINT ex_u = I_AS_U_F (ex);

      VFLOAT tab0 = RVV_VLUXEI32_F (trig_data->rempitabsp, ex_u, vlen);
      VFLOAT tab1 = RVV_VLUXEI32_F (trig_data->rempitabsp + 1, ex_u, vlen);
      VFLOAT tab2 = RVV_VLUXEI32_F (trig_data->rempitabsp + 2, ex_u, vlen);
      VFLOAT tab3 = RVV_VLUXEI32_F (trig_data->rempitabsp + 3, ex_u, vlen);

      VFLOAT r_hi, r_lo;
      DFMUL_VF2_VF_VF (a, tab0, r_hi, r_lo, vlen);

      /* rempisubf on r_hi */
      VFLOAT c = __riscv_vfsgnjx
        (VFMV_VF_F (RVV_SLEEF_TRIG_REMPI_TWO23, vlen), r_hi, vlen);
      VFLOAT r4 = __riscv_vfmul
        (r_hi, RVV_SLEEF_TRIG_REMPI_FOUR, vlen);
      VFLOAT abs_r4 = __riscv_vfabs (r4, vlen);
      VFLOAT abs_r = __riscv_vfabs (r_hi, vlen);
      VBOOL r4_large = __riscv_vmfgt
        (abs_r4, RVV_SLEEF_TRIG_REMPI_TWO23, vlen);
      VBOOL r_large = __riscv_vmfgt
        (abs_r, RVV_SLEEF_TRIG_REMPI_TWO23, vlen);
      VFLOAT rint4x = VFMADD_VV
        (VFMV_VF_F (RVV_SLEEF_TRIG_REMPI_FOUR, vlen), r_hi, c, vlen);
      rint4x = __riscv_vfsub (rint4x, c, vlen);
      rint4x = U_AS_F_F (__riscv_vor
        (F_AS_U_F (rint4x),
         __riscv_vand (F_AS_U_F (r_hi), 0x80000000u, vlen), vlen));
      rint4x = __riscv_vmerge (rint4x, r4, r4_large, vlen);
      VFLOAT rintx = __riscv_vfsub (__riscv_vfadd (r_hi, c, vlen), c, vlen);
      rintx = U_AS_F_F (__riscv_vor
        (F_AS_U_F (rintx),
         __riscv_vand (F_AS_U_F (r_hi), 0x80000000u, vlen), vlen));
      rintx = __riscv_vmerge (rintx, r_hi, r_large, vlen);
      VFLOAT rem = VFMADD_VV
        (VFMV_VF_F (RVV_SLEEF_TRIG_REMPI_NEG_QUARTER, vlen), rint4x,
         r_hi, vlen);
      VFLOAT qf = __riscv_vfmadd
        (rintx, RVV_SLEEF_TRIG_REMPI_NEG_FOUR, rint4x, vlen);
      VINT q_rem = __riscv_vfcvt_rtz_x (qf, vlen);
      r_hi = rem;
      DFNORMALIZE_VF2_VF2 (r_hi, r_lo, r_hi, r_lo, vlen);

      VFLOAT y_hi, y_lo;
      DFMUL_VF2_VF_VF (a, tab1, y_hi, y_lo, vlen);
      DFADD2_VF2_VF2_VF2 (r_hi, r_lo, y_hi, y_lo, r_hi, r_lo, vlen);

      /* rempisubf on r_hi again */
      c = __riscv_vfsgnjx
        (VFMV_VF_F (RVV_SLEEF_TRIG_REMPI_TWO23, vlen), r_hi, vlen);
      r4 = __riscv_vfmul (r_hi, RVV_SLEEF_TRIG_REMPI_FOUR, vlen);
      abs_r4 = __riscv_vfabs (r4, vlen);
      abs_r = __riscv_vfabs (r_hi, vlen);
      r4_large = __riscv_vmfgt
        (abs_r4, RVV_SLEEF_TRIG_REMPI_TWO23, vlen);
      r_large = __riscv_vmfgt
        (abs_r, RVV_SLEEF_TRIG_REMPI_TWO23, vlen);
      rint4x = VFMADD_VV
        (VFMV_VF_F (RVV_SLEEF_TRIG_REMPI_FOUR, vlen), r_hi, c, vlen);
      rint4x = __riscv_vfsub (rint4x, c, vlen);
      rint4x = U_AS_F_F (__riscv_vor
        (F_AS_U_F (rint4x),
         __riscv_vand (F_AS_U_F (r_hi), 0x80000000u, vlen), vlen));
      rint4x = __riscv_vmerge (rint4x, r4, r4_large, vlen);
      rintx = __riscv_vfsub (__riscv_vfadd (r_hi, c, vlen), c, vlen);
      rintx = U_AS_F_F (__riscv_vor
        (F_AS_U_F (rintx),
         __riscv_vand (F_AS_U_F (r_hi), 0x80000000u, vlen), vlen));
      rintx = __riscv_vmerge (rintx, r_hi, r_large, vlen);
      rem = VFMADD_VV
        (VFMV_VF_F (RVV_SLEEF_TRIG_REMPI_NEG_QUARTER, vlen), rint4x,
         r_hi, vlen);
      qf = __riscv_vfmadd
        (rintx, RVV_SLEEF_TRIG_REMPI_NEG_FOUR, rint4x, vlen);
      q_rem = __riscv_vadd (q_rem, __riscv_vfcvt_rtz_x (qf, vlen), vlen);
      r_hi = rem;
      DFNORMALIZE_VF2_VF2 (r_hi, r_lo, r_hi, r_lo, vlen);

      /* Add the last table term (as double-float) and scale by 2*pi.  */
      y_hi = tab2;
      y_lo = tab3;
      DFMUL_VF2_VF2_VF (y_hi, y_lo, a, y_hi, y_lo, vlen);
      DFADD2_VF2_VF2_VF2 (r_hi, r_lo, y_hi, y_lo, r_hi, r_lo, vlen);
      DFNORMALIZE_VF2_VF2 (r_hi, r_lo, r_hi, r_lo, vlen);
      RVV_SLEEF_ORIG_DFMUL_VF2_VF2
        (r_hi, r_lo,
         VFMV_VF_F (RVV_SLEEF_TRIG_REMPI_TWO_PI_HI, vlen),
         VFMV_VF_F (RVV_SLEEF_TRIG_REMPI_TWO_PI_LO, vlen),
         r_hi, r_lo, vlen);

      /* For small |a|, rempif returns {a, 0}.  */
      VBOOL a_small = __riscv_vmflt
        (__riscv_vfabs (a, vlen), RVV_SLEEF_TRIG_REMPI_SMALL, vlen);
      r_hi = __riscv_vmerge (r_hi, a, a_small, vlen);
      r_lo = VFMERGE_VFM (r_lo, fp_posZero_f, a_small, vlen);

      /* SLEEF xcosf_u1: post-process rempif result into (q2, t).  */
      VINT q2 = __riscv_vand (q_rem, 3, vlen);
      q2 = __riscv_vadd (q2, q2, vlen);
      VBOOL r_gt0 = __riscv_vmfgt (r_hi, fp_posZero_f, vlen);
      VINT add = VMVI_VX_F (7, vlen);
      add = VMERGE_VXM_I (add, 8, r_gt0, vlen);
      q2 = __riscv_vadd (q2, add, vlen);
      q2 = __riscv_vsra (q2, 1, vlen);

      VINT q_and1 = __riscv_vand (q_rem, 1, vlen);
      VBOOL o = __riscv_vmseq (q_and1, 0, vlen);
      VFLOAT y = VFMV_VF_F (fp_negOne_f, vlen);
      y = VFMERGE_VFM (y, fp_posZero_f, r_gt0, vlen);

      VFLOAT corr_hi = __riscv_vfsgnjx
        (VFMV_VF_F (RVV_SLEEF_COSF_REMPI_HALF_PI_HI, vlen), y, vlen);
      VFLOAT corr_lo = __riscv_vfsgnjx
        (VFMV_VF_F (RVV_SLEEF_COSF_REMPI_HALF_PI_LO, vlen), y, vlen);
      VFLOAT t_hi, t_lo;
      DFADD2_VF2_VF2_VF2 (r_hi, r_lo, corr_hi, corr_lo, t_hi, t_lo, vlen);
      r_hi = __riscv_vmerge (r_hi, t_hi, o, vlen);
      r_lo = __riscv_vmerge (r_lo, t_lo, o, vlen);
      DFNORMALIZE_VF2_VF2 (r_hi, r_lo, r_hi, r_lo, vlen);

      q = __riscv_vmerge (q, q2, g_large, vlen);
      s_hi = __riscv_vmerge (s_hi, r_hi, g_large, vlen);
      s_lo = __riscv_vmerge (s_lo, r_lo, g_large, vlen);
    }

  /* Polynomial evaluation (SLEEF xcosf_u1).  */
  VFLOAT t_hi = s_hi, t_lo = s_lo;
  RVV_SLEEF_ORIG_DFSQU_VF2 (s_hi, s_lo, s_hi, s_lo, vlen);

  VFLOAT u = VFMV_VF_F (RVV_SLEEF_COSF_POLY_0, vlen);
  u = __riscv_vfmadd (u, s_hi,
                      VFMV_VF_F (RVV_SLEEF_COSF_POLY_1, vlen), vlen);
  u = __riscv_vfmadd (u, s_hi,
                      VFMV_VF_F (RVV_SLEEF_COSF_POLY_2, vlen), vlen);

  VFLOAT x_hi, x_lo;
  VFLOAT ux = __riscv_vfmul (u, s_hi, vlen);
  RVV_SLEEF_ORIG_DFADD_VF_VF
    (VFMV_VF_F (RVV_SLEEF_COSF_POLY_3, vlen), ux,
     x_hi, x_lo, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (x_hi, x_lo, s_hi, s_lo, x_hi, x_lo, vlen);
  RVV_SLEEF_ORIG_DFADD_VF_VF2
    (VFMV_VF_F (fp_posOne_f, vlen), x_hi, x_lo, x_hi, x_lo, vlen);

  DFMUL_VF_VF2_VF2 (t_hi, t_lo, x_hi, x_lo, vy, vlen);

  VINT q_and2 = __riscv_vand (q, 2, vlen);
  VBOOL flip = __riscv_vmseq (q_and2, 0, vlen);
  VUINT sign_mask = __riscv_vmerge (VMVU_VX_F (0, vlen), 0x80000000, flip, vlen);
  vy = U_AS_F_F (__riscv_vxor (F_AS_U_F (vy), sign_mask, vlen));

  if (nb_special_args > 0)
    vy = VMERGE_VVM_F (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
