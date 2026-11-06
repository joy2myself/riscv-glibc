/* Single-precision vector sincosf function.
   sincosf(x) computes both sin(x) and cos(x).

   The finite path follows SLEEF XSINCOSF_U1 (Sleef_sincosf_u10).

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

/* The finite path below follows the pinned SLEEF XSINCOSF_U1 graph.  These
   helpers are target-owned translations of SLEEF's rempif/rempisubf support;
   the clean entry keeps its own classify/prepare/sanitize/merge envelope.  */

/* SLEEF vorsign: preserve X's bits and OR in Y's sign bit.  */
static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (sincosf_or_sign, LMUL)
  (VFLOAT x, VFLOAT y, size_t vlen)
{
  VUINT sign = __riscv_vand (F_AS_U_F (y), 0x80000000u, vlen);
  return U_AS_F_F (__riscv_vor (F_AS_U_F (x), sign, vlen));
}

/* SLEEF rempisubf, including its source comparison and RTZ quotient
   conversion.  */
static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (sincosf_reduce_quarter, LMUL)
  (VFLOAT x, VINT *q_out, size_t vlen)
{
  VFLOAT c = __riscv_vfsgnjx
    (VFMV_VF_F (RVV_SLEEF_TRIG_REMPI_TWO23, vlen), x, vlen);
  VFLOAT x4 = VFMUL_VF (x, RVV_SLEEF_TRIG_REMPI_FOUR, vlen);
  VFLOAT abs_x4 = __riscv_vfabs (x4, vlen);
  VFLOAT abs_x = __riscv_vfabs (x, vlen);
  VBOOL x4_large = __riscv_vmfgt
    (abs_x4, RVV_SLEEF_TRIG_REMPI_TWO23, vlen);
  VBOOL x_large = __riscv_vmfgt
    (abs_x, RVV_SLEEF_TRIG_REMPI_TWO23, vlen);

  VFLOAT rint4x = VFMADD_VV
    (VFMV_VF_F (RVV_SLEEF_TRIG_REMPI_FOUR, vlen), x, c, vlen);
  rint4x = VFSUB_VV (rint4x, c, vlen);
  rint4x = RVV_LIBMVEC_IMPL_NAME (sincosf_or_sign, LMUL)
    (rint4x, x, vlen);
  rint4x = VMERGE_VVM_F (rint4x, x4, x4_large, vlen);

  VFLOAT rintx = VFADD_VV (x, c, vlen);
  rintx = VFSUB_VV (rintx, c, vlen);
  rintx = RVV_LIBMVEC_IMPL_NAME (sincosf_or_sign, LMUL)
    (rintx, x, vlen);
  rintx = VMERGE_VVM_F (rintx, x, x_large, vlen);

  VFLOAT rem = VFMADD_VV
    (VFMV_VF_F (RVV_SLEEF_TRIG_REMPI_NEG_QUARTER, vlen), rint4x, x, vlen);
  *q_out = __riscv_vfcvt_rtz_x
    (VFMADD_VV (VFMV_VF_F (RVV_SLEEF_TRIG_REMPI_NEG_FOUR, vlen), rintx,
                rint4x, vlen), vlen);
  return rem;
}

/* This version computes both sinf and cosf simultaneously.  */

static void
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT vx_orig,
                                           VFLOAT *vsin_out,
                                           VFLOAT *vcos_out, size_t vlen)
{
  VFLOAT vx, vsin, vcos;
  VFLOAT vsin_special, vcos_special;
  VBOOL special_args;

  SET_ROUNDTONEAREST;

  vx = vx_orig;

  /* Special cases: NaN/Inf and |x| very small.  */
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
      /* For tiny |x|, sinf(x) ~= x.  */
      vsin_special = vx_orig;

      /* For tiny |x|, cosf(x) ~= 1. Use 1 + x so nonzero x raises inexact.  */
      vcos_special = VFMV_VF_F (fp_posOne_f, vlen);
      VFLOAT vcos_tiny = __riscv_vfadd (small_x, vx_orig, fp_posOne_f, vlen);
      vcos_special = VMERGE_VVM_F (vcos_special, vcos_tiny, small_x, vlen);

      /* Substitute Inf with sNaN and propagate NaNs via x+x.  */
      VBOOL inf_mask;
      IDENTIFY_F (vclass, class_Inf, inf_mask, vlen);
      VFLOAT vy_NaN_Inf
        = __riscv_vfmerge (vx_orig, fp_sNaN_f, inf_mask, vlen);
      vy_NaN_Inf = __riscv_vfadd (NaN_Inf, vy_NaN_Inf, vy_NaN_Inf, vlen);
      vsin_special = VMERGE_VVM_F (vsin_special, vy_NaN_Inf, NaN_Inf, vlen);
      vcos_special = VMERGE_VVM_F (vcos_special, vy_NaN_Inf, NaN_Inf, vlen);

      /* Mask special lanes out from the main path to avoid spurious
         exceptions (e.g. NaN->int conversion).  */
      vx = VFMERGE_VFM (vx, fp_posZero_f, special_args, vlen);
    }

  /* Pinned SLEEF XSINCOSF_U1 reduction.  Special lanes were sanitized above,
     so this source-shaped finite graph is safe to run over every lane.  */
  VFLOAT u = __riscv_vfcvt_f
    (__riscv_vfcvt_x
       (VFMUL_VF (vx, RVV_SLEEF_TRIG_TWO_OVER_PI, vlen), vlen), vlen);
  VINT q = __riscv_vfcvt_x (u, vlen);

  VFLOAT v = VFMADD_VV
    (u, VFMV_VF_F (RVV_SLEEF_TRIG_NEG_PIBY2_A2_F, vlen), vx, vlen);
  VFLOAT s_hi, s_lo;
  RVV_SLEEF_ORIG_DFADD2_VF_VF
    (v, VFMUL_VF (u, RVV_SLEEF_TRIG_NEG_PIBY2_B2_F, vlen), s_hi, s_lo,
     vlen);
  RVV_SLEEF_ORIG_DFADD_VF2_VF2_VF
    (s_hi, s_lo, VFMUL_VF (u, RVV_SLEEF_TRIG_NEG_PIBY2_C2_F, vlen), s_hi,
     s_lo, vlen);

  VBOOL g = __riscv_vmflt
    (__riscv_vfabs (vx, vlen), RVV_SLEEF_TRIG_RANGE_MAX2, vlen);
  if (__riscv_vcpop (g, vlen) != vlen)
    {
      VFLOAT rem_hi, rem_lo;
      VINT rem_q;
      /* Source rempif reduction, with RVV byte-indexed table loads.  */
      {
        VFLOAT input_a = vx;

        const struct v_sleef_trig_data *trig_data =
          ptr_barrier (&__v_sleef_trig_data);

        /* Extract the normalized exponent before table selection.  */
        VINT ex;
        {
          VUINT bits = F_AS_U_F (input_a);
          VUINT exponent = __riscv_vsrl (bits, MAN_LEN_F, vlen);
          exponent = __riscv_vand (exponent, 0xff, vlen);
          ex = __riscv_vsub (U_AS_I_F (exponent), EXP_BIAS_F, vlen);
        }
        ex = __riscv_vsub (ex, 25, vlen);

        VBOOL scale_mask = __riscv_vmsgt (ex, 90 - 25, vlen);
        VINT scale = VMERGE_VXM_I (VMVI_VX_F (0, vlen), -64, scale_mask, vlen);
        /* Apply the exact exponent-field scaling.  */
        input_a = I_AS_F_F (__riscv_vadd
          (F_AS_I_F (input_a), __riscv_vsll (scale, MAN_LEN_F, vlen), vlen));

        VBOOL ex_negative = __riscv_vmslt (ex, 0, vlen);
        ex = VMERGE_VXM_I (ex, 0, ex_negative, vlen);
        ex = __riscv_vsll (ex, 4, vlen);
        VUINT ex_u = I_AS_U_F (ex);

        VFLOAT tab0 = RVV_VLUXEI32_F (trig_data->rempitabsp, ex_u, vlen);
        VFLOAT tab1 = RVV_VLUXEI32_F (trig_data->rempitabsp + 1, ex_u, vlen);
        VFLOAT tab2 = RVV_VLUXEI32_F (trig_data->rempitabsp + 2, ex_u, vlen);
        VFLOAT tab3 = RVV_VLUXEI32_F (trig_data->rempitabsp + 3, ex_u, vlen);

        VFLOAT x_hi, x_lo;
        RVV_SLEEF_ORIG_DFMUL_VF_VF (input_a, tab0, x_hi, x_lo, vlen);
        VINT q;
        x_hi = RVV_LIBMVEC_IMPL_NAME (sincosf_reduce_quarter, LMUL)
          (x_hi, &q, vlen);
        RVV_SLEEF_ORIG_DFNORMALIZE_VF2 (x_hi, x_lo, x_hi, x_lo, vlen);

        VFLOAT y_hi, y_lo;
        RVV_SLEEF_ORIG_DFMUL_VF_VF (input_a, tab1, y_hi, y_lo, vlen);
        RVV_SLEEF_ORIG_DFADD2_VF2_VF2
          (x_hi, x_lo, y_hi, y_lo, x_hi, x_lo, vlen);
        VINT q_part;
        x_hi = RVV_LIBMVEC_IMPL_NAME (sincosf_reduce_quarter, LMUL)
          (x_hi, &q_part, vlen);
        q = __riscv_vadd (q, q_part, vlen);
        RVV_SLEEF_ORIG_DFNORMALIZE_VF2 (x_hi, x_lo, x_hi, x_lo, vlen);

        y_hi = tab2;
        y_lo = tab3;
        RVV_SLEEF_ORIG_DFMUL_VF2_VF (y_hi, y_lo, input_a, y_hi, y_lo, vlen);
        RVV_SLEEF_ORIG_DFADD2_VF2_VF2
          (x_hi, x_lo, y_hi, y_lo, x_hi, x_lo, vlen);
        RVV_SLEEF_ORIG_DFNORMALIZE_VF2 (x_hi, x_lo, x_hi, x_lo, vlen);
        RVV_SLEEF_ORIG_DFMUL_VF2_VF2
          (x_hi, x_lo,
           VFMV_VF_F ((RVV_PI_F * 2.0f), vlen),
           VFMV_VF_F ((RVV_PI_LO_F * 2.0f), vlen),
           x_hi, x_lo, vlen);

        VBOOL a_small = __riscv_vmflt
          (__riscv_vfabs (input_a, vlen), RVV_SLEEF_TRIG_REMPI_SMALL, vlen);
        x_hi = VMERGE_VVM_F (x_hi, input_a, a_small, vlen);
        x_lo = VFMERGE_VFM (x_lo, fp_posZero_f, a_small, vlen);

        rem_hi = x_hi;
        rem_lo = x_lo;
        rem_q = q;
      }
      VBOOL large = __riscv_vmnot (g, vlen);
      q = VMERGE_VVM_I (q, rem_q, large, vlen);
      s_hi = VMERGE_VVM_F (s_hi, rem_hi, large, vlen);
      s_lo = VMERGE_VVM_F (s_lo, rem_lo, large, vlen);
    }

  /* Source keeps the unreduced pair in T and replaces only S.x with
     dfsqu_vf_vf2(S).  */
  VFLOAT t_hi = s_hi;
  VFLOAT t_lo = s_lo;
  /* Replace only the high component with the source pair square.  */
  {
    VFLOAT xhxl = VFMUL_VV (s_hi, s_lo, vlen);
    s_hi = VFMADD_VV (s_hi, s_hi, VFADD_VV (xhxl, xhxl, vlen), vlen);
  }

  VFLOAT poly = VFMV_VF_F (RVV_SINCOSF_SIN_POLY0_F, vlen);
  poly = VFMADD_VV
    (poly, s_hi, VFMV_VF_F (RVV_SINCOSF_SIN_POLY1_F, vlen), vlen);
  poly = VFMADD_VV
    (poly, s_hi, VFMV_VF_F (RVV_SINCOSF_SIN_POLY2_F, vlen), vlen);
  VFLOAT sin_term = VFMUL_VV (s_hi, t_hi, vlen);
  sin_term = VFMUL_VV (poly, sin_term, vlen);
  VFLOAT sin_hi, sin_lo;
  RVV_SLEEF_ORIG_DFADD_VF2_VF2_VF
    (t_hi, t_lo, sin_term, sin_hi, sin_lo, vlen);
  VFLOAT sin_r = VFADD_VV (sin_hi, sin_lo, vlen);
  VBOOL negzero = __riscv_vmseq
    (F_AS_U_F (vx_orig), 0x80000000u, vlen);
  sin_r = VMERGE_VVM_F (sin_r, VFMV_VF_F (-0.0f, vlen), negzero, vlen);

  poly = VFMV_VF_F (RVV_SINCOSF_COS_POLY0_F, vlen);
  poly = VFMADD_VV
    (poly, s_hi, VFMV_VF_F (RVV_SINCOSF_COS_POLY1_F, vlen), vlen);
  poly = VFMADD_VV
    (poly, s_hi, VFMV_VF_F (RVV_SINCOSF_COS_POLY2_F, vlen), vlen);
  poly = VFMADD_VV
    (poly, s_hi, VFMV_VF_F (RVV_SINCOSF_COS_POLY3_F, vlen), vlen);
  poly = VFMADD_VV
    (poly, s_hi, VFMV_VF_F (RVV_SINCOSF_COS_POLY4_F, vlen), vlen);
  VFLOAT cos_term_hi, cos_term_lo;
  RVV_SLEEF_ORIG_DFMUL_VF_VF
    (s_hi, poly, cos_term_hi, cos_term_lo, vlen);
  VFLOAT cos_hi, cos_lo;
  RVV_SLEEF_ORIG_DFADD_VF_VF2
    (VFMV_VF_F (fp_posOne_f, vlen), cos_term_hi, cos_term_lo,
     cos_hi, cos_lo, vlen);
  VFLOAT cos_r = VFADD_VV (cos_hi, cos_lo, vlen);

  /* Source quadrant reconstruction and bitwise sign changes.  */
  VBOOL even = __riscv_vmseq
    (__riscv_vand (q, 1, vlen), 0, vlen);
  VFLOAT sin_pre = VMERGE_VVM_F (cos_r, sin_r, even, vlen);
  VFLOAT cos_pre = VMERGE_VVM_F (sin_r, cos_r, even, vlen);

  VBOOL sin_flip = __riscv_vmseq
    (__riscv_vand (q, 2, vlen), 2, vlen);
  VUINT sin_sign = __riscv_vmerge
    (VMVU_VX_F (0, vlen), VMVU_VX_F (0x80000000u, vlen), sin_flip, vlen);
  VBOOL cos_flip = __riscv_vmseq
    (__riscv_vand (__riscv_vadd (q, 1, vlen), 2, vlen), 2, vlen);
  VUINT cos_sign = __riscv_vmerge
    (VMVU_VX_F (0, vlen), VMVU_VX_F (0x80000000u, vlen), cos_flip, vlen);
  vsin = U_AS_F_F (__riscv_vxor (F_AS_U_F (sin_pre), sin_sign, vlen));
  vcos = U_AS_F_F (__riscv_vxor (F_AS_U_F (cos_pre), cos_sign, vlen));

  if (nb_special_args > 0)
    {
      vsin = VMERGE_VVM_F (vsin, vsin_special, special_args, vlen);
      vcos = VMERGE_VVM_F (vcos, vcos_special, special_args, vlen);
    }

  *vsin_out = vsin;
  *vcos_out = vcos;

  RESTORE_FRM;
}

#undef LMUL
