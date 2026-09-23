/* Single-precision vector sinf function.
   The finite path follows SLEEF xsinf_u1 (Sleef_sinf_u10).

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
  VFLOAT vx_orig = x;
  VFLOAT vx, vy, vy_special;
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
  VBOOL small_x = __riscv_vmsleu
    (expo_x, EXP_BIAS_F - MAN_LEN_F - 5, vlen);
  special_args = __riscv_vmor (NaN_Inf, small_x, vlen);

  UINT nb_special_args = __riscv_vcpop (special_args, vlen);
  if (nb_special_args > 0)
    {
      /* For tiny |x|, sinf(x) ~= x.  */
      vy_special = vx_orig;

      /* Scalar glibc forces underflow for nonzero subnormal inputs even
         though the final result is x.  Preserve that side effect locally.  */
      VBOOL subnormal_x = __riscv_vmand
        (small_x, __riscv_vmfne (ax, fp_posZero_f, vlen), vlen);
      if (__riscv_vcpop (subnormal_x, vlen) > 0)
        {
          VFLOAT force_input = VFMERGE_VFM
            (vx_orig, fp_posZero_f,
             __riscv_vmnot (subnormal_x, vlen), vlen);
          volatile VFLOAT force_underflow =
            __riscv_vfmul (force_input, force_input, vlen);
          (void) force_underflow;
        }

      /* Substitute Inf with sNaN and propagate NaNs via x+x.  */
      VBOOL inf_mask;
      IDENTIFY_F (vclass, class_Inf, inf_mask, vlen);
      VFLOAT vy_NaN_Inf = __riscv_vfmerge
        (vx_orig, fp_sNaN_f, inf_mask, vlen);
      vy_NaN_Inf = __riscv_vfadd
        (NaN_Inf, vy_NaN_Inf, vy_NaN_Inf, vlen);
      vy_special = VMERGE_VVM_F
        (vy_special, vy_NaN_Inf, NaN_Inf, vlen);

      /* Mask special lanes out from the main path to avoid spurious
         exceptions (e.g. NaN->int conversion).  */
      vx = VFMERGE_VFM (vx, fp_posZero_f, special_args, vlen);
    }

  /* Direct SLEEF xsinf_u1 finite graph.  */
  VFLOAT rounded_product = RVV_LIBMVEC_IMPL_NAME
    (rvv_vrint_f32, LMUL)
    (VFMUL_VF (vx, RVV_ONE_OVER_PI_F, vlen), vlen);
  VINT q = RVV_VFCVT_X_F_RNE (rounded_product, vlen);
  VFLOAT u = __riscv_vfcvt_f (q, vlen);

  VFLOAT v = VFMADD_VV
    (u, VFMV_VF_F (RVV_SLEEF_TRIG_NEG_PI_A2_F, vlen), vx, vlen);
  VFLOAT s_hi, s_lo;
  RVV_SLEEF_ORIG_DFADD2_VF_VF
    (v, VFMUL_VF (u, RVV_SLEEF_TRIG_NEG_PI_B2_F, vlen), s_hi, s_lo, vlen);
  RVV_SLEEF_ORIG_DFADD_VF2_VF2_VF
    (s_hi, s_lo, VFMUL_VF (u, RVV_SLEEF_TRIG_NEG_PI_C2_F, vlen),
     s_hi, s_lo, vlen);

  VBOOL g = __riscv_vmflt
    (__riscv_vfabs (vx, vlen), RVV_SLEEF_TRIG_RANGE_MAX2, vlen);
  if (__riscv_vcpop (g, vlen) != vlen)
    {
      VFLOAT rem_hi, rem_lo;
      VINT rem_q;
      /* Source-shaped large-argument reduction with exact exponent adapters.  */
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
        ex = __riscv_vsub (ex, RVV_SLEEF_SINF_REMPI_EXPONENT_BIAS, vlen);

        VBOOL scale_mask = __riscv_vmsgt
          (ex, RVV_SLEEF_SINF_REMPI_SCALE_LIMIT, vlen);
        VINT scale = VMERGE_VXM_I
          (VMVI_VX_F (0, vlen), RVV_SLEEF_SINF_REMPI_SCALE, scale_mask, vlen);
        /* Apply the exact exponent-field scaling.  */
        input_a = I_AS_F_F (__riscv_vadd
          (F_AS_I_F (input_a), __riscv_vsll (scale, MAN_LEN_F, vlen), vlen));

        VBOOL ex_negative = __riscv_vmslt (ex, 0, vlen);
        ex = VMERGE_VXM_I (ex, 0, ex_negative, vlen);
        ex = __riscv_vsll (ex, RVV_SLEEF_SINF_REMPI_INDEX_SHIFT, vlen);
        VUINT ex_u = I_AS_U_F (ex);

        VFLOAT tab0 = RVV_VLUXEI32_F (trig_data->rempitabsp, ex_u, vlen);
        VFLOAT tab1 = RVV_VLUXEI32_F (trig_data->rempitabsp + 1, ex_u, vlen);
        VFLOAT tab2 = RVV_VLUXEI32_F (trig_data->rempitabsp + 2, ex_u, vlen);
        VFLOAT tab3 = RVV_VLUXEI32_F (trig_data->rempitabsp + 3, ex_u, vlen);

        VFLOAT x_hi, x_lo;
        RVV_SLEEF_ORIG_DFMUL_VF_VF (input_a, tab0, x_hi, x_lo, vlen);
        VINT q;
        x_hi = RVV_LIBMVEC_IMPL_NAME (sinf_reduce_quarter, LMUL) (x_hi, &q, vlen);
        RVV_SLEEF_ORIG_DFNORMALIZE_VF2 (x_hi, x_lo, x_hi, x_lo, vlen);

        VFLOAT y_hi, y_lo;
        RVV_SLEEF_ORIG_DFMUL_VF_VF (input_a, tab1, y_hi, y_lo, vlen);
        RVV_SLEEF_ORIG_DFADD2_VF2_VF2
          (x_hi, x_lo, y_hi, y_lo, x_hi, x_lo, vlen);
        VINT q_part;
        x_hi = RVV_LIBMVEC_IMPL_NAME (sinf_reduce_quarter, LMUL) (x_hi, &q_part, vlen);
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

      VINT q2 = __riscv_vand (rem_q, 3, vlen);
      VBOOL rem_positive = __riscv_vmfgt
        (rem_hi, fp_posZero_f, vlen);
      VINT q2_add = VMERGE_VXM_I
        (VMVI_VX_F (1, vlen), 2, rem_positive, vlen);
      q2 = __riscv_vadd (__riscv_vadd (q2, q2, vlen), q2_add, vlen);
      q2 = __riscv_vsra (q2, 2, vlen);

      VBOOL rem_q_odd = __riscv_vmseq
        (__riscv_vand (rem_q, 1, vlen), 1, vlen);
      VFLOAT corr_hi = __riscv_vfsgnjx
        (VFMV_VF_F ((-RVV_PIBY2_F), vlen), rem_hi, vlen);
      VFLOAT corr_lo = __riscv_vfsgnjx
        (VFMV_VF_F ((-RVV_PIBY2_LO_F), vlen), rem_hi, vlen);
      VFLOAT alt_hi, alt_lo;
      RVV_SLEEF_ORIG_DFADD2_VF2_VF2
        (rem_hi, rem_lo, corr_hi, corr_lo, alt_hi, alt_lo, vlen);
      VFLOAT selected_hi = VMERGE_VVM_F
        (rem_hi, alt_hi, rem_q_odd, vlen);
      VFLOAT selected_lo = VMERGE_VVM_F
        (rem_lo, alt_lo, rem_q_odd, vlen);
      RVV_SLEEF_ORIG_DFNORMALIZE_VF2
        (selected_hi, selected_lo, selected_hi, selected_lo, vlen);

      VBOOL large = __riscv_vmnot (g, vlen);
      q = __riscv_vmerge (q, q2, large, vlen);
      s_hi = VMERGE_VVM_F (s_hi, selected_hi, large, vlen);
      s_lo = VMERGE_VVM_F (s_lo, selected_lo, large, vlen);
    }

  VFLOAT t_hi = s_hi;
  VFLOAT t_lo = s_lo;
  VFLOAT s2_hi, s2_lo;
  RVV_SLEEF_ORIG_DFSQU_VF2
    (s_hi, s_lo, s2_hi, s2_lo, vlen);

  u = VFMV_VF_F (RVV_SLEEF_TRIG_SIN_POLY0_F, vlen);
  u = VFMADD_VV
    (u, s2_hi, VFMV_VF_F (RVV_SLEEF_TRIG_SIN_POLY1_F, vlen), vlen);
  u = VFMADD_VV
    (u, s2_hi, VFMV_VF_F (RVV_SLEEF_TRIG_SIN_POLY2_F, vlen), vlen);

  VFLOAT inner_hi, inner_lo;
  RVV_SLEEF_ORIG_DFADD_VF_VF
    (VFMV_VF_F (RVV_SLEEF_TRIG_SIN_POLY3_F, vlen),
     VFMUL_VV (u, s2_hi, vlen), inner_hi, inner_lo, vlen);
  VFLOAT prod_hi, prod_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (inner_hi, inner_lo, s2_hi, s2_lo, prod_hi, prod_lo, vlen);
  VFLOAT x_hi, x_lo;
  RVV_SLEEF_ORIG_DFADD_VF_VF2
    (VFMV_VF_F (1.0f, vlen), prod_hi, prod_lo, x_hi, x_lo, vlen);

  VFLOAT result;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2_SCALAR
    (t_hi, t_lo, x_hi, x_lo, result, vlen);
  VBOOL q_odd = __riscv_vmseq
    (__riscv_vand (q, 1, vlen), 1, vlen);
  VUINT sign = __riscv_vmerge
    (VMVU_VX_F (0, vlen), VMVU_VX_F (0x80000000u, vlen), q_odd, vlen);
  vy = U_AS_F_F (__riscv_vxor (F_AS_U_F (result), sign, vlen));

  if (nb_special_args > 0)
    vy = VMERGE_VVM_F (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
