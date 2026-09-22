/* Double-precision vector hypot function.
   hypot(x, y) computes sqrt(x^2 + y^2) without overflow

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

/* The generic DF helpers in v_math.h are deliberately convenient shared
   building blocks, but their reciprocal and square-root correction trees do
   not preserve SLEEF xhypot_u05's source operation order.  Keep a
   routine-local translation of the source FMA branch here.  The public
   function below still retains the clean EH quarantine and final merge.  */
/* SLEEF ddadd2_vd2_vd2_vd2, common/dd.h:166--170.  */
static inline void
RVV_LIBMVEC_IMPL_NAME (rvv_sleef_clean_hypot_ddadd2_pair_pair, LMUL)
  (VFLOAT x_hi, VFLOAT x_lo, VFLOAT y_hi, VFLOAT y_lo,
   VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT s = VFADD_VV (x_hi, y_hi, vlen);
  VFLOAT v = VFSUB_VV (s, x_hi, vlen);
  VFLOAT t = VFADD_VV
    (VFSUB_VV (x_hi, VFSUB_VV (s, v, vlen), vlen),
     VFSUB_VV (y_hi, v, vlen), vlen);
  *r_hi = s;
  *r_lo = VFADD_VV (t, VFADD_VV (x_lo, y_lo, vlen), vlen);
}

/* SLEEF ddadd2_vd2_vd, common/dd.h:140--144.  */
static inline void
RVV_LIBMVEC_IMPL_NAME (rvv_sleef_clean_hypot_ddadd2_pair_scalar, LMUL)
  (VFLOAT x_hi, VFLOAT x_lo, VFLOAT y, VFLOAT *r_hi, VFLOAT *r_lo,
   size_t vlen)
{
  VFLOAT s = VFADD_VV (x_hi, y, vlen);
  VFLOAT v = VFSUB_VV (s, x_hi, vlen);
  VFLOAT t = VFADD_VV
    (VFSUB_VV (x_hi, VFSUB_VV (s, v, vlen), vlen),
     VFSUB_VV (y, v, vlen), vlen);
  *r_hi = s;
  *r_lo = VFADD_VV (t, x_lo, vlen);
}

/* SLEEF ddmul_vd_vd, common/dd.h:199--202.  */
static inline void
RVV_LIBMVEC_IMPL_NAME (rvv_sleef_clean_hypot_ddmul_scalar_scalar, LMUL)
  (VFLOAT x, VFLOAT y, VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT s = VFMUL_VV (x, y, vlen);
  *r_hi = s;
  *r_lo = VFMSUB_VV (x, y, s, vlen);
}

/* SLEEF ddmul_vd2_vd2_vd2, common/dd.h:209--212.  */
static inline void
RVV_LIBMVEC_IMPL_NAME (rvv_sleef_clean_hypot_ddmul_pair_pair, LMUL)
  (VFLOAT x_hi, VFLOAT x_lo, VFLOAT y_hi, VFLOAT y_lo,
   VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT s = VFMUL_VV (x_hi, y_hi, vlen);
  VFLOAT t = VFMSUB_VV (x_hi, y_hi, s, vlen);
  t = VFMADD_VV (x_lo, y_hi, t, vlen);
  t = VFMADD_VV (x_hi, y_lo, t, vlen);
  *r_hi = s;
  *r_lo = t;
}

/* SLEEF ddmul_vd2_vd, common/dd.h:222--225.  */
static inline void
RVV_LIBMVEC_IMPL_NAME (rvv_sleef_clean_hypot_ddmul_pair_scalar, LMUL)
  (VFLOAT x_hi, VFLOAT x_lo, VFLOAT y, VFLOAT *r_hi, VFLOAT *r_lo,
   size_t vlen)
{
  VFLOAT s = VFMUL_VV (x_hi, y, vlen);
  VFLOAT t = VFMSUB_VV (x_hi, y, s, vlen);
  t = VFMADD_VV (x_lo, y, t, vlen);
  *r_hi = s;
  *r_lo = t;
}

/* SLEEF ddrec_vd2_vd, common/dd.h:227--230.  */
static inline void
RVV_LIBMVEC_IMPL_NAME (rvv_sleef_clean_hypot_ddrec_scalar, LMUL)
  (VFLOAT d, VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT t = VFDIV_VV (VFMV_VF (1.0, vlen), d, vlen);
  VFLOAT u = VFNMSUB_VV (d, t, VFMV_VF (1.0, vlen), vlen);
  *r_hi = t;
  *r_lo = VFMUL_VV (t, u, vlen);
}

/* SLEEF dddiv_vd2_vd2_vd2, common/dd.h:191--197.  */
static inline void
RVV_LIBMVEC_IMPL_NAME (rvv_sleef_clean_hypot_dddiv_pair_pair, LMUL)
  (VFLOAT n_hi, VFLOAT n_lo, VFLOAT d_hi, VFLOAT d_lo,
   VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT t = VFDIV_VV (VFMV_VF (1.0, vlen), d_hi, vlen);
  VFLOAT s = VFMUL_VV (n_hi, t, vlen);
  VFLOAT u = VFMSUB_VV (t, n_hi, s, vlen);
  VFLOAT v = VFNMSUB_VV
    (d_lo, t, VFNMSUB_VV
       (d_hi, t, VFMV_VF (1.0, vlen), vlen), vlen);
  VFLOAT correction = VFMADD_VV (n_lo, t, u, vlen);
  *r_hi = s;
  *r_lo = VFMADD_VV (s, v, correction, vlen);
}

/* SLEEF ddsqu_vd2_vd2, common/dd.h:204--207.  */
static inline void
RVV_LIBMVEC_IMPL_NAME (rvv_sleef_clean_hypot_ddsqu_pair, LMUL)
  (VFLOAT x_hi, VFLOAT x_lo, VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT s = VFMUL_VV (x_hi, x_hi, vlen);
  VFLOAT t = VFMSUB_VV (x_hi, x_hi, s, vlen);
  t = VFMADD_VV (VFADD_VV (x_hi, x_hi, vlen), x_lo, t, vlen);
  *r_hi = s;
  *r_lo = t;
}

/* SLEEF ddsqrt_vd2_vd2, common/dd.h:312--315.  In particular, retain the
   low part of t*t when forming d + t*t before multiplying by 1/t.  */
static inline void
RVV_LIBMVEC_IMPL_NAME (rvv_sleef_clean_hypot_ddsqrt_pair, LMUL)
  (VFLOAT d_hi, VFLOAT d_lo, VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT t = __riscv_vfsqrt (VFADD_VV (d_hi, d_lo, vlen), vlen);

  VFLOAT tsq_hi, tsq_lo;
  RVV_LIBMVEC_IMPL_NAME (rvv_sleef_clean_hypot_ddmul_scalar_scalar, LMUL)
    (t, t, &tsq_hi, &tsq_lo, vlen);

  VFLOAT sum_hi, sum_lo;
  RVV_LIBMVEC_IMPL_NAME (rvv_sleef_clean_hypot_ddadd2_pair_pair, LMUL)
    (d_hi, d_lo, tsq_hi, tsq_lo, &sum_hi, &sum_lo, vlen);

  VFLOAT rec_hi, rec_lo;
  RVV_LIBMVEC_IMPL_NAME (rvv_sleef_clean_hypot_ddrec_scalar, LMUL)
    (t, &rec_hi, &rec_lo, vlen);

  VFLOAT q_hi, q_lo;
  RVV_LIBMVEC_IMPL_NAME (rvv_sleef_clean_hypot_ddmul_pair_pair, LMUL)
    (sum_hi, sum_lo, rec_hi, rec_lo, &q_hi, &q_lo, vlen);

  *r_hi = VFMUL_VF (q_hi, 0.5, vlen);
  *r_lo = VFMUL_VF (q_lo, 0.5, vlen);
}

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT vx_orig, VFLOAT vy_orig,
                                           size_t vlen)
{
  VFLOAT vx, vy, vz, vz_special;
  VBOOL special_args;

  SET_ROUNDTONEAREST;
  vx = vx_orig;
  vy = vy_orig;

  /* Handle NaN/Inf early to match scalar libm behavior.
     For lanes with Inf and no sNaN inputs, hypot(x, y) returns +Inf even if
     the other operand is a qNaN.  Otherwise, propagate via x + y.  */
  VUINT vclass_x = __riscv_vfclass (vx_orig, vlen);
  VUINT vclass_y = __riscv_vfclass (vy_orig, vlen);
  VBOOL special_x, special_y;
  IDENTIFY (vclass_x, class_NaN | class_Inf, special_x, vlen);
  IDENTIFY (vclass_y, class_NaN | class_Inf, special_y, vlen);
  special_args = __riscv_vmor (special_x, special_y, vlen);

  UINT nb_special_args = __riscv_vcpop (special_args, vlen);
  vz_special = VFMV_VF (fp_posZero, vlen);
  if (nb_special_args > 0)
    {
      /* Default: propagate NaNs and raise invalid for sNaN.  */
      vz_special = VFADD_VV_M (special_args, vx_orig, vy_orig, vlen);

      VBOOL x_snan, y_snan, any_snan;
      IDENTIFY (vclass_x, class_sNaN, x_snan, vlen);
      IDENTIFY (vclass_y, class_sNaN, y_snan, vlen);
      any_snan = __riscv_vmor (x_snan, y_snan, vlen);

      VBOOL x_inf, y_inf, any_inf;
      IDENTIFY (vclass_x, class_Inf, x_inf, vlen);
      IDENTIFY (vclass_y, class_Inf, y_inf, vlen);
      any_inf = __riscv_vmor (x_inf, y_inf, vlen);

      VBOOL inf_no_snan
          = __riscv_vmand (any_inf, __riscv_vmnot (any_snan, vlen), vlen);
      vz_special = VFMERGE_VFM (vz_special, fp_posInf, inf_no_snan, vlen);

      /* Sanitize lanes that will be overridden later, avoiding 0/0 in the
         ratio core.  */
      vx = VFMERGE_VFM (vx, fp_posOne, special_args, vlen);
      vy = VFMERGE_VFM (vy, fp_posZero, special_args, vlen);
    }

  /* Take absolute values and order them.  */
  VFLOAT ax = __riscv_vfabs (vx, vlen);
  VFLOAT ay = __riscv_vfabs (vy, vlen);
  VFLOAT vmin = __riscv_vfmin (ax, ay, vlen);
  VFLOAT vmax = __riscv_vfmax (ax, ay, vlen);

  VFLOAT d = vmax;
  VBOOL all_zero = __riscv_vmfeq (vmax, fp_posZero, vlen);
  d = VFMERGE_VFM (d, fp_posOne, all_zero, vlen);

  /* Scale both operands before entering the double-double graph.  At the
     low end this prevents the final low product from underflowing before it
     is combined with a normal high result.  At the high end it prevents the
     reciprocal of max from underflowing and lets the final rescale raise
     only the public overflow/inexact pair when the result really overflows.
     The scale factors are exact powers of two, so they do not themselves
     add an exception.  */
  VBOOL is_tiny
    = __riscv_vmflt (vmax, RVV_HYPOT_CLEAN_SCALE_UP_LIMIT_D, vlen);
  VBOOL is_large
    = __riscv_vmfge (vmax, RVV_HYPOT_CLEAN_SCALE_DOWN_LIMIT_D, vlen);
  VBOOL scale_mask = __riscv_vmor (is_tiny, is_large, vlen);
  UINT nb_scaled = __riscv_vcpop (scale_mask, vlen);

  VFLOAT core_scale = VFMV_VF (1.0, vlen);
  core_scale = VFMERGE_VFM (core_scale, RVV_HYPOT_CLEAN_SCALE_UP_D, is_tiny,
                            vlen);
  core_scale = VFMERGE_VFM (core_scale, RVV_HYPOT_CLEAN_SCALE_DOWN_D, is_large,
                            vlen);

  VFLOAT n = __riscv_vfmul (vmin, core_scale, vlen);
  d = __riscv_vfmul (d, core_scale, vlen);
  VFLOAT vmax_core = __riscv_vfmul (vmax, core_scale, vlen);

  VFLOAT t_hi, t_lo;
  RVV_LIBMVEC_IMPL_NAME (rvv_sleef_clean_hypot_dddiv_pair_pair, LMUL)
    (n, VFMV_VF (fp_posZero, vlen), d, VFMV_VF (fp_posZero, vlen),
     &t_hi, &t_lo, vlen);

  VFLOAT t2_hi, t2_lo;
  RVV_LIBMVEC_IMPL_NAME (rvv_sleef_clean_hypot_ddsqu_pair, LMUL)
    (t_hi, t_lo, &t2_hi, &t2_lo, vlen);

  VFLOAT sum_hi, sum_lo;
  RVV_LIBMVEC_IMPL_NAME (rvv_sleef_clean_hypot_ddadd2_pair_scalar, LMUL)
    (t2_hi, t2_lo, VFMV_VF (fp_posOne, vlen), &sum_hi, &sum_lo, vlen);

  VFLOAT sqrt_hi, sqrt_lo;
  RVV_LIBMVEC_IMPL_NAME (rvv_sleef_clean_hypot_ddsqrt_pair, LMUL)
    (sum_hi, sum_lo, &sqrt_hi, &sqrt_lo, vlen);

  VFLOAT result_hi, result_lo;
  RVV_LIBMVEC_IMPL_NAME (rvv_sleef_clean_hypot_ddmul_pair_scalar, LMUL)
    (sqrt_hi, sqrt_lo, vmax_core, &result_hi, &result_lo, vlen);

  if (nb_scaled > 0)
    {
      VFLOAT result_scale = VFMV_VF (1.0, vlen);
      result_scale = VFMERGE_VFM
        (result_scale, RVV_HYPOT_CLEAN_SCALE_DOWN_D, is_tiny, vlen);
      result_scale = VFMERGE_VFM
        (result_scale, RVV_HYPOT_CLEAN_SCALE_UP_D, is_large, vlen);

      /* Do not multiply the low correction by the inverse scale separately:
         that product may underflow even when the final rescaled result is
         normal.  Conversely, when the high result overflows, adding a
         separately scaled low term could create Inf-Inf and leak invalid.
         One FMA preserves the public scalar exception behavior.  */
      VFLOAT result_hi_scaled
        = __riscv_vfmul (result_hi, result_scale, vlen);
      VUINT result_hi_class = __riscv_vfclass (result_hi_scaled, vlen);
      VBOOL result_hi_inf;
      IDENTIFY (result_hi_class, class_Inf, result_hi_inf, vlen);
      VFLOAT low_for_scale
        = VFMERGE_VFM (result_lo, fp_posZero, result_hi_inf, vlen);
      VFLOAT result_scaled
        = VFMADD_VV (low_for_scale, result_scale, result_hi_scaled, vlen);

      result_hi = __riscv_vmerge (result_hi, result_scaled, scale_mask,
                                  vlen);
      result_lo = VFMERGE_VFM (result_lo, fp_posZero, scale_mask, vlen);
    }
  vz = __riscv_vfadd (result_hi, result_lo, vlen);

  VBOOL is_nan_result = __riscv_vmfne (vz, vz, vlen);
  vz = __riscv_vfmerge (vz, fp_posInf, is_nan_result, vlen);

  VBOOL min_is_zero = __riscv_vmfeq (vmin, fp_posZero, vlen);
  vz = __riscv_vmerge (vz, vmax, min_is_zero, vlen);

  /* Merge special lanes back in.  */
  vz = VMERGE_VVM_F (vz, vz_special, special_args, vlen);

  RESTORE_FRM;
  return vz;
}
