/* Double-precision source-original SLEEF hypot skeleton.

   This is a target-owned translation of SLEEF xhypot_u05.  The source
   double-float operation order is kept independent from the clean entry's
   exceptional-lane quarantine, denominator guard, and final merge.  SLEEF
   source files are evidence only; this file has no SLEEF include or link
   dependency.

   The public argument order is adapted by v_math_importer_arg2.h: the glibc
   wrapper passes (y, x), while the source kernel below receives those values
   as its two operands and computes the symmetric hypot operation.

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

/* SLEEF ddadd2_vd2_vd2_vd, source common/dd.h:140--144.  */
static inline void
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_hypot_ddadd2_pair_scalar, LMUL)
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

/* SLEEF ddadd2_vd2_vd2_vd2, source common/dd.h:166--170.  */
static inline void
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_hypot_ddadd2_pair_pair, LMUL)
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

/* SLEEF ddmul_vd2_vd_vd, source common/dd.h:199--202.  */
static inline void
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_hypot_ddmul_scalar_scalar, LMUL)
  (VFLOAT x, VFLOAT y, VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT s = VFMUL_VV (x, y, vlen);
  *r_hi = s;
  *r_lo = VFMSUB_VV (x, y, s, vlen);
}

/* SLEEF ddmul_vd2_vd2_vd2, source common/dd.h:209--212.  */
static inline void
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_hypot_ddmul_pair_pair, LMUL)
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

/* SLEEF ddmul_vd2_vd2_vd, source common/dd.h:222--225.  */
static inline void
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_hypot_ddmul_pair_scalar, LMUL)
  (VFLOAT x_hi, VFLOAT x_lo, VFLOAT y, VFLOAT *r_hi, VFLOAT *r_lo,
   size_t vlen)
{
  VFLOAT s = VFMUL_VV (x_hi, y, vlen);
  VFLOAT t = VFMSUB_VV (x_hi, y, s, vlen);
  t = VFMADD_VV (x_lo, y, t, vlen);
  *r_hi = s;
  *r_lo = t;
}

/* SLEEF ddrec_vd2_vd, source common/dd.h:227--230.  */
static inline void
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_hypot_ddrec_scalar, LMUL)
  (VFLOAT d, VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT t = VFDIV_VV (VFMV_VF (1.0, vlen), d, vlen);
  VFLOAT u = VFNMSUB_VV (d, t, VFMV_VF (1.0, vlen), vlen);
  *r_hi = t;
  *r_lo = VFMUL_VV (t, u, vlen);
}

/* SLEEF dddiv_vd2_vd2_vd2, source common/dd.h:191--197.  The current
   target is the RVV FMA configuration, so retain that source branch's
   association and operation order.  */
static inline void
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_hypot_dddiv_pair_pair, LMUL)
  (VFLOAT n_hi, VFLOAT n_lo, VFLOAT d_hi, VFLOAT d_lo,
   VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT t = VFDIV_VV (VFMV_VF (1.0, vlen), d_hi, vlen);
  VFLOAT s = VFMUL_VV (n_hi, t, vlen);
  VFLOAT u = VFMSUB_VV (t, n_hi, s, vlen);
  VFLOAT v = VFNMSUB_VV
    (d_lo, t, VFNMSUB_VV (d_hi, t, VFMV_VF (1.0, vlen), vlen), vlen);
  VFLOAT correction = VFMADD_VV
    (n_lo, t, u, vlen);
  *r_hi = s;
  *r_lo = VFMADD_VV (s, v, correction, vlen);
}

/* SLEEF ddsqu_vd2_vd2, source common/dd.h:204--207.  */
static inline void
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_hypot_ddsqu_pair, LMUL)
  (VFLOAT x_hi, VFLOAT x_lo, VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT s = VFMUL_VV (x_hi, x_hi, vlen);
  VFLOAT t = VFMSUB_VV (x_hi, x_hi, s, vlen);
  t = VFMADD_VV (VFADD_VV (x_hi, x_hi, vlen), x_lo, t, vlen);
  *r_hi = s;
  *r_lo = t;
}

/* SLEEF ddsqrt_vd2_vd2, source common/dd.h:312--315.  */
static inline void
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_hypot_ddsqrt_pair, LMUL)
  (VFLOAT d_hi, VFLOAT d_lo, VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT t = __riscv_vfsqrt (VFADD_VV (d_hi, d_lo, vlen), vlen);

  VFLOAT tsq_hi, tsq_lo;
  RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_hypot_ddmul_scalar_scalar, LMUL)
    (t, t, &tsq_hi, &tsq_lo, vlen);

  VFLOAT sum_hi, sum_lo;
  RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_hypot_ddadd2_pair_pair, LMUL)
    (d_hi, d_lo, tsq_hi, tsq_lo, &sum_hi, &sum_lo, vlen);

  VFLOAT rec_hi, rec_lo;
  RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_hypot_ddrec_scalar, LMUL)
    (t, &rec_hi, &rec_lo, vlen);

  VFLOAT q_hi, q_lo;
  RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_hypot_ddmul_pair_pair, LMUL)
    (sum_hi, sum_lo, rec_hi, rec_lo, &q_hi, &q_lo, vlen);

  *r_hi = VFMUL_VF (q_hi, 0.5, vlen);
  *r_lo = VFMUL_VF (q_lo, 0.5, vlen);
}

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT vx, VFLOAT vy,
                                           size_t vlen)
{
  /* The source RVV path uses the target's normal RNE/FMA arithmetic; no
     clean-side special-lane quarantine is inserted into this comparison.  */
  VFLOAT x = __riscv_vfabs (vx, vlen);
  VFLOAT y = __riscv_vfabs (vy, vlen);
  VFLOAT vmin = __riscv_vfmin (x, y, vlen);
  VFLOAT vmax = __riscv_vfmax (x, y, vlen);

  VBOOL scale = __riscv_vmflt (vmax, RVV_HYPOT_SOURCE_MIN_D, vlen);
  VFLOAT n_scaled = VFMUL_VF
    (vmin, RVV_HYPOT_SOURCE_SCALE_D, vlen);
  VFLOAT d_scaled = VFMUL_VF
    (vmax, RVV_HYPOT_SOURCE_SCALE_D, vlen);
  VFLOAT n = VMERGE_VVM_F (vmin, n_scaled, scale, vlen);
  VFLOAT d = VMERGE_VVM_F (vmax, d_scaled, scale, vlen);

  VFLOAT zero = VFMV_VF (0.0, vlen);
  VFLOAT t_hi, t_lo;
  RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_hypot_dddiv_pair_pair, LMUL)
    (n, zero, d, zero, &t_hi, &t_lo, vlen);

  VFLOAT t2_hi, t2_lo;
  RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_hypot_ddsqu_pair, LMUL)
    (t_hi, t_lo, &t2_hi, &t2_lo, vlen);

  VFLOAT sum_hi, sum_lo;
  RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_hypot_ddadd2_pair_scalar, LMUL)
    (t2_hi, t2_lo, VFMV_VF (1.0, vlen), &sum_hi, &sum_lo, vlen);

  VFLOAT root_hi, root_lo;
  RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_hypot_ddsqrt_pair, LMUL)
    (sum_hi, sum_lo, &root_hi, &root_lo, vlen);

  VFLOAT result_hi, result_lo;
  RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_hypot_ddmul_pair_scalar, LMUL)
    (root_hi, root_lo, vmax, &result_hi, &result_lo, vlen);
  VFLOAT result = VFADD_VV (result_hi, result_lo, vlen);

  /* Preserve xhypot_u05's late selection order.  Exceptional behavior here
     is intentionally source-native; clean's EH path belongs to the clean
     implementation and is not copied into _original.  */
  VBOOL result_nan = __riscv_vmfne (result, result, vlen);
  result = VFMERGE_VFM (result, fp_posInf, result_nan, vlen);

  VBOOL min_zero = __riscv_vmfeq (vmin, fp_posZero, vlen);
  result = VMERGE_VVM_F (result, vmax, min_zero, vlen);

  VBOOL any_nan = __riscv_vmor
    (__riscv_vmfne (x, x, vlen), __riscv_vmfne (y, y, vlen), vlen);
  result = VFMERGE_VFM (result, fp_qNaN, any_nan, vlen);

  VBOOL any_inf = __riscv_vmor
    (__riscv_vmfeq (x, fp_posInf, vlen),
     __riscv_vmfeq (y, fp_posInf, vlen), vlen);
  result = VFMERGE_VFM (result, fp_posInf, any_inf, vlen);

  return result;
}

#undef LMUL
#undef SIMDLEN
