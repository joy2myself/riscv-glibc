/* Single-precision source-original SLEEF acosf skeleton.

   This is a target-owned translation of SLEEF xacosf_u1.  The source
   finite-domain predicate, double-float operations, polynomial constants,
   and operation order are kept independent from the clean implementation.
   The SLEEF source files are evidence only; this translation has no SLEEF
   include or link dependency.

   In particular, this comparison entry intentionally does not import the
   clean-side EH quarantine or endpoint repair.  xacosf_u1's native vector
   behavior is the comparison contract; the clean entry is responsible for
   the glibc/VECLIBM exception path.

   Copyright (C) 2026 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public License as
   published by the Free Software Foundation; either version 2.1 of the
   License, or (at your option) any later version.

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

/* SLEEF dfadd_vf2_vf_vf.  */
static inline void
RVV_SLEEF_ORIG_ACOSF_MANGLE (rvv_sleef_orig_acosf_dfadd_vf_vf_m)
  (VFLOAT x, VFLOAT y, VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT s = VFADD_VV (x, y, vlen);
  *r_hi = s;
  *r_lo = VFADD_VV (VFSUB_VV (x, s, vlen), y, vlen);
}

/* SLEEF dfadd2_vf2_vf_vf2.  */
static inline void
RVV_SLEEF_ORIG_ACOSF_MANGLE (rvv_sleef_orig_acosf_dfadd2_vf_vf2_m)
  (VFLOAT x, VFLOAT y_hi, VFLOAT y_lo, VFLOAT *r_hi, VFLOAT *r_lo,
   size_t vlen)
{
  VFLOAT s = VFADD_VV (x, y_hi, vlen);
  VFLOAT v = VFSUB_VV (s, x, vlen);
  VFLOAT t = VFADD_VV
    (VFSUB_VV (x, VFSUB_VV (s, v, vlen), vlen),
     VFSUB_VV (y_hi, v, vlen), vlen);
  *r_hi = s;
  *r_lo = VFADD_VV (t, y_lo, vlen);
}

/* SLEEF dfsub_vf2_vf2_vf2.  The source contract supplies |x| >= |y|.  */
static inline void
RVV_SLEEF_ORIG_ACOSF_MANGLE (rvv_sleef_orig_acosf_dfsub_vf2_vf2_m)
  (VFLOAT x_hi, VFLOAT x_lo, VFLOAT y_hi, VFLOAT y_lo,
   VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT s = VFSUB_VV (x_hi, y_hi, vlen);
  VFLOAT t = VFSUB_VV (x_hi, s, vlen);
  t = VFSUB_VV (t, y_hi, vlen);
  t = VFADD_VV (t, x_lo, vlen);
  *r_hi = s;
  *r_lo = VFSUB_VV (t, y_lo, vlen);
}

/* SLEEF dfmul_vf2_vf_vf in the pinned ENABLE_FMA_SP branch.  */
static inline void
RVV_SLEEF_ORIG_ACOSF_MANGLE (rvv_sleef_orig_acosf_dfmul_vf_vf_m)
  (VFLOAT x, VFLOAT y, VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  RVV_SLEEF_ORIG_DFMUL_VF_VF (x, y, *r_hi, *r_lo, vlen);
}

/* SLEEF dfrec_vf2_vf.  */
static inline void
RVV_SLEEF_ORIG_ACOSF_MANGLE (rvv_sleef_orig_acosf_dfrec_vf_m)
  (VFLOAT d, VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT t = VFDIV_VV (VFMV_VF (1.0f, vlen), d, vlen);
  VFLOAT u = VFNMSUB_VV (d, t, VFMV_VF (1.0f, vlen), vlen);
  *r_hi = t;
  *r_lo = VFMUL_VV (t, u, vlen);
}

/* Forward declaration because dfsqrt uses the pair-by-pair source helper.  */
static inline void
RVV_SLEEF_ORIG_ACOSF_MANGLE (rvv_sleef_orig_acosf_dfmul_vf2_vf2_m)
  (VFLOAT x_hi, VFLOAT x_lo, VFLOAT y_hi, VFLOAT y_lo,
   VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen);

/* SLEEF dfsqrt_vf2_vf, translated without importing SLEEF df.h.  */
static inline void
RVV_SLEEF_ORIG_ACOSF_MANGLE (rvv_sleef_orig_acosf_dfsqrt_m)
  (VFLOAT d, VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT t = __riscv_vfsqrt (d, vlen);

  VFLOAT tsq_hi, tsq_lo;
  RVV_SLEEF_ORIG_ACOSF_MANGLE (rvv_sleef_orig_acosf_dfmul_vf_vf_m)
    (t, t, &tsq_hi, &tsq_lo, vlen);

  VFLOAT sum_hi, sum_lo;
  RVV_SLEEF_ORIG_ACOSF_MANGLE (rvv_sleef_orig_acosf_dfadd2_vf_vf2_m)
    (d, tsq_hi, tsq_lo, &sum_hi, &sum_lo, vlen);

  VFLOAT rec_hi, rec_lo;
  RVV_SLEEF_ORIG_ACOSF_MANGLE (rvv_sleef_orig_acosf_dfrec_vf_m)
    (t, &rec_hi, &rec_lo, vlen);

  VFLOAT q_hi, q_lo;
  RVV_SLEEF_ORIG_ACOSF_MANGLE (rvv_sleef_orig_acosf_dfmul_vf2_vf2_m)
    (sum_hi, sum_lo, rec_hi, rec_lo, &q_hi, &q_lo, vlen);

  RVV_SLEEF_ORIG_DFSCALE_VF2_VF
    (q_hi, q_lo, VFMV_VF (0.5f, vlen), *r_hi, *r_lo, vlen);
}

/* SLEEF dfmul_vf2_vf2_vf2 in the pinned ENABLE_FMA_SP branch.  */
static inline void
RVV_SLEEF_ORIG_ACOSF_MANGLE (rvv_sleef_orig_acosf_dfmul_vf2_vf2_m)
  (VFLOAT x_hi, VFLOAT x_lo, VFLOAT y_hi, VFLOAT y_lo,
   VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (x_hi, x_lo, y_hi, y_lo, *r_hi, *r_lo, vlen);
}

/* Keep the translated xacosf_u1 finite kernel reusable by SLEEF-component
   composites such as acospif.  The helper-only include path below imports
   this function without creating a second public implementation.  */
static VFLOAT
RVV_SLEEF_ORIG_ACOSF_MANGLE (rvv_sleef_orig_acosf_kernel_m)
  (VFLOAT d, size_t vlen)
{
  /* xacosf_u1 does not establish an internal rounding-mode envelope.  Keep
     the source operation in the caller's mode; no clean EH state is added.  */
  VFLOAT abs_d = __riscv_vfabs (d, vlen);
  VBOOL small = __riscv_vmflt (abs_d, 0.5f, vlen);

  VFLOAT d2 = VFMUL_VV (d, d, vlen);
  VFLOAT z = VFMUL_VV
    (VFSUB_VV (VFMV_VF (1.0f, vlen), abs_d, vlen),
     VFMV_VF (0.5f, vlen), vlen);
  VFLOAT x2 = __riscv_vmerge (z, d2, small, vlen);

  VFLOAT sqrt_hi, sqrt_lo;
  RVV_SLEEF_ORIG_ACOSF_MANGLE (rvv_sleef_orig_acosf_dfsqrt_m)
    (x2, &sqrt_hi, &sqrt_lo, vlen);

  VFLOAT x_hi = __riscv_vmerge (sqrt_hi, abs_d, small, vlen);
  VFLOAT x_lo = __riscv_vmerge (sqrt_lo, VFMV_VF (0.0f, vlen), small, vlen);

  VBOOL abs_one = __riscv_vmfeq (abs_d, 1.0f, vlen);
  x_hi = __riscv_vfmerge (x_hi, fp_posZero, abs_one, vlen);
  x_lo = __riscv_vfmerge (x_lo, fp_posZero, abs_one, vlen);

  VFLOAT u = VFMV_VF (RVV_SLEEF_ORIG_ACOSF_C4, vlen);
  u = VFMADD_VV (u, x2,
                 VFMV_VF (RVV_SLEEF_ORIG_ACOSF_C3, vlen), vlen);
  u = VFMADD_VV (u, x2,
                 VFMV_VF (RVV_SLEEF_ORIG_ACOSF_C2, vlen), vlen);
  u = VFMADD_VV (u, x2,
                 VFMV_VF (RVV_SLEEF_ORIG_ACOSF_C1, vlen), vlen);
  u = VFMADD_VV (u, x2,
                 VFMV_VF (RVV_SLEEF_ORIG_ACOSF_C0, vlen), vlen);
  u = VFMUL_VV (u, VFMUL_VV (x2, x_hi, vlen), vlen);

  VFLOAT signed_x_hi = __riscv_vfsgnjx (x_hi, d, vlen);
  VFLOAT signed_u = __riscv_vfsgnjx (u, d, vlen);
  VFLOAT small_term_hi, small_term_lo;
  RVV_SLEEF_ORIG_ACOSF_MANGLE (rvv_sleef_orig_acosf_dfadd_vf_vf_m)
    (signed_x_hi, signed_u, &small_term_hi, &small_term_lo, vlen);

  VFLOAT y_small_hi, y_small_lo;
  RVV_SLEEF_ORIG_ACOSF_MANGLE (rvv_sleef_orig_acosf_dfsub_vf2_vf2_m)
    (VFMV_VF (RVV_SLEEF_ORIG_ACOSF_PI2_HI, vlen),
     VFMV_VF (RVV_SLEEF_ORIG_ACOSF_PI2_LO, vlen),
     small_term_hi, small_term_lo, &y_small_hi, &y_small_lo, vlen);

  VFLOAT x_plus_u_hi, x_plus_u_lo;
  RVV_SLEEF_ORIG_DFADD_VF2_VF2_VF
    (x_hi, x_lo, u, x_plus_u_hi, x_plus_u_lo, vlen);

  VFLOAT y_large_hi, y_large_lo;
  RVV_SLEEF_ORIG_DFSCALE_VF2_VF
    (x_plus_u_hi, x_plus_u_lo, VFMV_VF (2.0f, vlen),
     y_large_hi, y_large_lo, vlen);

  VBOOL neg_large = __riscv_vmandn
    (__riscv_vmflt (d, fp_posZero, vlen), small, vlen);
  VFLOAT y_neg_hi, y_neg_lo;
  RVV_SLEEF_ORIG_ACOSF_MANGLE (rvv_sleef_orig_acosf_dfsub_vf2_vf2_m)
    (VFMV_VF (RVV_SLEEF_ORIG_ACOSF_PI_HI, vlen),
     VFMV_VF (RVV_SLEEF_ORIG_ACOSF_PI_LO, vlen),
     y_large_hi, y_large_lo, &y_neg_hi, &y_neg_lo, vlen);
  y_large_hi = __riscv_vmerge (y_large_hi, y_neg_hi, neg_large, vlen);
  y_large_lo = __riscv_vmerge (y_large_lo, y_neg_lo, neg_large, vlen);

  VFLOAT y_hi = __riscv_vmerge (y_large_hi, y_small_hi, small, vlen);
  VFLOAT y_lo = __riscv_vmerge (y_large_lo, y_small_lo, small, vlen);
  return VFADD_VV (y_hi, y_lo, vlen);
}

#if !RVV_MATH_IS_ACOSF_HELPERS_ONLY (V_MATH_FUNC)
static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT d, size_t vlen)
{
  return RVV_SLEEF_ORIG_ACOSF_MANGLE
    (rvv_sleef_orig_acosf_kernel_m) (d, vlen);
}
#endif

#if !RVV_MATH_IS_ACOSF_HELPERS_ONLY (V_MATH_FUNC)
#undef LMUL
#undef SIMDLEN
#endif
