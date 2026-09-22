/* Single-precision source-original SLEEF asin skeleton.

   This is a target-owned translation of SLEEF xasinf_u1.  Its finite-domain
   reduction, polynomial constants, double-float subtraction, and source
   operation order are kept separate from the clean implementation's
   exceptional-lane quarantine, endpoint repair, and fenv policy.  SLEEF
   source files are evidence only; this file has no SLEEF include or link
   dependency.

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

/* SLEEF dfsub_vf2_vf2_vf2, source common/df.h:172--180.  */
static inline void
RVV_SLEEF_ORIG_ASINF_MANGLE (rvv_sleef_orig_asinf_dfsub_pair_pair)
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

/* SLEEF dfsub_vf2_vf2_vf, source common/df.h:134--137.  */
static inline void
RVV_SLEEF_ORIG_ASINF_MANGLE (rvv_sleef_orig_asinf_dfsub_pair_scalar)
  (VFLOAT x_hi, VFLOAT x_lo, VFLOAT y, VFLOAT *r_hi, VFLOAT *r_lo,
   size_t vlen)
{
  VFLOAT s = VFSUB_VV (x_hi, y, vlen);
  VFLOAT t = VFSUB_VV (x_hi, s, vlen);
  t = VFSUB_VV (t, y, vlen);
  *r_hi = s;
  *r_lo = VFADD_VV (t, x_lo, vlen);
}

/* SLEEF dfrec_vf2_vf, source common/df.h:219--222.  */
static inline void
RVV_SLEEF_ORIG_ASINF_MANGLE (rvv_sleef_orig_asinf_dfrec_scalar)
  (VFLOAT d, VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT t = VFDIV_VV (VFMV_VF (1.0f, vlen), d, vlen);
  VFLOAT u = VFNMSUB_VV (d, t, VFMV_VF (1.0f, vlen), vlen);
  *r_hi = t;
  *r_lo = VFMUL_VV (t, u, vlen);
}

/* SLEEF dfsqrt_vf2_vf, source common/df.h:366--369.  */
static inline void
RVV_SLEEF_ORIG_ASINF_MANGLE (rvv_sleef_orig_asinf_dfsqrt_scalar)
  (VFLOAT d, VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT t = __riscv_vfsqrt (d, vlen);

  VFLOAT tsq_hi, tsq_lo;
  RVV_SLEEF_ORIG_DFMUL_VF_VF (t, t, tsq_hi, tsq_lo, vlen);

  VFLOAT sum_hi, sum_lo;
  RVV_SLEEF_ORIG_DFADD2_VF_VF2
    (d, tsq_hi, tsq_lo, sum_hi, sum_lo, vlen);

  VFLOAT rec_hi, rec_lo;
  RVV_SLEEF_ORIG_ASINF_MANGLE (rvv_sleef_orig_asinf_dfrec_scalar)
    (t, &rec_hi, &rec_lo, vlen);

  VFLOAT q_hi, q_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (sum_hi, sum_lo, rec_hi, rec_lo, q_hi, q_lo, vlen);

  *r_hi = VFMUL_VF (q_hi, 0.5f, vlen);
  *r_lo = VFMUL_VF (q_lo, 0.5f, vlen);
}

/* Keep the translated xasinf_u1 finite kernel reusable by SLEEF-component
   composites such as asinpif.  */
static VFLOAT
RVV_SLEEF_ORIG_ASINF_MANGLE (rvv_sleef_orig_asinf_kernel_m)
  (VFLOAT d, size_t vlen)
{
  VFLOAT abs_d = __riscv_vfabs (d, vlen);
  VBOOL small = __riscv_vmflt (abs_d, 0.5f, vlen);

  VFLOAT d2 = VFMUL_VV (d, d, vlen);
  VFLOAT z = VFMUL_VF
    (VFSUB_VV (VFMV_VF (1.0f, vlen), abs_d, vlen), 0.5f, vlen);
  VFLOAT x2 = VMERGE_VVM_F (z, d2, small, vlen);

  VFLOAT sqrt_hi, sqrt_lo;
  RVV_SLEEF_ORIG_ASINF_MANGLE (rvv_sleef_orig_asinf_dfsqrt_scalar)
    (x2, &sqrt_hi, &sqrt_lo, vlen);

  VFLOAT x_hi = VMERGE_VVM_F (sqrt_hi, abs_d, small, vlen);
  VFLOAT x_lo = VMERGE_VVM_F
    (sqrt_lo, VFMV_VF (0.0f, vlen), small, vlen);

  VBOOL abs_one = __riscv_vmfeq (abs_d, 1.0f, vlen);
  x_hi = VFMERGE_VFM (x_hi, fp_posZero, abs_one, vlen);
  x_lo = VFMERGE_VFM (x_lo, fp_posZero, abs_one, vlen);

  VFLOAT u = VFMV_VF (RVV_SLEEF_ORIG_ASINF_C4, vlen);
  u = VFMADD_VV (u, x2,
                 VFMV_VF (RVV_SLEEF_ORIG_ASINF_C3, vlen), vlen);
  u = VFMADD_VV (u, x2,
                 VFMV_VF (RVV_SLEEF_ORIG_ASINF_C2, vlen), vlen);
  u = VFMADD_VV (u, x2,
                 VFMV_VF (RVV_SLEEF_ORIG_ASINF_C1, vlen), vlen);
  u = VFMADD_VV (u, x2,
                 VFMV_VF (RVV_SLEEF_ORIG_ASINF_C0, vlen), vlen);
  u = VFMUL_VV (u, VFMUL_VV (x2, x_hi, vlen), vlen);

  VFLOAT pi4_hi = VFMV_VF (RVV_SLEEF_ORIG_ASINF_PI4_HI, vlen);
  VFLOAT pi4_lo = VFMV_VF (RVV_SLEEF_ORIG_ASINF_PI4_LO, vlen);
  VFLOAT y1_hi, y1_lo;
  RVV_SLEEF_ORIG_ASINF_MANGLE (rvv_sleef_orig_asinf_dfsub_pair_pair)
    (pi4_hi, pi4_lo, x_hi, x_lo, &y1_hi, &y1_lo, vlen);

  VFLOAT y_hi, y_lo;
  RVV_SLEEF_ORIG_ASINF_MANGLE (rvv_sleef_orig_asinf_dfsub_pair_scalar)
    (y1_hi, y1_lo, u, &y_hi, &y_lo, vlen);

  VFLOAT small_result = VFADD_VV (u, x_hi, vlen);
  VFLOAT large_result = VFMUL_VF
    (VFADD_VV (y_hi, y_lo, vlen), 2.0f, vlen);
  VFLOAT result = VMERGE_VVM_F (large_result, small_result, small, vlen);

  /* SLEEF vmulsign is sign-xor; result is positive for the finite source
     domain, but preserve the source operation explicitly.  */
  return __riscv_vfsgnjx (result, d, vlen);
}

#if !RVV_MATH_IS_ASINF_HELPERS_ONLY (V_MATH_FUNC)
static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT d, size_t vlen)
{
  return RVV_SLEEF_ORIG_ASINF_MANGLE
    (rvv_sleef_orig_asinf_kernel_m) (d, vlen);
}
#endif

#if !RVV_MATH_IS_ASINF_HELPERS_ONLY (V_MATH_FUNC)
#undef LMUL
#undef SIMDLEN
#endif
