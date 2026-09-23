/* Single-precision source-original SLEEF sqrtf skeleton.

   This is a target-owned translation of SLEEF xsqrtf_u05.  It keeps the
   source scaling, reciprocal-square-root bit seed, FMA Newton refinements,
   final reconstruction, and source special selections.  The active clean
   implementation uses RVV hardware vfsqrt; that is deliberately not reused
   here because this file is the source-algorithm comparison implementation.

   SLEEF requests -ffp-contract=off for its source.  The wrapper is compiled
   with the same option by the RVV mathvec Makefile.

   Copyright (C) 2026 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public License as
   published by the Free Software Foundation; either version 2.1 of the
   License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with the GNU C Library; if not, see
   <https://www.gnu.org/licenses/>.  */

#include <vecmath_config.h>
#include <v_math.h>

#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

/* Private scaling and reciprocal-square-root seed from the source graph.  */
#define RVV_SQRTF_ORIGINAL_SCALE 1.8889465931478580e+22f
#define RVV_SQRTF_ORIGINAL_SMALL_SCALE 7.2759576141834260e-12f
#define RVV_SQRTF_ORIGINAL_SMALL_LIMIT 5.2939559203393770e-23f
#define RVV_SQRTF_ORIGINAL_MAGIC 0x5f3759dfu

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT d, size_t vlen)
{
  VFLOAT q, w, x, y, z;

  /* Source xsqrtf_u05 begins with a numerical selection of a quiet NaN for
     negative inputs; it does not quarantine those lanes from later vector
     operations.  */
  VBOOL negative = __riscv_vmflt (d, fp_posZero_f, vlen);
  d = VMERGE_VVM_F (d, VFMV_VF (fp_qNaN_f, vlen), negative, vlen);

  VBOOL small = __riscv_vmflt
    (d, RVV_SQRTF_ORIGINAL_SMALL_LIMIT, vlen);
  VFLOAT scaled = VFMUL_VF
    (d, RVV_SQRTF_ORIGINAL_SCALE, vlen);
  d = VMERGE_VVM_F (d, scaled, small, vlen);
  q = VFMERGE_VFM
    (VFMV_VF (1.0f, vlen), RVV_SQRTF_ORIGINAL_SMALL_SCALE,
     small, vlen);

  /* y = reinterpret(0x5f3759df - (reinterpret(d) >> 1)).  */
  VUINT d_bits = F_AS_U (d);
  VUINT seed_bits = __riscv_vsub
    (VMVU_VX (RVV_SQRTF_ORIGINAL_MAGIC, vlen),
     __riscv_vsrl (d_bits, 1, vlen), vlen);
  y = U_AS_F (seed_bits);

  /* First two reciprocal-square-root Newton refinements.  The vfmanp/vfmapn
     source operations are mapped to target-owned fused operations with the
     same operand order.  */
  x = VFMUL_VV (d, y, vlen);
  w = VFMUL_VF (y, 0.5f, vlen);
  y = VFNMSUB_VV (x, w, VFMV_VF (0.5f, vlen), vlen);
  x = VFMADD_VV (x, y, x, vlen);
  w = VFMADD_VV (w, y, w, vlen);
  y = VFNMSUB_VV (x, w, VFMV_VF (0.5f, vlen), vlen);
  x = VFMADD_VV (x, y, x, vlen);
  w = VFMADD_VV (w, y, w, vlen);

  /* Source final reconstruction.  */
  y = VFNMSUB_VV (x, w, VFMV_VF (1.5f, vlen), vlen);
  w = VFADD_VV (w, w, vlen);
  w = VFMUL_VV (w, y, vlen);
  x = VFMUL_VV (w, d, vlen);
  y = VFMSUB_VV (w, d, x, vlen);
  z = VFNMSUB_VV (w, x, VFMV_VF (1.0f, vlen), vlen);

  z = VFNMSUB_VV (w, y, z, vlen);
  w = VFMUL_VF (x, 0.5f, vlen);
  w = VFMADD_VV (w, z, y, vlen);
  w = VFADD_VV (w, x, vlen);
  w = VFMUL_VV (w, q, vlen);

  /* Keep the source's post-core selection order.  d is intentionally the
     source-updated value here (small inputs have already been rescaled and
     negative inputs have already become NaN).  */
  VBOOL zero_or_pos_inf = __riscv_vmor
    (__riscv_vmfeq (d, fp_posZero_f, vlen),
     __riscv_vmfeq (d, fp_posInf_f, vlen), vlen);
  w = VMERGE_VVM_F (w, d, zero_or_pos_inf, vlen);
  w = VFMERGE_VFM (w, fp_qNaN_f,
                   __riscv_vmflt (d, fp_posZero_f, vlen), vlen);

  return w;
}

#undef RVV_SQRTF_ORIGINAL_SCALE
#undef RVV_SQRTF_ORIGINAL_SMALL_SCALE
#undef RVV_SQRTF_ORIGINAL_SMALL_LIMIT
#undef RVV_SQRTF_ORIGINAL_MAGIC

#undef LMUL
#undef SIMDLEN
