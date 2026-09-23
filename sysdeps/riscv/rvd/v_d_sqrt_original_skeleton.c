/* Double-precision source-original SLEEF sqrt skeleton.

   This is a target-owned translation of the stock RVV SLEEF xsqrt_u05
   kernel.  It is deliberately separate from the clean implementation's
   hardware-sqrt substitution and its exceptional-lane adaptation.  The
   source finite operation order, constants, and late selections are kept in
   this comparison implementation; SLEEF source files are evidence only.

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
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT d, size_t vlen)
{
  /* SLEEF xsqrt_u05 begins by selecting NaN for negative inputs.  In
     particular, negative zero is not selected and remains available for the
     source late zero selection below.  */
  VBOOL negative = __riscv_vmflt (d, fp_posZero, vlen);
  d = VFMERGE_VFM (d, fp_qNaN, negative, vlen);

  /* Keep the source tiny-input scaling and its corresponding output scale.  */
  VBOOL tiny = __riscv_vmflt (d, RVV_SQRT_TINY, vlen);
  VFLOAT d_scaled = VFMUL_VF (d, RVV_SQRT_SCALE, vlen);
  d = VMERGE_VVM_F (d, d_scaled, tiny, vlen);
  VFLOAT q = VFMERGE_VFM (VFMV_VF (fp_posOne, vlen),
                          RVV_SQRT_TINY_Q, tiny, vlen);

  /* Bit-level reciprocal-square-root seed from SLEEF's RVV helper.  */
  VUINT seed_bits = __riscv_vsub (
      VMVU_VX (RVV_SQRT_MAGIC, vlen),
      __riscv_vsrl (F_AS_U (d), 1, vlen), vlen);
  VFLOAT y = U_AS_F (seed_bits);

  VFLOAT x = VFMUL_VV (d, y, vlen);
  VFLOAT w = VFMUL_VV (VFMV_VF (0.5, vlen), y, vlen);
  y = VFNMSUB_VV (x, w, VFMV_VF (0.5, vlen), vlen);
  x = VFMADD_VV (x, y, x, vlen);
  w = VFMADD_VV (w, y, w, vlen);

  y = VFNMSUB_VV (x, w, VFMV_VF (0.5, vlen), vlen);
  x = VFMADD_VV (x, y, x, vlen);
  w = VFMADD_VV (w, y, w, vlen);

  y = VFNMSUB_VV (x, w, VFMV_VF (0.5, vlen), vlen);
  x = VFMADD_VV (x, y, x, vlen);
  w = VFMADD_VV (w, y, w, vlen);

  y = VFNMSUB_VV (x, w, VFMV_VF (1.5, vlen), vlen);
  w = VFADD_VV (w, w, vlen);
  w = VFMUL_VV (w, y, vlen);

  x = VFMUL_VV (w, d, vlen);
  y = VFMSUB_VV (w, d, x, vlen);
  VFLOAT z = VFNMSUB_VV (w, x, VFMV_VF (1.0, vlen), vlen);

  z = VFNMSUB_VV (w, y, z, vlen);
  w = VFMUL_VV (VFMV_VF (0.5, vlen), x, vlen);
  w = VFMADD_VV (w, z, y, vlen);
  w = VFADD_VV (w, x, vlen);
  w = VFMUL_VV (w, q, vlen);

  /* Preserve the source late selection order after the finite kernel.  */
  VBOOL zero = __riscv_vmfeq (d, fp_posZero, vlen);
  VBOOL inf = __riscv_vmfeq (d, fp_posInf, vlen);
  VBOOL zero_or_inf = __riscv_vmor (zero, inf, vlen);
  w = VMERGE_VVM_F (w, d, zero_or_inf, vlen);

  VBOOL negative_late = __riscv_vmflt (d, fp_posZero, vlen);
  w = VFMERGE_VFM (w, fp_qNaN, negative_late, vlen);

  return w;
}

#undef LMUL
#undef SIMDLEN
