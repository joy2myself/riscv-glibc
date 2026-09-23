/* Single-precision source-original SLEEF sinhf skeleton.

   This is a target-owned translation of SLEEF expk2f and xsinhf.  The
   source double-float reduction, polynomial, reconstruction, scaling, and
   late special-result selections remain separate from clean's EH path.
   Only the RVV representation, LMUL/VL plumbing, and glibc importer ABI are
   adapted; no SLEEF source file is included or linked.

   SLEEF requests -ffp-contract=off for this source.  The corresponding
   build flag is attached to the wrapper in the RVV mathvec Makefile.

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
  VFLOAT y = __riscv_vfabs (x, vlen);
  VFLOAT d_hi, d_lo;
  RVV_LIBMVEC_IMPL_NAME (sleef_f32_exp_pair, LMUL)
    (y, VFMV_VF (0.0f, vlen), &d_hi, &d_lo, vlen);

  VFLOAT e_hi, e_lo;
  RVV_SLEEF_ORIG_DFREC_VF2_VF2 (d_hi, d_lo, e_hi, e_lo, vlen);

  VFLOAT diff_hi, diff_lo;
  RVV_SLEEF_ORIG_DFSUB_VF2_VF2
    (d_hi, d_lo, e_hi, e_lo, diff_hi, diff_lo, vlen);
  y = VFMUL_VF (VFADD_VV (diff_hi, diff_lo, vlen),
                RVV_HALF_F, vlen);

  VBOOL result_special = __riscv_vmor
    (__riscv_vmfgt (__riscv_vfabs (x, vlen),
                    RVV_SLEEF_HYPERBOLIC_SOURCE_THRESHOLD, vlen),
     __riscv_vmfne (y, y, vlen), vlen);
  y = VFMERGE_VFM (y, fp_posInf_f, result_special, vlen);
  /* SLEEF vmulsign is sign-XOR (vfsgnjx), not sign-copy (vfsgnj).  The
     magnitude is normally nonnegative, but preserve the source primitive
     exactly for NaN/signed-zero/intermediate edge behavior.  */
  y = __riscv_vfsgnjx (y, x, vlen);

  /* Source uses vor_vm_vo32_vm(visnan(x), reinterpret(y)): a true mask
     contributes all-one bits, rather than selecting a clean qNaN value.  */
  VBOOL x_nan = __riscv_vmfne (x, x, vlen);
  VINT nan_bits = VMERGE_VXM_I (VMVI_VX (0, vlen), -1, x_nan, vlen);
  y = I_AS_F (__riscv_vor (F_AS_I (y), nan_bits, vlen));

  return y;
}

#undef LMUL
#undef SIMDLEN
