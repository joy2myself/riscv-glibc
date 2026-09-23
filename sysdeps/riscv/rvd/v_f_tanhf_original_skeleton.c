/* Single-precision source-original SLEEF tanhf skeleton.

   This is a target-owned translation of SLEEF expk2f and xtanhf.  The
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

  VFLOAT neg_e_hi, neg_e_lo;
  RVV_SLEEF_ORIG_DFNEG_VF2 (e_hi, e_lo,
                            neg_e_hi, neg_e_lo, vlen);

  VFLOAT num_hi, num_lo;
  RVV_SLEEF_ORIG_DFADD_VF2_VF2
    (d_hi, d_lo, neg_e_hi, neg_e_lo, num_hi, num_lo, vlen);
  VFLOAT den_hi, den_lo;
  RVV_SLEEF_ORIG_DFADD_VF2_VF2
    (d_hi, d_lo, e_hi, e_lo, den_hi, den_lo, vlen);

  VFLOAT result_hi, result_lo;
  RVV_SLEEF_ORIG_DFDIV_VF2_VF2
    (num_hi, num_lo, den_hi, den_lo, result_hi, result_lo, vlen);
  y = VFADD_VV (result_hi, result_lo, vlen);

  VBOOL result_special = __riscv_vmor
    (__riscv_vmfgt (__riscv_vfabs (x, vlen),
                    RVV_TANHF_SATURATION_THRESHOLD_F, vlen),
     __riscv_vmfne (y, y, vlen), vlen);
  y = VFMERGE_VFM (y, fp_posOne_f, result_special, vlen);
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
