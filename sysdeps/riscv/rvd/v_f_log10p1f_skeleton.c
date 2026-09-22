/* Single-precision vector log10p1f function.
   log10p1f(x) computes log10(1 + x) accurately for small x

   Based on log10f implementation with log10(e) scaling.

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

/* SLEEF has no native log10p1f entry.  Reuse the target-owned xlog10f
   translation for the finite composite, while keeping this public routine's
   clean EH envelope separate from the source-native comparator.  */
#define RVV_SLEEF_ORIG_LOG10F_HELPERS_ONLY
#define RVV_SLEEF_ORIG_LOG10F_CLEAN_KERNEL_ONLY
#include "v_f_log10f_original_skeleton.c"
#undef RVV_SLEEF_ORIG_LOG10F_CLEAN_KERNEL_ONLY
#undef RVV_SLEEF_ORIG_LOG10F_HELPERS_ONLY

#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

/* Match libm reference: log10f(1.0 + x) with float rounding.  */
static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT vxp1, vy_special;
  VBOOL special_args;

  SET_ROUNDTONEAREST;

  /* Compute 1.0f + x with float rounding, then apply log10f kernel.  */
  vxp1 = __riscv_vfadd (x, 1.0f, vlen);

  VUINT vclass = __riscv_vfclass (vxp1, vlen);
  IDENTIFY (vclass, 0x3BF, special_args, vlen);
  UINT nb_special_args = __riscv_vcpop (special_args, vlen);
  vy_special = VFMV_VF_F (0.0f, vlen);

  if (nb_special_args > 0)
    {
      VBOOL id_mask;
      IDENTIFY (vclass, class_negative, id_mask, vlen);
      vxp1 = VFMERGE_VFM (vxp1, fp_sNaN_f, id_mask, vlen);
      IDENTIFY (vclass, class_posZero, id_mask, vlen);
      vxp1 = VFMERGE_VFM (vxp1, fp_negZero_f, id_mask, vlen);
      IDENTIFY (vclass, 0x39F, special_args, vlen);
      vy_special = VFREC7_M (special_args, vxp1, vlen);
      vy_special = VFADD_VV_M (special_args, vy_special, vxp1, vlen);
      vxp1 = VFMERGE_VFM (vxp1, fp_posOne_f, special_args, vlen);
    }

  /* The finite algorithm is the same target-owned source-shaped xlog10f
     component used by _original.  Only the EH envelope above and the final
     merge below are clean-owned adaptations.  */
  VFLOAT vy = RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_log10f_clean_kernel_m, LMUL) (vxp1, vlen);

  vy = VMERGE_VVM_F (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
