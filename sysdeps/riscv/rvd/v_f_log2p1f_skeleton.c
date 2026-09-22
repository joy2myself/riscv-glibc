/* Single-precision vector log2p1f function.
   log2p1f(x) computes log2(1 + x) accurately for small x.

   The pinned SLEEF revision has no native log2p1f entry.  The finite
   comparator contract is therefore the explicitly rounded transform
   dp1 = round_RN(x + 1), followed by the target-owned xlog2f kernel.

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

/* Reuse the target-owned translation of SLEEF xlog2f for the finite core.
   The clean EH envelope below remains separate from this source-shaped
   kernel.  */
#include "v_f_log2f_original_skeleton.c"

/* The source kernel selects its clean masked subnormal-scale adapter from the
   importing function name; no routine-local compile-control macro is needed.  */

/* Match libm reference: log2f(1.0 + x) with float rounding.  */
static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT vx, vy, vy_special;
  VBOOL special_args;

  SET_ROUNDTONEAREST;

  /* The composite source contract deliberately rounds the transformed
     argument before classification and finite evaluation.  */
  vx = x;
  vx = __riscv_vfadd (vx, 1.0f, vlen);

  VUINT vclass = __riscv_vfclass (vx, vlen);
  IDENTIFY (vclass, 0x3BF, special_args, vlen);
  UINT nb_special_args = __riscv_vcpop (special_args, vlen);
  vy_special = VFMV_VF_F (0.0f, vlen);

  if (nb_special_args > 0)
    {
      VBOOL id_mask;
      IDENTIFY (vclass, class_negative, id_mask, vlen);
      vx = VFMERGE_VFM (vx, fp_sNaN_f, id_mask, vlen);
      IDENTIFY (vclass, class_posZero, id_mask, vlen);
      vx = VFMERGE_VFM (vx, fp_negZero_f, id_mask, vlen);
      IDENTIFY (vclass, 0x39F, special_args, vlen);
      vy_special = VFREC7_M (special_args, vx, vlen);
      vy_special = VFADD_VV_M (special_args, vy_special, vx, vlen);
      vx = VFMERGE_VFM (vx, fp_posOne_f, special_args, vlen);
    }

  /* Keep the source xlog2f finite graph, but quarantine a large finite
     transformed input before its d*(4/3) reduction.  Scaling by 2^-4
     lowers the extracted exponent by four; adding four after the source
     kernel is exact for this result range and avoids a spurious overflow.  */
  VBOOL adjust_large = __riscv_vmfge (vx, RVV_LOG2F_CLEAN_BOUND_F, vlen);
  VFLOAT vx_scaled = VFMUL_VF_M
    (adjust_large, vx, RVV_LOG2F_CLEAN_SCALE_DOWN_F, vlen);
  VFLOAT vx_core = VMERGE_VVM_F (vx, vx_scaled, adjust_large, vlen);
  vy = RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_log2f_kernel_m, LMUL) (vx_core, vlen);
  VFLOAT vy_large = VFADD_VF_M
    (adjust_large, vy, RVV_LOG2F_CLEAN_LARGE_EXP_ADJUST_F, vlen);
  vy = VMERGE_VVM_F (vy, vy_large, adjust_large, vlen);

  vy = VMERGE_VVM_F (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
