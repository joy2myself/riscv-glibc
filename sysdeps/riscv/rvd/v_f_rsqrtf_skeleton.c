/* Single-precision vector rsqrtf function.
   rsqrtf(x) computes 1/sqrt(x)

   Based on aarch64 implementation.

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

#include <math.h>
#include <vecmath_config.h>
#include <v_math.h>

#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

/* Single-precision implementation of vector rsqrtf(x).  */
static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT vx = x;
  VFLOAT vy, vy_special;
  VUINT vclass;
  VBOOL zero_args, x_nan, x_pos_inf, x_neg_nonzero, unsafe_args,
      special_args;

  SET_ROUNDTONEAREST;

  vclass = __riscv_vfclass (vx, vlen);
  IDENTIFY_F (vclass, class_Zero, zero_args, vlen);
  IDENTIFY_F (vclass, class_NaN, x_nan, vlen);
  IDENTIFY_F (vclass, class_posInf, x_pos_inf, vlen);
  IDENTIFY_F (vclass, class_negative, x_neg_nonzero, vlen);
  unsafe_args = __riscv_vmor (x_nan, x_pos_inf, vlen);
  unsafe_args = __riscv_vmor (unsafe_args, x_neg_nonzero, vlen);
  special_args = __riscv_vmor (unsafe_args, zero_args, vlen);
  UINT nb_special = __riscv_vcpop (special_args, vlen);

  /* Prepare all exceptional results before sanitizing the core.  */
  if (nb_special > 0)
    {
      vy_special = VFMV_VF (fp_posZero_f, vlen);

      VBOOL id_mask;
      IDENTIFY_F (vclass, class_posZero, id_mask, vlen);
      vy_special = VFMERGE_VFM (vy_special, fp_posInf_f, id_mask, vlen);
      IDENTIFY_F (vclass, class_negZero, id_mask, vlen);
      vy_special = VFMERGE_VFM (vy_special, fp_negInf_f, id_mask, vlen);

      /* Quiet NaNs preserve the scalar qNaN/sNaN exception behavior.  */
      VFLOAT vy_nan = VFADD_VV_M (x_nan, vx, vx, vlen);
      vy_special = VMERGE_VVM_F (vy_special, vy_nan, x_nan, vlen);

      /* Negative nonzero inputs are domain errors.  */
      VFLOAT invalid_input
          = VFMERGE_VFM (vx, fp_sNaN_f, x_neg_nonzero, vlen);
      VFLOAT invalid_result
          = VFADD_VV_M (x_neg_nonzero, invalid_input, invalid_input, vlen);
      vy_special
          = VMERGE_VVM_F (vy_special, invalid_result, x_neg_nonzero, vlen);
    }

  /* Full precision path: rely on hardware sqrt and division.  */
  VFLOAT vx_core = VFMERGE_VFM (vx, fp_posOne_f, unsafe_args, vlen);
  vy = __riscv_vfrdiv (__riscv_vfsqrt (vx_core, vlen), 1.0f, vlen);

  if (nb_special > 0)
    vy = VMERGE_VVM_F (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
