/* Single-precision vector sqrtf function.
   sqrtf(x) computes the square root of x

   Uses RVV hardware sqrt with a clean exceptional-lane envelope.
   The separate source-original entry retains SLEEF Newton refinement.

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
  VFLOAT vx = x;
  VFLOAT vy, vy_special;
  VUINT vclass = __riscv_vfclass (vx, vlen);
  VBOOL x_nan, x_inf, x_pos_inf, x_neg_nonzero, unsafe_args;

  IDENTIFY_F (vclass, class_NaN, x_nan, vlen);
  IDENTIFY_F (vclass, class_Inf, x_inf, vlen);
  IDENTIFY_F (vclass, class_posInf, x_pos_inf, vlen);
  IDENTIFY_F (vclass, class_negative, x_neg_nonzero, vlen);

  unsafe_args = __riscv_vmor (x_nan, x_inf, vlen);
  unsafe_args = __riscv_vmor (unsafe_args, x_neg_nonzero, vlen);
  UINT nb_unsafe = __riscv_vcpop (unsafe_args, vlen);

  /* Prepare exceptional results before sanitizing the hardware core.  */
  if (nb_unsafe > 0)
    {
      vy_special = VFMV_VF (fp_posZero_f, vlen);

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

      /* +Inf is the only infinity with a non-error sqrt result.  */
      vy_special = VFMERGE_VFM (vy_special, fp_posInf_f, x_pos_inf, vlen);
    }

  /* Keep finite values and signed zero on the hardware sqrt path.  */
  VFLOAT vx_core = VFMERGE_VFM (vx, fp_posOne_f, unsafe_args, vlen);
  vy = __riscv_vfsqrt (vx_core, vlen);

  if (nb_unsafe > 0)
    vy = VMERGE_VVM_F (vy, vy_special, unsafe_args, vlen);

  return vy;
}

#undef LMUL
#undef SIMDLEN
