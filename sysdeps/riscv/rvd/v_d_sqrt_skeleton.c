/* Double-precision vector sqrt function.
   sqrt(x) computes the square root of x

   Based on the stock RVV SLEEF xsqrt_u05 finite kernel.  The public clean
   entry keeps the target's exception-handling quarantine around that kernel.

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
  VBOOL x_nan, x_inf, x_pos_inf, x_neg_nonzero, x_zero, unsafe_args;

  IDENTIFY (vclass, class_NaN, x_nan, vlen);
  IDENTIFY (vclass, class_Inf, x_inf, vlen);
  IDENTIFY (vclass, class_posInf, x_pos_inf, vlen);
  IDENTIFY (vclass, class_negative, x_neg_nonzero, vlen);
  IDENTIFY (vclass, class_Zero, x_zero, vlen);

  unsafe_args = __riscv_vmor (x_nan, x_inf, vlen);
  unsafe_args = __riscv_vmor (unsafe_args, x_neg_nonzero, vlen);
  unsafe_args = __riscv_vmor (unsafe_args, x_zero, vlen);
  UINT nb_unsafe = __riscv_vcpop (unsafe_args, vlen);

  /* Prepare exceptional results before sanitizing the hardware core.  */
  if (nb_unsafe > 0)
    {
      vy_special = VFMV_VF (fp_posZero, vlen);

      /* Quiet NaNs preserve the scalar qNaN/sNaN exception behavior.  */
      VFLOAT vy_nan = VFADD_VV_M (x_nan, vx, vx, vlen);
      vy_special = VMERGE_VVM_F (vy_special, vy_nan, x_nan, vlen);

      /* Negative nonzero inputs are domain errors.  */
      VFLOAT invalid_input
          = VFMERGE_VFM (vx, fp_sNaN, x_neg_nonzero, vlen);
      VFLOAT invalid_result
          = VFADD_VV_M (x_neg_nonzero, invalid_input, invalid_input, vlen);
      vy_special
          = VMERGE_VVM_F (vy_special, invalid_result, x_neg_nonzero, vlen);

      /* +Inf is the only infinity with a non-error sqrt result.  */
      vy_special = VFMERGE_VFM (vy_special, fp_posInf, x_pos_inf, vlen);

      /* Preserve both signs of zero without exposing them to the source
         bit-level reciprocal-square-root seed.  */
      vy_special = VMERGE_VVM_F (vy_special, vx, x_zero, vlen);
    }

  /* Run the source finite kernel only on positive finite values.  Every
     exceptional or zero lane is replaced before the bit-level seed and the
     fused refinement sequence are evaluated.  */
  VFLOAT vx_core = VFMERGE_VFM (vx, fp_posOne, unsafe_args, vlen);
  VBOOL core_args = __riscv_vmnot (unsafe_args, vlen);

  VFLOAT q, w, y, z;

  VBOOL tiny = __riscv_vmflt
    (vx_core, RVV_SQRT_TINY, vlen);
  VFLOAT scaled;
  scaled = VFMUL_VF_M (core_args, vx_core, RVV_SQRT_SCALE,
                       vlen);
  vx_core = VMERGE_VVM_F (vx_core, scaled, tiny, vlen);
  q = VFMERGE_VFM
    (VFMV_VF (fp_posOne, vlen), RVV_SQRT_TINY_Q, tiny, vlen);

  VUINT seed_bits = __riscv_vsub
    (VMVU_VX (RVV_SQRT_MAGIC, vlen),
     __riscv_vsrl (F_AS_U (vx_core), 1, vlen), vlen);
  y = U_AS_F (seed_bits);

  VFLOAT core_half = VFMV_VF (0.5, vlen);
  VFLOAT core_one = VFMV_VF (1.0, vlen);
  VFLOAT core_three_halves = VFMV_VF (1.5, vlen);

  VFLOAT core_x = VFMUL_VV_M (core_args, vx_core, y, vlen);
  w = VFMUL_VV_M (core_args, core_half, y, vlen);
  y = VFNMSUB_VV_M (core_args, core_x, w, core_half, vlen);
  core_x = VFMADD_VV_M (core_args, core_x, y, core_x, vlen);
  w = VFMADD_VV_M (core_args, w, y, w, vlen);

  y = VFNMSUB_VV_M (core_args, core_x, w, core_half, vlen);
  core_x = VFMADD_VV_M (core_args, core_x, y, core_x, vlen);
  w = VFMADD_VV_M (core_args, w, y, w, vlen);

  y = VFNMSUB_VV_M (core_args, core_x, w, core_half, vlen);
  core_x = VFMADD_VV_M (core_args, core_x, y, core_x, vlen);
  w = VFMADD_VV_M (core_args, w, y, w, vlen);

  y = VFNMSUB_VV_M (core_args, core_x, w, core_three_halves, vlen);
  w = VFADD_VV_M (core_args, w, w, vlen);
  w = VFMUL_VV_M (core_args, w, y, vlen);

  core_x = VFMUL_VV_M (core_args, w, vx_core, vlen);
  y = VFMSUB_VV_M (core_args, w, vx_core, core_x, vlen);
  z = VFNMSUB_VV_M (core_args, w, core_x, core_one, vlen);

  z = VFNMSUB_VV_M (core_args, w, y, z, vlen);
  w = VFMUL_VV_M (core_args, core_half, core_x, vlen);
  w = VFMADD_VV_M (core_args, w, z, y, vlen);
  w = VFADD_VV_M (core_args, w, core_x, vlen);
  vy = VFMUL_VV_M (core_args, w, q, vlen);

  if (nb_unsafe > 0)
    vy = VMERGE_VVM_F (vy, vy_special, unsafe_args, vlen);

  return vy;
}

#undef LMUL
#undef SIMDLEN
