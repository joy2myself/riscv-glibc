/* Single-precision vector exp2 function.
   exp2f(x) computes 2^x

   Based on SLEEF implementation.

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

/* Keep the source xexp2f vpow2i/vldexp2 reconstruction in the finite core.
   The surrounding classify/sanitize/final-merge path remains the clean EH
   adaptation; only the clean-only FAST_LDEXP operation sequence is removed.  */
static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (exp2f_source_vpow2i, LMUL) (VINT q, size_t vlen)
{
  VINT biased = __riscv_vadd (q, VMVI_VX (EXP_BIAS_F, vlen), vlen);
  return I_AS_F (__riscv_vsll (biased, MAN_LEN_F, vlen));
}

static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (exp2f_source_vldexp2, LMUL)
  (VFLOAT d, VINT q, size_t vlen)
{
  VINT q_half = __riscv_vsra (q, 1, vlen);
  VINT q_rest = __riscv_vsub (q, q_half, vlen);
  VFLOAT p_half = RVV_LIBMVEC_IMPL_NAME
    (exp2f_source_vpow2i, LMUL) (q_half, vlen);
  VFLOAT p_rest = RVV_LIBMVEC_IMPL_NAME
    (exp2f_source_vpow2i, LMUL) (q_rest, vlen);
  return VFMUL_VV (VFMUL_VV (d, p_half, vlen), p_rest, vlen);
}

/* This version computes 2^x directly using argument reduction.
   For x = n + r where n is integer and |r| <= 0.5,
   2^x = 2^n * 2^r. */

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT vx = x;
  VFLOAT vy, vy_special;
  VBOOL special_args, range_args, overflow, underflow;

  SET_ROUNDTONEAREST;

  VFLOAT vx_orig = vx;
  vy_special = VFMV_VF (fp_posZero, vlen);

  /* Set results for input of NaN and Inf; substitute them with zero */
  VUINT vclass = __riscv_vfclass (vx, vlen);
  IDENTIFY (vclass, class_NaN | class_Inf, special_args, vlen);

  UINT nb_special_args = __riscv_vcpop (special_args, vlen);
  if (nb_special_args > 0)
    {
      /* Substitute -Inf with +0 */
      VBOOL id_mask;
      IDENTIFY (vclass, class_negInf, id_mask, vlen);
      vx = __riscv_vfmerge (vx, fp_posZero, id_mask, vlen);
      vy_special = __riscv_vfadd (special_args, vx, vx, vlen);
      vx = __riscv_vfmerge (vx, fp_posZero, special_args, vlen);
    }

  VFLOAT vx_cmp = __riscv_vfmerge (vx_orig, fp_posZero, special_args, vlen);
  overflow = __riscv_vmfge (vx_cmp, RVV_EXP2F_CLEAN_MAX_INPUT_F, vlen);
  underflow = __riscv_vmflt (vx_cmp, RVV_EXP2F_CLEAN_MIN_INPUT_F, vlen);
  range_args = __riscv_vmor (overflow, underflow, vlen);
  vy_special = __riscv_vfmerge (vy_special, fp_posInf, overflow, vlen);
  vx = __riscv_vfmerge (vx, fp_posZero, range_args, vlen);

  /* Argument reduction: 2^x = 2^n * 2^r where n is integer, |r| <= 0.5 */
  VFLOAT u = __riscv_vfcvt_f (__riscv_vfcvt_x (vx, vlen), vlen);
  VINT n = __riscv_vfcvt_x (u, vlen);
  VFLOAT s = __riscv_vfsub (vx, u, vlen);

  /* Polynomial approximation - SLEEF xexp2f u10 */
  u = VFMV_VF (RVV_EXP2F_CLEAN_POLY0_F, vlen);
  u = __riscv_vfmadd (u, s, VFMV_VF (RVV_EXP2F_CLEAN_POLY1_F, vlen), vlen);
  u = __riscv_vfmadd (u, s, VFMV_VF (RVV_EXP2F_CLEAN_POLY2_F, vlen), vlen);
  u = __riscv_vfmadd (u, s, VFMV_VF (RVV_EXP2F_CLEAN_POLY3_F, vlen), vlen);
  u = __riscv_vfmadd (u, s, VFMV_VF (RVV_EXP2F_CLEAN_POLY4_F, vlen), vlen);
  u = __riscv_vfmadd (u, s, VFMV_VF (RVV_EXP2F_CLEAN_POLY5_F, vlen), vlen);

  vy = __riscv_vfmadd (u, s, VFMV_VF (fp_posOne, vlen), vlen);

  /* Scale by 2^n with the source's split exponent and two multiplies.  */
  vy = RVV_LIBMVEC_IMPL_NAME (exp2f_source_vldexp2, LMUL) (vy, n, vlen);

  /* Keep finite range exceptions visible even though the unsafe input lanes
     were sanitized before the finite core.  The probe values are discarded
     by the final special merge; only their fflags effects are retained.  */
  VBOOL special_total = __riscv_vmor (special_args, range_args, vlen);
  UINT nb_special_total = __riscv_vcpop (special_total, vlen);
  if (nb_special_total > 0)
    {
      VFLOAT overflow_probe = VFMV_VF (RVV_EXP2F_CLEAN_MAX_FINITE_F, vlen);
      overflow_probe = VFMUL_VF_M (overflow, overflow_probe, 2.0f, vlen);
      vy = VMERGE_VVM_F (vy, overflow_probe, overflow, vlen);

      /* Use a genuinely inexact tiny product so the underflow flag is
         visible on implementations with tininess-after-rounding.  */
      VFLOAT underflow_probe = VFMV_VF (RVV_EXP2F_CLEAN_MIN_NORMAL_F, vlen);
      underflow_probe = VFMUL_VF_M (underflow, underflow_probe, 0.3f, vlen);
      vy = VMERGE_VVM_F (vy, underflow_probe, underflow, vlen);
    }

  /* Incorporate results of exceptional inputs */
  vy = __riscv_vmerge (vy, vy_special, special_total, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
