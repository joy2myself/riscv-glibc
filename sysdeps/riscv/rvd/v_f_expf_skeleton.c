/* Single-precision vector exp function.
   expf(x) computes e^x

   Based on SLEEF xexpf implementation (u10 precision).

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

/* Keep the source xexpf vpow2i/vldexp2 reconstruction in the finite core.
   The surrounding classify/sanitize/final-merge path remains the clean EH
   adaptation; only the clean-only FAST_LDEXP operation sequence is removed.  */
static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (expf_source_vpow2i, LMUL) (VINT e, size_t vlen)
{
  VINT biased = __riscv_vadd (e, EXP_BIAS_F, vlen);
  return I_AS_F (__riscv_vsll (biased, MAN_LEN_F, vlen));
}

static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (expf_source_vldexp2, LMUL)
  (VFLOAT d, VINT e, size_t vlen)
{
  VINT e_half = __riscv_vsra (e, 1, vlen);
  VINT e_rest = __riscv_vsub (e, e_half, vlen);
  VFLOAT p_half = RVV_LIBMVEC_IMPL_NAME
    (expf_source_vpow2i, LMUL) (e_half, vlen);
  VFLOAT p_rest = RVV_LIBMVEC_IMPL_NAME
    (expf_source_vpow2i, LMUL) (e_rest, vlen);
  return VFMUL_VV (VFMUL_VV (d, p_half, vlen), p_rest, vlen);
}

/* This version uses SLEEF-style argument reduction to [-log(2)/2, log(2)/2]
   and computes exp using a simplified polynomial for single precision. */

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT vx = x;
  VFLOAT vy, vy_special;
  VBOOL special_args, range_args, overflow, underflow;

  SET_ROUNDTONEAREST;

  VFLOAT vx_orig = vx;
  vy_special = VFMV_VF (0.0f, vlen);

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
  overflow = __riscv_vmfgt (vx_cmp, RVV_EXPF_CLEAN_MAX_INPUT_F, vlen);
  underflow = __riscv_vmflt (vx_cmp, RVV_EXPF_CLEAN_MIN_INPUT_F, vlen);
  range_args = __riscv_vmor (overflow, underflow, vlen);
  vy_special = __riscv_vfmerge (vy_special, fp_posInf, overflow, vlen);
  vx = __riscv_vfmerge (vx, fp_posZero, range_args, vlen);

  /* Argument reduction: exp(x) = 2^n * exp(r) where |r| <= log(2)/2 */
  VFLOAT flt_n = __riscv_vfmul (vx, RVV_EXPF_CLEAN_R_LN2_F, vlen);
  VINT n = __riscv_vfcvt_x (flt_n, vlen);
  flt_n = __riscv_vfcvt_f (n, vlen);
  VFLOAT r = __riscv_vfnmsac (vx, RVV_EXPF_CLEAN_L2U_F, flt_n, vlen);
  r = __riscv_vfnmsac (r, RVV_EXPF_CLEAN_L2L_F, flt_n, vlen);

  /* Polynomial computation - SLEEF xexpf u10 (6 coefficients) */
  VFLOAT poly = VFMV_VF (RVV_EXPF_CLEAN_POLY0_F, vlen);
  poly = __riscv_vfmadd
    (poly, r, VFMV_VF (RVV_EXPF_CLEAN_POLY1_F, vlen), vlen);
  poly = __riscv_vfmadd
    (poly, r, VFMV_VF (RVV_EXPF_CLEAN_POLY2_F, vlen), vlen);
  poly = __riscv_vfmadd
    (poly, r, VFMV_VF (RVV_EXPF_CLEAN_POLY3_F, vlen), vlen);
  poly = __riscv_vfmadd
    (poly, r, VFMV_VF (RVV_EXPF_CLEAN_POLY4_F, vlen), vlen);
  poly = __riscv_vfmadd (poly, r, VFMV_VF (RVV_EXPF_CLEAN_POLY5_F, vlen),
                         vlen);

  VFLOAT r_sq = __riscv_vfmul (r, r, vlen);
  vy = __riscv_vfmadd (r_sq, poly, r, vlen);
  vy = __riscv_vfadd (vy, fp_posOne, vlen);

  /* at this point, vy is the polynomial approximation
     vy ~=~ exp(r) */

  /* Need to compute 2^n * exp(r). */
  vy = RVV_LIBMVEC_IMPL_NAME (expf_source_vldexp2, LMUL) (vy, n, vlen);

  /* Keep finite range exceptions visible even though the unsafe input lanes
     were sanitized before the finite core.  The probe values are discarded
     by the final special merge; only their fflags effects are retained.  */
  VBOOL special_total = __riscv_vmor (special_args, range_args, vlen);
  UINT nb_special_total = __riscv_vcpop (special_total, vlen);
  if (nb_special_total > 0)
    {
      VFLOAT overflow_probe = VFMV_VF (RVV_EXPF_CLEAN_MAX_FINITE_F, vlen);
      overflow_probe = VFMUL_VF_M (overflow, overflow_probe, 2.0f, vlen);
      vy = VMERGE_VVM_F (vy, overflow_probe, overflow, vlen);

      /* Use a genuinely inexact tiny product so the underflow flag is
         visible on implementations with tininess-after-rounding.  */
      VFLOAT underflow_probe = VFMV_VF (RVV_EXPF_CLEAN_MIN_NORMAL_F, vlen);
      underflow_probe = VFMUL_VF_M (underflow, underflow_probe, 0.3f, vlen);
      vy = VMERGE_VVM_F (vy, underflow_probe, underflow, vlen);
    }

  /* Incorporate results of exceptional inputs.  */
  vy = __riscv_vmerge (vy, vy_special, special_total, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
