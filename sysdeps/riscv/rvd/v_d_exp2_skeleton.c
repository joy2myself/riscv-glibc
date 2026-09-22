/* Double-precision vector exp2 function.
   exp2(x) computes 2^x

   Based on SLEEF implementation adapted for double precision.

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

/* Keep the source vpow2i/vldexp2 operation order local to exp2.  The
   public implementation still uses the clean EH shell below; these helpers
   only replace the clean-only FAST_LDEXP finite reconstruction.  */
static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (exp2_source_vpow2i, LMUL) (VINT q, size_t vlen)
{
  VINT biased = __riscv_vadd (q, VMVI_VX (EXP_BIAS, vlen), vlen);
  return I_AS_F (__riscv_vsll (biased, MAN_LEN, vlen));
}

static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (exp2_source_vldexp2, LMUL)
  (VFLOAT d, VINT q, size_t vlen)
{
  VINT q_half = __riscv_vsra (q, 1, vlen);
  VINT q_rest = __riscv_vsub (q, q_half, vlen);
  VFLOAT p_half = RVV_LIBMVEC_IMPL_NAME
    (exp2_source_vpow2i, LMUL) (q_half, vlen);
  VFLOAT p_rest = RVV_LIBMVEC_IMPL_NAME
    (exp2_source_vpow2i, LMUL) (q_rest, vlen);
  return VFMUL_VV (VFMUL_VV (d, p_half, vlen), p_rest, vlen);
}

/* SLEEF xexp2 finite core with the existing RVV special-lane
   sanitization retained around it.  */

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT vx = x;
  VFLOAT vy, vy_special;
  VBOOL special_args, subnormal_args, front_args;

  SET_ROUNDTONEAREST;

  /* Set results for input of NaN and Inf; substitute them with zero */
  vy_special = VFMV_VF (fp_posZero, vlen);
  VUINT vclass = __riscv_vfclass (vx, vlen);
  IDENTIFY (vclass, class_NaN | class_Inf, special_args, vlen);
  IDENTIFY (vclass, class_Subnormal, subnormal_args, vlen);
  front_args = VMOR (special_args, subnormal_args, vlen);

  UINT nb_front_args = __riscv_vcpop (front_args, vlen);
  if (nb_front_args > 0)
    {
      /* Substitute -Inf with +0 */
      VBOOL id_mask;
      IDENTIFY (vclass, class_negInf, id_mask, vlen);
      vx = __riscv_vfmerge (vx, fp_posZero, id_mask, vlen);
      vy_special = __riscv_vfadd (special_args, vx, vx, vlen);

      /* exp2(x) rounds to 1 for every binary64 subnormal x, while the
         scalar operation still raises FE_INEXACT.  Do not send subnormal
         lanes through s^2/s^4/s^8: those products would spuriously raise
         FE_UNDERFLOW.  The controlled 1+x provider preserves the scalar
         inexact side effect before the lane is sanitized out of the core.  */
      VFLOAT subnormal_result = __riscv_vfadd
        (subnormal_args, VFMV_VF (fp_posOne, vlen), vx, vlen);
      vy_special = __riscv_vmerge
        (vy_special, subnormal_result, subnormal_args, vlen);

      vx = __riscv_vfmerge (vx, fp_posZero, front_args, vlen);
    }

  VBOOL overflow = __riscv_vmfge (vx, RVV_EXP2_OVERFLOW_BOUND_D, vlen);
  VBOOL underflow = __riscv_vmflt (vx, RVV_EXP2_UNDERFLOW_BOUND_D, vlen);

  VFLOAT vx_core
      = FCLIP (vx, RVV_EXP2_UNDERFLOW_BOUND_D, RVV_EXP2_OVERFLOW_BOUND_D,
               vlen);
  /* SLEEF performs vrint_vd_vd first and then vrint_vi_vd.  Keep both
     conversions explicit even though the enclosing RNE mode makes the
     second conversion exact for this bounded range.  */
  VINT rounded_once = __riscv_vfcvt_x (vx_core, vlen);
  VFLOAT flt_n = __riscv_vfcvt_f (rounded_once, vlen);
  VINT n = __riscv_vfcvt_x (flt_n, vlen);
  VFLOAT s = __riscv_vfsub (vx_core, flt_n, vlen);

  VFLOAT s2 = __riscv_vfmul (s, s, vlen);
  VFLOAT s4 = __riscv_vfmul (s2, s2, vlen);
  VFLOAT s8 = __riscv_vfmul (s4, s4, vlen);

  VFLOAT p_hi
      = PSTEP (RVV_EXP2_POLY_HI0_D, s,
               VFMV_VF (RVV_EXP2_POLY_HI1_D, vlen), vlen);

  VFLOAT p_mid_hi_a
      = PSTEP (RVV_EXP2_POLY_MID_HI_A0_D, s,
               VFMV_VF (RVV_EXP2_POLY_MID_HI_A1_D, vlen), vlen);
  VFLOAT p_mid_hi_b
      = PSTEP (RVV_EXP2_POLY_MID_HI_B0_D, s,
               VFMV_VF (RVV_EXP2_POLY_MID_HI_B1_D, vlen), vlen);
  VFLOAT p_mid_hi
      = VFMADD_VV (s2, p_mid_hi_a, p_mid_hi_b, vlen);

  VFLOAT p_mid_lo_a
      = PSTEP (RVV_EXP2_POLY_MID_LO_A0_D, s,
               VFMV_VF (RVV_EXP2_POLY_MID_LO_A1_D, vlen), vlen);
  VFLOAT p_mid_lo_b
      = PSTEP (RVV_EXP2_POLY_MID_LO_B0_D, s,
               VFMV_VF (RVV_EXP2_POLY_MID_LO_B1_D, vlen), vlen);
  VFLOAT p_mid_lo
      = VFMADD_VV (s2, p_mid_lo_a, p_mid_lo_b, vlen);

  VFLOAT poly8 = VFMADD_VV (s4, p_mid_hi, p_mid_lo, vlen);
  vy = VFMADD_VV (s8, p_hi, poly8, vlen);
  vy = VFMADD_VV (vy, s, VFMV_VF (RVV_LN2_HI_D, vlen), vlen);
  vy = VFMADD_VV (vy, s, VFMV_VF (fp_posOne, vlen), vlen);

  vy = RVV_LIBMVEC_IMPL_NAME
    (exp2_source_vldexp2, LMUL) (vy, n, vlen);
  vy = __riscv_vfmerge (vy, fp_posInf, overflow, vlen);
  vy = __riscv_vfmerge (vy, fp_posZero, underflow, vlen);

  /* Incorporate results of exceptional inputs */
  vy = __riscv_vmerge (vy, vy_special, front_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
