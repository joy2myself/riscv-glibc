/* Single-precision vector logp1f function.
   logp1f(x) computes log(1 + x) accurately for small x

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

#include <vecmath_config.h>
#include <v_math.h>

#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

/* Match libm reference: logf(1.0 + x) with float rounding.  */
static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT vx = x;
  VFLOAT vy, vy_special;
  VBOOL special_args;
  VINT n_adjust, e_adjust;

  SET_ROUNDTONEAREST;

  vx = __riscv_vfadd (vx, 1.0f, vlen);

  VUINT vclass = __riscv_vfclass (vx, vlen);
  IDENTIFY (vclass, 0x3BF, special_args, vlen);
  UINT nb_special_args = __riscv_vcpop (special_args, vlen);
  n_adjust = VMVI_VX_F (0, vlen);
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
      IDENTIFY (vclass, class_posDenorm, id_mask, vlen);
      n_adjust = VMERGE_VXM_I (n_adjust, 64, id_mask, vlen);
      VFLOAT vx_normalized = VFMUL_VF_M
        (id_mask, vx, RVV_LOGP1F_SCALE_UP_F, vlen);
      vx = VMERGE_VVM_F (vx, vx_normalized, id_mask, vlen);
    }

  /* Scale large finite rounded inputs to avoid spurious overflow in
     vx * (4/3). */
  VBOOL adjust_large = __riscv_vmfge (vx, RVV_LOGP1F_BOUND_F, vlen);
  e_adjust = VMVI_VX_F (0, vlen);
  e_adjust = __riscv_vmerge (e_adjust, 4, adjust_large, vlen);
  vx = __riscv_vfmul (
      vx,
      __riscv_vfmerge (VFMV_VF_F (1.0f, vlen), RVV_LOGP1F_SCALE_DOWN_F,
                       adjust_large, vlen),
      vlen);

  VFLOAT vx_adj = __riscv_vfmul (vx, RVV_LOGP1F_INV_075_F, vlen);
  VUINT ix_adj = F_AS_U_F (vx_adj);

  VINT e = U_AS_I_F (__riscv_vsrl (ix_adj, MAN_LEN_F, vlen));
  e = __riscv_vand (e, 0xff, vlen);
  e = __riscv_vsub (e, EXP_BIAS_F, vlen);
  e = __riscv_vadd (e, e_adjust, vlen);

  VINT e_core = __riscv_vsub (e, e_adjust, vlen);
  VINT neg_e = __riscv_vneg (e_core, vlen);
  VFLOAT m = I_AS_F_F (
      __riscv_vadd (F_AS_I_F (vx), __riscv_vsll (neg_e, MAN_LEN_F, vlen), vlen));
  e = __riscv_vsub (e, n_adjust, vlen);

  VFLOAT num_hi, num_lo;
  VFLOAT den_hi, den_lo;
  DFADD2_VF2_VF_S (m, VFMV_VF_F (0.0f, vlen), -1.0f, num_hi, num_lo, vlen);
  DFADD2_VF2_VF_S (m, VFMV_VF_F (0.0f, vlen), 1.0f, den_hi, den_lo, vlen);

  VFLOAT x_hi, x_lo;
  DFDIV_VF2_VF2_VF2 (num_hi, num_lo, den_hi, den_lo, x_hi, x_lo, vlen);

  VFLOAT x2 = __riscv_vfmul (x_hi, x_hi, vlen);

  VFLOAT poly = VFMV_VF_F (RVV_SLEEF_LOG_U1_POLY0_F, vlen);
  poly = __riscv_vfmacc
    (VFMV_VF_F (RVV_SLEEF_LOG_U1_POLY1_F, vlen), poly, x2, vlen);
  poly = __riscv_vfmacc
    (VFMV_VF_F (RVV_SLEEF_LOG_U1_POLY2_F, vlen), poly, x2, vlen);

  VFLOAT poly_tail = __riscv_vfmul (__riscv_vfmul (x2, x_hi, vlen), poly, vlen);

  VFLOAT ef = __riscv_vfcvt_f (e, vlen);

  VFLOAT s_hi, s_lo;
  DFMUL_VF2_VF2_VF (VFMV_VF_F (RVV_LOGP1F_LN2_HI_F, vlen),
                    VFMV_VF_F (RVV_LOGP1F_LN2_LO_F, vlen), ef,
                    s_hi, s_lo, vlen);

  VFLOAT two_x_hi, two_x_lo;
  DFSCALE_VF2_VF2_S (x_hi, x_lo, 2.0f, two_x_hi, two_x_lo, vlen);
  DFADD_VF2_VF2_VF2 (s_hi, s_lo, two_x_hi, two_x_lo, s_hi, s_lo, vlen);
  DFADD_VF2_VF2_VF (s_hi, s_lo, poly_tail, s_hi, s_lo, vlen);
  vy = __riscv_vfadd (s_hi, s_lo, vlen);

  vy = VMERGE_VVM_F (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
