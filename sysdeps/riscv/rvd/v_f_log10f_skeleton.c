/* Single-precision vector log10f function.
   log10f(x) computes the base-10 logarithm of x

   Based on logf implementation with log10(e) scaling.

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

/* Shared float-float helpers live in v_math.h.  */

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT vx, vy, vy_special;
  VBOOL special_args;
  VINT n_adjust;

  SET_ROUNDTONEAREST;
  vx = x;

  /* Handle special cases */
  VUINT vclass = __riscv_vfclass(vx, vlen);
  IDENTIFY(vclass, 0x3BF, special_args, vlen);
  UINT nb_special_args = __riscv_vcpop(special_args, vlen);
  n_adjust = VMVI_VX(0, vlen);
  vy_special = VFMV_VF(0.0f, vlen);

  if (nb_special_args > 0)
    {
      VBOOL id_mask;
      IDENTIFY(vclass, class_negative, id_mask, vlen);
      vx = VFMERGE_VFM(vx, fp_sNaN, id_mask, vlen);
      IDENTIFY(vclass, class_posZero, id_mask, vlen);
      vx = VFMERGE_VFM(vx, fp_negZero, id_mask, vlen);
      IDENTIFY(vclass, 0x39F, special_args, vlen);
      vy_special = VFREC7_M(special_args, vx, vlen);
      vy_special = VFADD_VV_M(special_args, vy_special, vx, vlen);
      vx = VFMERGE_VFM(vx, fp_posOne, special_args, vlen);
      IDENTIFY(vclass, class_posDenorm, id_mask, vlen);
      n_adjust = VMERGE_VXM_I(n_adjust, 64, id_mask, vlen);
      VFLOAT vx_normalized
          = VFMUL_VF_M (id_mask, vx, RVV_LOG10F_CLEAN_SUBNORMAL_SCALE_F,
                        vlen);
      vx = VMERGE_VVM_F(vx, vx_normalized, id_mask, vlen);
    }

  /* Range reduction: multiply by 4/3 */
  VFLOAT vx_adj = __riscv_vfmul (vx, RVV_LOG10F_CLEAN_INV_075_F, vlen);
  VUINT ix_adj = F_AS_U(vx_adj);

  /* Extract exponent e */
  VINT e = U_AS_I(__riscv_vsrl(ix_adj, MAN_LEN, vlen));
  e = __riscv_vand(e, 0xff, vlen);
  e = __riscv_vsub(e, EXP_BIAS, vlen);

  /* Compute m = vx * 2^(-e) */
  VINT neg_e = __riscv_vneg(e, vlen);
  VFLOAT m = I_AS_F(__riscv_vadd(F_AS_I(vx), __riscv_vsll(neg_e, MAN_LEN, vlen), vlen));
  e = __riscv_vsub (e, n_adjust, vlen);

  /* Float-float division: x = (m-1)/(m+1) as (x_hi, x_lo).  Keep the
     source-shaped DF association here; the surrounding EH envelope is the
     clean adaptation, not a reason to change the finite kernel.  */
  VFLOAT num_hi, num_lo;
  VFLOAT den_hi, den_lo;
  RVV_SLEEF_ORIG_DFADD2_VF_VF (VFMV_VF (-1.0f, vlen), m,
                               num_hi, num_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF_VF (VFMV_VF (1.0f, vlen), m,
                               den_hi, den_lo, vlen);

  VFLOAT x_hi, x_lo;
  RVV_SLEEF_ORIG_DFDIV_VF2_VF2 (num_hi, num_lo, den_hi, den_lo,
                                x_hi, x_lo, vlen);

  /* x2 = x_hi^2, poly = C2 + x2*(C1 + x2*C0), poly_tail = x2 * x_hi * poly */
  VFLOAT x2 = __riscv_vfmul(x_hi, x_hi, vlen);

  VFLOAT poly = VFMV_VF (RVV_LOG10F_CLEAN_POLY0_F, vlen);
  poly = __riscv_vfmacc (VFMV_VF (RVV_LOG10F_CLEAN_POLY1_F, vlen), poly, x2,
                         vlen);
  poly = __riscv_vfmacc (VFMV_VF (RVV_LOG10F_CLEAN_POLY2_F, vlen), poly, x2,
                         vlen);

  VFLOAT poly_tail = __riscv_vfmul(__riscv_vfmul(x2, x_hi, vlen), poly, vlen);

  VFLOAT ef = __riscv_vfcvt_f(e, vlen);

  /* Compute log10(x) in double-float: s = e*log10(2) +
     x*(2*log10(e)) + poly_tail.  These pair operations retain the pinned
     SLEEF operation order while the caller-facing EH remains clean-owned.  */
  VFLOAT x_scale_hi1, x_scale_lo1;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (x_hi, x_lo, VFMV_VF (RVV_LOG10F_CLEAN_TWO_LOG10_E_HI_F, vlen),
     VFMV_VF (RVV_LOG10F_CLEAN_TWO_LOG10_E_LO_F, vlen), x_scale_hi1,
     x_scale_lo1,
     vlen);

  VFLOAT s_hi1, s_lo1;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (VFMV_VF (RVV_LOG10F_CLEAN_LOG10_2_HI_F, vlen),
     VFMV_VF (RVV_LOG10F_CLEAN_LOG10_2_LO_F, vlen), ef, s_hi1, s_lo1, vlen);

  VFLOAT sum_hi, sum_lo;
  RVV_SLEEF_ORIG_DFADD_VF2_VF2
    (s_hi1, s_lo1, x_scale_hi1, x_scale_lo1, sum_hi, sum_lo, vlen);
  RVV_SLEEF_ORIG_DFADD_VF2_VF2_VF
    (sum_hi, sum_lo, poly_tail, sum_hi, sum_lo, vlen);

  vy = __riscv_vfadd (sum_hi, sum_lo, vlen);

  vy = VMERGE_VVM_F(vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
