/* Single-precision vector log function.
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

#ifndef V_F_LOGF_VFREC7_OPT_SKELETON_ONCE
#define V_F_LOGF_VFREC7_OPT_SKELETON_ONCE

#endif

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
    VFLOAT vx_normalized = VFMUL_VF_M(id_mask, vx, 0x1.0p64f, vlen);
    vx = VMERGE_VVM_F(vx, vx_normalized, id_mask, vlen);
  }

  /* Range reduction: multiply by 4/3 */
  VFLOAT vx_adj = __riscv_vfmul(vx, 0x1.555556p0f, vlen);
  VUINT ix_adj = F_AS_U(vx_adj);
  
  /* Extract exponent e - keep as integer, don't convert to float yet! */
  VINT e = U_AS_I(__riscv_vsrl(ix_adj, MAN_LEN, vlen));
  e = __riscv_vand(e, 0xff, vlen);
  e = __riscv_vsub(e, EXP_BIAS, vlen);
  
  /* Compute m = vx * 2^(-e) */
  VINT neg_e = __riscv_vneg(e, vlen);
  VFLOAT m = I_AS_F(__riscv_vadd(F_AS_I(vx), __riscv_vsll(neg_e, MAN_LEN, vlen), vlen));
  e = __riscv_vsub(e, n_adjust, vlen);

  /* Float-float division: x = (m-1)/(m+1) as (x_hi, x_lo) */
  VFLOAT num = __riscv_vfsub(m, 1.0f, vlen);
  VFLOAT den = __riscv_vfadd(m, 1.0f, vlen);

  VFLOAT r = __riscv_vfrec7(den, vlen);
  VFLOAT two = VFMV_VF(2.0f, vlen);
  r = __riscv_vfmul(__riscv_vfnmsac(two, den, r, vlen), r, vlen);
  two = VFMV_VF(2.0f, vlen);
  r = __riscv_vfmul(__riscv_vfnmsac(two, den, r, vlen), r, vlen);
  
  VFLOAT x_hi = __riscv_vfmul(num, r, vlen);
  VFLOAT residual = __riscv_vfnmsac(num, x_hi, den, vlen);
  VFLOAT x_lo = __riscv_vfmul(residual, r, vlen);

  /* x2 = x_hi*x_hi
     poly = C2 + x2*(C1 + x2*C0)
     poly_tail = x2 * x_hi * poly
  */
  VFLOAT x2 = __riscv_vfmul(x_hi, x_hi, vlen);
  
  VFLOAT poly = VFMV_VF(RVV_SLEEF_LOG_U1_POLY0_F, vlen);
  poly = __riscv_vfmacc(VFMV_VF(RVV_SLEEF_LOG_U1_POLY1_F, vlen), poly,
                        x2, vlen);
  poly = __riscv_vfmacc(VFMV_VF(RVV_SLEEF_LOG_U1_POLY2_F, vlen), poly,
                        x2, vlen);
  
  VFLOAT poly_tail = __riscv_vfmul(__riscv_vfmul(x2, x_hi, vlen), poly, vlen);

  VFLOAT ef = __riscv_vfcvt_f(e, vlen);

  VFLOAT s_hi = __riscv_vfmul(ef, RVV_SLEEF_LOG_U1_LN2_HI_F, vlen);
  VFLOAT s_lo = __riscv_vfmul(ef, RVV_SLEEF_LOG_U1_LN2_LO_F, vlen);

  /* Compact dfadd: s + 2*x + poly_tail
     sum_hi = s_hi + 2*x_hi
     sum_lo = (s_hi - sum_hi) + 2*x_hi + s_lo + 2*x_lo + poly_tail
  */
  VFLOAT sum_hi = __riscv_vfadd(s_hi, __riscv_vfadd(x_hi, x_hi, vlen), vlen);
  VFLOAT err = __riscv_vfsub(s_hi, sum_hi, vlen);
  VFLOAT sum_lo = __riscv_vfadd(err, __riscv_vfadd(x_hi, x_hi, vlen), vlen);
  sum_lo = __riscv_vfadd(sum_lo, s_lo, vlen);
  sum_lo = __riscv_vfadd(sum_lo, __riscv_vfadd(x_lo, x_lo, vlen), vlen);
  sum_lo = __riscv_vfadd(sum_lo, poly_tail, vlen);
  vy = __riscv_vfadd(sum_hi, sum_lo, vlen);

  vy = VMERGE_VVM_F(vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
