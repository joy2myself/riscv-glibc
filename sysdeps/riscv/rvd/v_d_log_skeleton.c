/* Double-precision vector log function.
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

#ifndef V_D_LOG_SKELETON_ONCE
#define V_D_LOG_SKELETON_ONCE

#define LOGB_2_HI RVV_LN2_HI_D
#define LOGB_2_LO RVV_LN2_LO_D

#endif

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  const struct v_log_data *d = ptr_barrier(&__v_log_data);
  VFLOAT vx, vy, vy_special;
  VBOOL special_args;
  VINT n_adjust;

  SET_ROUNDTONEAREST;
  /* stripmining over input arguments: vlen chosen by wrapper */
  vx = x;

  /* NaN, Inf, and -ve handling, as well as scaling denormal input by 2^64
   */
  VUINT vclass = __riscv_vfclass((vx), (vlen));
  /* special handling except positive normal number */
  IDENTIFY(vclass, 0x3BF, (special_args), (vlen));
  UINT nb_special_args = __riscv_vcpop((special_args), (vlen));
  n_adjust = VMVI_VX(0, (vlen));
  if (nb_special_args > 0)
    {
      VBOOL id_mask;
      /* substitute negative arguments with sNaN */
      IDENTIFY(vclass, class_negative, id_mask, (vlen));
      vx = VFMERGE_VFM(vx, fp_sNaN, id_mask, vlen);
      /* substitute +0 argument with -0 */
      IDENTIFY(vclass, class_posZero, id_mask, vlen);
      vx = VFMERGE_VFM(vx, fp_negZero, id_mask, vlen);
      /* eliminate positive denorm input from special_args */
      IDENTIFY(vclass, 0x39F, (special_args), (vlen));
      /* for narrowed set of special arguments, compute vx+vfrec7(vx) */
      vy_special = VFREC7_M((special_args), (vx), (vlen));
      vy_special
          = VFADD_VV_M((special_args), vy_special, (vx), (vlen));
      vx = VFMERGE_VFM((vx), fp_posOne, (special_args), (vlen));
      /* scale up input for positive denormals */
      IDENTIFY(vclass, class_posDenorm, id_mask, (vlen));
      n_adjust = VMERGE_VXM_I(n_adjust, 64, id_mask, vlen);
      VFLOAT vx_normalized = VFMUL_VF_M(id_mask, vx, 0x1.0p64, vlen);
      vx = VMERGE_VVM_F(vx, vx_normalized, id_mask, vlen);
    }

  /* in_arg at this point are positive, finite and not subnormal
     Decompose in_arg into n, B, r: in_arg = 2^n (1/B) (1 + r)
     B is equivalently defined by ind, 0 <= ind < 128 */
  VINT n = U_AS_I(__riscv_vadd(
      __riscv_vsrl(F_AS_U (vx), MAN_LEN - 1, vlen), 1, vlen));
  n = __riscv_vsra(n, 1, vlen);
  n = __riscv_vsub(n, EXP_BIAS, vlen);
  vx = U_AS_F(
      __riscv_vsrl(__riscv_vsll (F_AS_U (vx), BIT_WIDTH - MAN_LEN, vlen),
                    BIT_WIDTH - MAN_LEN, vlen));
  vx = U_AS_F(
      __riscv_vadd (F_AS_U (vx), (uint64_t)EXP_BIAS << MAN_LEN, vlen));
  n = __riscv_vsub (n, n_adjust, vlen);
  VFLOAT n_flt = __riscv_vfcvt_f(n, vlen);
  VFLOAT B = __riscv_vfrec7(vx, vlen);
  /* get 7 msb of mantissa, and left shift by 3 to get address */
  VUINT ind = __riscv_vand(__riscv_vsrl (F_AS_U (vx), MAN_LEN - 10, vlen),
                            0x3F8, vlen);
  /* adjust B to be 1.0 if ind == 0 */
  VBOOL adjust_B = __riscv_vmseq(ind, 0, vlen);
  B = VFMERGE_VFM(B, fp_posOne, adjust_B, vlen);
  /* finally get r = B * in_arg - 1.0 */
  VFLOAT r = VFMV_VF(fp_posOne, vlen);
  r = __riscv_vfmsac(r, vx, B, vlen);

  /* Base-B log is logB(in_arg) = logB(2^n * 1/B) + logB(1 + r)
     (n + log2(1/B))*logB(2) + log(1+r)*logB(e)
     log2(1/B) is stored in a table
     and log(1+r) is approximated by r + poly
     poly is a polynomial in r in the form r^2 * (p0 + p1 r + ... )
     To deliver this result accurately, one uses logB(2) and logB(e)
     with extra precision and sums the various terms in an appropriate
     order */
  VFLOAT rsq = __riscv_vfmul(r, r, vlen);
  VFLOAT rcube = __riscv_vfmul(rsq, r, vlen);

  VFLOAT poly_right = PSTEP(
      0x1.9999998877038p-3, r,
      PSTEP(-0x1.555c54f8b7c6cp-3, 0x1.2499765b3c27ap-3, r, vlen), vlen);

  VFLOAT poly_left = PSTEP(
      -0x1.000000000001cp-1, r,
      PSTEP(0x1.55555555555a9p-2, -0x1.fffffff2018cfp-3, r, vlen), vlen);

  VFLOAT poly = __riscv_vfmadd(poly_right, rcube, poly_left, vlen);
  poly = __riscv_vfmul(rsq, poly, vlen);
  /* log_e(1+r) is r + poly */

  /* Load table values and get n_flt + T to be A + a */
  VINT T = __riscv_vluxei64(d->table, ind, vlen);
  VINT T_hi = __riscv_vsll(__riscv_vsra (T, 24, vlen), 24, vlen);
  VINT T_lo = __riscv_vsub(T, T_hi, vlen);
  VFLOAT T_hi_flt = __riscv_vfcvt_f(T_hi, vlen);
  VFLOAT A = __riscv_vfmadd(T_hi_flt, 0x1.0p-63, n_flt, vlen);
  VFLOAT a = __riscv_vfcvt_f(T_lo, vlen);
  a = __riscv_vfmul (a, 0x1.0p-63, vlen);

  /* Compute (A + a) * (logB_2_hi + logB_2_lo) + (r + P) * (logB_e_hi +
     logB_e_lo) */
  VFLOAT delta_1 = __riscv_vfmul(A, LOGB_2_LO, vlen);
  delta_1 = __riscv_vfmadd(a, LOGB_2_HI, delta_1, vlen);
  poly = __riscv_vfadd(poly, delta_1, vlen);

  poly = __riscv_vfadd(poly, r, vlen);

  vy = __riscv_vfmadd(A, LOGB_2_HI, poly, vlen);
  vy = VMERGE_VVM_F(vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
