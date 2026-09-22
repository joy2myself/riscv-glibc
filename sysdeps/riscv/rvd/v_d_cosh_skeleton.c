/* Double-precision vector cosh function.
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
  VFLOAT vx_orig = x;
  VFLOAT vx, vy, vy_special = VFMV_VF (fp_posZero, vlen);
  VBOOL special_args;

  SET_ROUNDTONEAREST;

  VUINT expo_x = __riscv_vand (__riscv_vsrl (F_AS_U (vx_orig), MAN_LEN, vlen),
                               0x7FF, vlen);

  VUINT vclass = __riscv_vfclass (vx_orig, vlen);
  VBOOL nan_inf;
  IDENTIFY (vclass, class_NaN | class_Inf, nan_inf, vlen);
  VBOOL small_x
    = __riscv_vmsleu (expo_x, EXP_BIAS - MAN_LEN - RVV_COSHF_SMALL_EXP_OFFSET,
                      vlen);
  special_args = __riscv_vmor (nan_inf, small_x, vlen);

  UINT nb_special_args = __riscv_vcpop (special_args, vlen);
  if (nb_special_args > 0)
    {
      vy_special
        = __riscv_vfmadd (special_args, vx_orig, vx_orig,
                          VFMV_VF (fp_posOne, vlen), vlen);
      vx_orig = __riscv_vfmerge (vx_orig, fp_posZero, special_args, vlen);
    }

  vx = __riscv_vfsgnj (vx_orig, fp_posOne, vlen);
  vx = __riscv_vfmin (vx, RVV_COSHF_X_MAX_D, vlen);

  VFLOAT n_flt = __riscv_vfmul (vx, RVV_LN2_INV_D, vlen);
  VINT n = __riscv_vfcvt_x (n_flt, vlen);
  n_flt = __riscv_vfcvt_f (n, vlen);
  VFLOAT r_delta = __riscv_vfnmsac (vx, RVV_LN2_HI_D, n_flt, vlen);
  VFLOAT r = __riscv_vfnmsac (r_delta, RVV_LN2_LO_D, n_flt, vlen);
  r_delta = __riscv_vfsub (r_delta, r, vlen);
  r_delta = __riscv_vfnmsac (r_delta, RVV_LN2_LO_D, n_flt, vlen);

  VFLOAT rsq = __riscv_vfmul (r, r, vlen);
  VFLOAT rcube = __riscv_vfmul (rsq, r, vlen);

  VFLOAT p_even = PSTEP (RVV_COSH_P_EVEN0_D, rsq,
                         PSTEP (RVV_COSH_P_EVEN1_D, rsq,
                                PSTEP (RVV_COSH_P_EVEN2_D, rsq,
                                       PSTEP (RVV_COSH_P_EVEN3_D,
                                              RVV_COSH_P_EVEN4_D, rsq, vlen),
                                       vlen),
                                vlen),
                         vlen);

  VFLOAT p_odd = PSTEP (RVV_COSH_P_ODD0_D, rsq,
                        PSTEP (RVV_COSH_P_ODD1_D, rsq,
                               PSTEP (RVV_COSH_P_ODD2_D, RVV_COSH_P_ODD3_D, rsq,
                                      vlen),
                               vlen),
                        vlen);

  VFLOAT p_pos = __riscv_vfmadd (p_odd, r, p_even, vlen);
  VFLOAT p_neg = __riscv_vfnmsub (p_odd, r, p_even, vlen);
  p_pos = __riscv_vfmul (p_pos, rcube, vlen);
  p_neg = __riscv_vfmul (p_neg, rcube, vlen);

  VINT m = __riscv_vmin (__riscv_vadd (n, n, vlen), RVV_COSHF_SCALE_LIMIT,
                         vlen);
  VFLOAT s = U_AS_F (__riscv_vsll (I_AS_U (__riscv_vrsub (m, EXP_BIAS, vlen)),
                                   MAN_LEN, vlen));

  VFLOAT poly = __riscv_vfnmsac (p_pos, s, p_neg, vlen);

  VFLOAT r_half = __riscv_vfmul (r, fp_posHalf, vlen);
  VFLOAT b_plus = __riscv_vfmadd (r, r_half, r, vlen);
  VFLOAT c_plus
    = __riscv_vfmacc (__riscv_vfsub (r, b_plus, vlen), r, r_half, vlen);
  VFLOAT delta_b_plus = __riscv_vfmadd (r, r_delta, r_delta, vlen);
  c_plus = __riscv_vfadd (c_plus, delta_b_plus, vlen);

  VFLOAT b_minus = __riscv_vfnmsub (r, r_half, r, vlen);
  VFLOAT c_minus
    = __riscv_vfnmsac (__riscv_vfsub (r, b_minus, vlen), r, r_half, vlen);
  VFLOAT delta_b_minus = __riscv_vfnmsub (r, r_delta, r_delta, vlen);
  c_minus = __riscv_vfadd (c_minus, delta_b_minus, vlen);

  VFLOAT b = __riscv_vfnmsub (b_minus, s, b_plus, vlen);
  VFLOAT c = __riscv_vfnmsac (__riscv_vfsub (b_plus, b, vlen), s, b_minus,
                              vlen);
  c = __riscv_vfadd (c, __riscv_vfnmsub (c_minus, s, c_plus, vlen), vlen);

  VBOOL n_large = __riscv_vmsge (n, RVV_COSHF_LARGE_N, vlen);
  VFLOAT s_hi = s;
  VFLOAT s_lo = VFMV_VF (0.0, vlen);
  s_hi = __riscv_vfmerge (s_hi, fp_posZero, n_large, vlen);
  s_lo = __riscv_vmerge (s_lo, s, n_large, vlen);

  VFLOAT a = __riscv_vfadd (s_hi, fp_posOne, vlen);
  c = __riscv_vfadd (c, s_lo, vlen);

  VFLOAT z_hi, z_lo;
  FAST2SUM (b, poly, z_hi, z_lo, vlen);
  c = __riscv_vfadd (c, z_lo, vlen);
  b = z_hi;
  FAST2SUM (a, b, z_hi, z_lo, vlen);
  c = __riscv_vfadd (c, z_lo, vlen);
  vy = __riscv_vfadd (z_hi, c, vlen);

  n = __riscv_vsub (n, 1, vlen);
  FAST_LDEXP (vy, n, vlen);

  vy = __riscv_vmerge (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
