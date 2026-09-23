/* Double-precision vector sinh function.
   sinh(x) = (exp(x) - exp(-x)) / 2

   Based on RISC-V Vector Library for Math (rvvlm) implementation.

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

/* This version reduces argument to [-log2/2, log2/2]
   Exploits common expressions exp(R) and exp(-R), and uses purely
   floating point method to preserve precision. */

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT vx_orig = x;
  VFLOAT vx, vy, vy_special;
  VBOOL special_args;
  VUINT expo_x;

  SET_ROUNDTONEAREST;

  expo_x = __riscv_vand (__riscv_vsrl (F_AS_U (vx_orig), MAN_LEN, vlen), 0x7FF, vlen);

  /* Set results for input of NaN and Inf and also for |x| very small */
  VUINT vclass = __riscv_vfclass (vx_orig, vlen);
  VBOOL NaN_Inf;
  IDENTIFY (vclass, class_NaN | class_Inf, NaN_Inf, vlen);
  VBOOL small_x = __riscv_vmsleu (expo_x, EXP_BIAS - MAN_LEN - 5, vlen);
  special_args = __riscv_vmor (NaN_Inf, small_x, vlen);

  vy_special = VFMV_VF (fp_posZero, vlen);
  UINT nb_special_args = __riscv_vcpop (special_args, vlen);
  if (nb_special_args > 0)
    {
      vy_special = __riscv_vfmadd (special_args, vx_orig, 0x1.0p-60, vx_orig,
                                   vlen);
      vx_orig = __riscv_vfmerge (vx_orig, fp_posZero, special_args, vlen);
    }

  /* Both sinh and cosh have sign symmetry; suffices to work on |x|.
     For sinh(x) = sign(x) * sinh(|x|). */
  vx = __riscv_vfsgnj (vx_orig, fp_posOne, vlen);

  /* Suffices to clip |x| to 714.0, which is bigger than 1030 log(2) */
  vx = __riscv_vfmin (vx, 0x1.65p9, vlen);

  /* Argument reduction: exp(x) = 2^n * exp(r) where |r| <= log(2)/2 */
  VFLOAT n_flt = __riscv_vfmul (vx, RVV_LN2_INV_D, vlen);
  VINT n = __riscv_vfcvt_x (n_flt, vlen);
  n_flt = __riscv_vfcvt_f (n, vlen);
  VFLOAT r_delta = __riscv_vfnmsac (vx, RVV_LN2_HI_D, n_flt, vlen);
  VFLOAT r = __riscv_vfnmsac (r_delta, RVV_LN2_LO_D, n_flt, vlen);
  r_delta = __riscv_vfsub (r_delta, r, vlen);
  r_delta = __riscv_vfnmsac (r_delta, RVV_LN2_LO_D, n_flt, vlen);

  /* At this point exp(x) = 2^n exp(r'), where r' = r + delta_r
     sinh(x) is 2^(n-1) ( exp(r') - 2^(-2n) exp(-r') )
     Note that n >= 0. Moreover, the factor 2^(-2n) can be replaced by
     s = 2^(-m), m = min(2n, 60)
     sinh(x) = 2^(n-1)(exp(r') - s exp(-r')) */

  /* exp(r') and exp(-r') will be computed purely in floating point
     using extra-precision simulation when needed
     Note exp(t) is approximated by
           1 + t + t^2/2 + t^3(p_even(t^2) + t*p_odd(t^2))
     and thus exp(-t) is approximated
           1 - t + t^2/2 - t^3(p_even(t^2) - t*p_odd(t^2))
     So we compute the common expressions p_even and p_odd separately.
     Moreover, they can be evaluated as r*r alone, not needing r_delta
     because they are at least a factor of (log(2)/2)^2/6 smaller than the
     final result of interest. */

  VFLOAT rsq = __riscv_vfmul (r, r, vlen);
  VFLOAT rcube = __riscv_vfmul (rsq, r, vlen);

  VFLOAT p_even = PSTEP (0x1.555555555555ap-3, rsq,
                        PSTEP (0x1.111111110ef6ap-7, rsq,
                              PSTEP (0x1.a01a01b32b633p-13, rsq,
                                    PSTEP (0x1.71ddef82f4beep-19,
                                          0x1.af6eacd796f0bp-26, rsq, vlen),
                                    vlen),
                              vlen),
                        vlen);

  VFLOAT p_odd = PSTEP (
      0x1.5555555553aefp-5, rsq,
      PSTEP (0x1.6c16c17a09506p-10, rsq,
            PSTEP (0x1.a019b37a2b3dfp-16, 0x1.289788d8bdadfp-22, rsq, vlen),
            vlen),
      vlen);

  VFLOAT p_pos = __riscv_vfmadd (p_odd, r, p_even, vlen);
  VFLOAT p_neg = __riscv_vfnmsub (p_odd, r, p_even, vlen);
  p_pos = __riscv_vfmul (p_pos, rcube, vlen);
  p_neg = __riscv_vfmul (p_neg, rcube, vlen);

  /* exp( r') is approximated by 1 + r' + (r')^2/2 + p_pos
     exp(-r') is approximated by 1 - r' + (r')^2/2 - p_neg */

  VINT m = __riscv_vmin (__riscv_vadd (n, n, vlen), 60, vlen);
  VFLOAT s = U_AS_F (
      __riscv_vsll (I_AS_U (__riscv_vrsub (m, EXP_BIAS, vlen)), MAN_LEN, vlen));

  /* For sinh: poly = p_pos + s * p_neg */
  VFLOAT poly = __riscv_vfmacc (p_pos, s, p_neg, vlen);

  /* sinh = (1 - s) + ([r' + (r')^2/2] + s [r' - (r')^2/2]) + poly
     We need r' + (r')^2/2 and r' - (r')^2/2 and their sum/diff to high precision
     and 1 - s to high precision */
  VFLOAT r_half = __riscv_vfmul (r, 0x1.0p-1, vlen);
  VFLOAT B_plus = __riscv_vfmadd (r, r_half, r, vlen);
  VFLOAT b_plus =
      __riscv_vfmacc (__riscv_vfsub (r, B_plus, vlen), r, r_half, vlen);
  VFLOAT delta_b_plus = __riscv_vfmadd (r, r_delta, r_delta, vlen);
  b_plus = __riscv_vfadd (b_plus, delta_b_plus, vlen);

  VFLOAT B_minus = __riscv_vfnmsub (r, r_half, r, vlen);
  VFLOAT b_minus =
      __riscv_vfnmsac (__riscv_vfsub (r, B_minus, vlen), r, r_half, vlen);
  VFLOAT delta_b_minus = __riscv_vfnmsub (r, r_delta, r_delta, vlen);
  b_minus = __riscv_vfadd (b_minus, delta_b_minus, vlen);

  /* For sinh: B = B_plus + s * B_minus */
  VFLOAT B = __riscv_vfmadd (B_minus, s, B_plus, vlen);
  VFLOAT b = __riscv_vfmacc (__riscv_vfsub (B_plus, B, vlen), s, B_minus, vlen);
  b = __riscv_vfadd (b, __riscv_vfmadd (b_minus, s, b_plus, vlen), vlen);

  VBOOL n_large = __riscv_vmsge (n, 50, vlen);
  VFLOAT s_hi = s;
  VFLOAT s_lo = VFMV_VF (0.0, vlen);
  s_hi = __riscv_vfmerge (s_hi, fp_posZero, n_large, vlen);
  s_lo = __riscv_vmerge (s_lo, s, n_large, vlen);

  /* For sinh: A = 1 - s_hi */
  VFLOAT A = __riscv_vfrsub (s_hi, fp_posOne, vlen);
  s_lo = __riscv_vfsgnjn (s_lo, s_lo, vlen);
  b = __riscv_vfadd (b, s_lo, vlen);

  VFLOAT Z_hi, Z_lo;
  FAST2SUM (B, poly, Z_hi, Z_lo, vlen);
  b = __riscv_vfadd (b, Z_lo, vlen);
  B = Z_hi;
  FAST2SUM (A, B, Z_hi, Z_lo, vlen);
  b = __riscv_vfadd (b, Z_lo, vlen);
  vy = __riscv_vfadd (Z_hi, b, vlen);

  /* scale vy by 2^(n-1) */
  n = __riscv_vsub (n, 1, vlen);
  FAST_LDEXP (vy, n, vlen);

  /* Apply sign for sinh */
  vy = __riscv_vfsgnj (vy, vx_orig, vlen);

  vy = __riscv_vmerge (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
