/* Double-precision vector tanh function.
   tanh(x) computes the hyperbolic tangent of x

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

/* This version uses tanh(x) = (1 - exp(-2x)) / (1 + exp(-2x)) with
   extra-precision arithmetic.  */

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT vx_orig = x;
  VFLOAT vx, vy, vy_special;
  VBOOL special_args;
  VUINT expo_x;

  SET_ROUNDTONEAREST;

  expo_x = __riscv_vand (__riscv_vsrl (F_AS_U (vx_orig), MAN_LEN, vlen), 0x7FF,
                         vlen);

  /* Handle special cases: NaN, Inf, and very small |x|.  */
  VUINT vclass = __riscv_vfclass (vx_orig, vlen);
  VBOOL NaN_Inf;
  IDENTIFY (vclass, class_NaN | class_Inf, NaN_Inf, vlen);
  VBOOL small_x = __riscv_vmsleu (expo_x, EXP_BIAS - MAN_LEN - 5, vlen);
  special_args = __riscv_vmor (NaN_Inf, small_x, vlen);

  vy_special = VFMV_VF (fp_posZero, vlen);
  UINT nb_special_args = __riscv_vcpop (special_args, vlen);
  if (nb_special_args > 0)
    {
      /* For small x: tanh(x) ~= x + x^3/3 ~= x * (1 + x^2 * small_constant).  */
      vy_special
        = __riscv_vfmadd (special_args, vx_orig, 0x1.0p-60, vx_orig, vlen);
      /* For Inf: tanh(Inf) = 1, tanh(-Inf) = -1.  */
      VBOOL Inf;
      IDENTIFY (vclass, class_Inf, Inf, vlen);
      VFLOAT one = VFMV_VF (fp_posOne, vlen);
      one = __riscv_vfsgnj (one, vx_orig, vlen);
      vy_special = __riscv_vmerge (vy_special, one, Inf, vlen);
      vx_orig = __riscv_vfmerge (vx_orig, fp_posZero, special_args, vlen);
    }

  /* tanh(x) = sign(x) * tanh(|x|); work on |x| for the main part.  */
  vx = __riscv_vfsgnj (vx_orig, fp_posOne, vlen);

  /* Clip |x| to 20, which is bigger than 28 log(2).  */
  vx = __riscv_vfmin (vx, 0x1.4p4, vlen);

  /* Compute exp(-2x) by replacing x by -2x.  */
  vx = __riscv_vfmul (vx, -0x1.0p1, vlen);

  /* Argument reduction: exp(vx) = 2^n * exp(r), where |r| <= log(2)/2.  */
  VINT n;
  VFLOAT r, r_delta;
  VFLOAT n_flt = __riscv_vfmul (vx, RVV_LN2_INV_D, vlen);
  n = __riscv_vfcvt_x (n_flt, vlen);
  n_flt = __riscv_vfcvt_f (n, vlen);
  r_delta = __riscv_vfnmsac (vx, RVV_LN2_HI_D, n_flt, vlen);
  r = __riscv_vfnmsac (r_delta, RVV_LN2_LO_D, n_flt, vlen);
  r_delta = __riscv_vfsub (r_delta, r, vlen);
  r_delta = __riscv_vfnmsac (r_delta, RVV_LN2_LO_D, n_flt, vlen);

  /* Compute exp(r) using polynomial approximation.  */
  VFLOAT rsq = __riscv_vfmul (r, r, vlen);
  VFLOAT rcube = __riscv_vfmul (rsq, r, vlen);

  VFLOAT p_even
    = PSTEP (RVV_COSH_P_EVEN0_D, rsq,
             PSTEP (RVV_COSH_P_EVEN1_D, rsq,
                    PSTEP (RVV_COSH_P_EVEN2_D, rsq,
                           PSTEP (RVV_COSH_P_EVEN3_D, RVV_COSH_P_EVEN4_D,
                                  rsq, vlen),
                           vlen),
                    vlen),
             vlen);

  VFLOAT p_odd
    = PSTEP (RVV_COSH_P_ODD0_D, rsq,
             PSTEP (RVV_COSH_P_ODD1_D, rsq,
                    PSTEP (RVV_COSH_P_ODD2_D, RVV_COSH_P_ODD3_D, rsq, vlen),
                    vlen),
             vlen);

  VFLOAT poly = __riscv_vfmadd (p_odd, r, p_even, vlen);

  /* Compute r' + (r')^2/2 extra precisely.  */
  VFLOAT r_prime = __riscv_vfmul (r, 0x1.0p-1, vlen);
  VFLOAT B = __riscv_vfmadd (r, r_prime, r, vlen);
  VFLOAT b = __riscv_vfsub (r, B, vlen);
  b = __riscv_vfmacc (b, r, r_prime, vlen);
  /* Incorporate r_delta in r' + (r')^2/2.  */
  VFLOAT c = __riscv_vfmadd (r, r_delta, r_delta, vlen);
  b = __riscv_vfadd (b, c, vlen);
  poly = __riscv_vfmadd (poly, rcube, b, vlen);
  /* B + poly is r' + (r')^2/2 + r^3(...) and exp(r') ~= s*(1 + B + poly).  */

  /* Compute s = 2^n.  */
  VFLOAT s = I_AS_F (__riscv_vsll (__riscv_vadd (n, EXP_BIAS, vlen), MAN_LEN,
                                   vlen));
  VBOOL s_is_small = __riscv_vmsle (n, -(MAN_LEN + 1), vlen);
  VBOOL s_not_small = __riscv_vmnot (s_is_small, vlen);
  VFLOAT s_head = __riscv_vfmerge (s, fp_posZero, s_is_small, vlen);
  VFLOAT s_tail = __riscv_vfmerge (s, fp_posZero, s_not_small, vlen);

  /* Compute denominator 1 + exp(r) = 1 + s*(1 + B + poly).  */
  VFLOAT Z = __riscv_vfadd (s_head, fp_posOne, vlen);
  VFLOAT D_tmp = __riscv_vfmadd (B, s, Z, vlen);
  VFLOAT d_tmp = __riscv_vfsub (Z, D_tmp, vlen);
  d_tmp = __riscv_vfmacc (d_tmp, s, B, vlen);
  d_tmp = __riscv_vfadd (d_tmp, s_tail, vlen);
  d_tmp = __riscv_vfmacc (d_tmp, s, poly, vlen);
  /* Normalize D_tmp + d_tmp to D + d.  */
  VFLOAT D, d;
  FAST2SUM (D_tmp, d_tmp, D, d, vlen);

  /* Compute 1/(D+d) as E + e using extra-precision division.  */
  VFLOAT One = VFMV_VF (fp_posOne, vlen);
  VFLOAT E, e;
  DIV_N1D2 (One, D, d, E, e, vlen);

  /* Compute numerator 1 - s*(1 + B + poly).  */
  Z = __riscv_vfrsub (s_head, fp_posOne, vlen);
  VFLOAT Numer = __riscv_vfnmsub (B, s, Z, vlen);
  VFLOAT numer = __riscv_vfsub (Z, Numer, vlen);
  numer = __riscv_vfnmsac (numer, s, B, vlen);
  numer = __riscv_vfsub (numer, s_tail, vlen);
  numer = __riscv_vfnmsac (numer, s, poly, vlen);

  /* Compute (Numer + numer) * (E + e) with extra precision.  */
  vy = __riscv_vfmul (e, numer, vlen);
  vy = __riscv_vfmacc (vy, Numer, e, vlen);
  vy = __riscv_vfmacc (vy, numer, E, vlen);
  vy = __riscv_vfmacc (vy, Numer, E, vlen);

  /* Apply sign from original input.  */
  vy = __riscv_vfsgnj (vy, vx_orig, vlen);

  /* Merge special case results.  */
  vy = __riscv_vmerge (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
