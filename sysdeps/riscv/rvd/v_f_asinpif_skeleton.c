/* Single-precision vector asinpif function.
   asinpif(x) = asinf(x) / pi, defined for |x| <= 1.

   Custom implementation based on asinf algorithm with 1/pi scaling.

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

/* The finite kernel below follows the pinned SLEEF xasinf_u1 operation
   sequence.  These helpers are source-shaped; the public wrapper still
   supplies the clean EH quarantine around them.  */

static inline void
RVV_ASINPIF_SOURCE_MANGLE (asinpif_source_dfsqrt)
  (VFLOAT d, VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT t = __riscv_vfsqrt (d, vlen);

  VFLOAT tsq_hi, tsq_lo;
  RVV_SLEEF_ORIG_DFMUL_VF_VF (t, t, tsq_hi, tsq_lo, vlen);

  VFLOAT sum_hi, sum_lo;
  RVV_SLEEF_ORIG_DFADD2_VF_VF2
    (d, tsq_hi, tsq_lo, sum_hi, sum_lo, vlen);

  VFLOAT rec_hi, rec_lo;
  RVV_SLEEF_ORIG_DFREC_VF2_VF (t, rec_hi, rec_lo, vlen);

  VFLOAT q_hi, q_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (sum_hi, sum_lo, rec_hi, rec_lo, q_hi, q_lo, vlen);

  RVV_SLEEF_ORIG_DFSCALE_VF2_VF
    (q_hi, q_lo, VFMV_VF (0.5f, vlen), *r_hi, *r_lo, vlen);
}

static inline void
RVV_ASINPIF_SOURCE_MANGLE (asinpif_source_dfsub_pair_pair)
  (VFLOAT x_hi, VFLOAT x_lo, VFLOAT y_hi, VFLOAT y_lo,
   VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  RVV_SLEEF_ORIG_DFSUB_VF2_VF2
    (x_hi, x_lo, y_hi, y_lo, *r_hi, *r_lo, vlen);
}

static inline void
RVV_ASINPIF_SOURCE_MANGLE (asinpif_source_dfsub_pair_scalar)
  (VFLOAT x_hi, VFLOAT x_lo, VFLOAT y,
   VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT s = VFSUB_VV (x_hi, y, vlen);
  VFLOAT t = VFSUB_VV (x_hi, s, vlen);
  t = VFSUB_VV (t, y, vlen);
  *r_hi = s;
  *r_lo = VFADD_VV (t, x_lo, vlen);
}

/* asinpif(x) = asinf(x) / pi for |x| <= 1:
   For |x| < 0.5: asinf(x) ~= x + x^3 * poly(x^2)
   For |x| >= 0.5: asinf(x) = pi/2 - 2*asin(sqrt((1-|x|)/2))  */

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT vx, vy, vy_special;
  VBOOL x_nan, x_neg, special_args;

  SET_ROUNDTONEAREST;
  vx = x;
  VUINT ux = F_AS_U (vx);
  VUINT vclass = __riscv_vfclass (vx, vlen);
  IDENTIFY_F (vclass, class_NaN, x_nan, vlen);
  x_neg = __riscv_vmsne (__riscv_vand (ux, 0x80000000u, vlen), 0, vlen);

  /* Compute |x| */
  VFLOAT abs_x = __riscv_vfabs (vx, vlen);
  VFLOAT abs_x_cmp = __riscv_vfmerge (abs_x, fp_posZero, x_nan, vlen);
  VBOOL abs_x_gt1 = __riscv_vmfgt (abs_x_cmp, fp_posOne, vlen);
  VBOOL abs_x_eq1 = __riscv_vmfeq (abs_x_cmp, fp_posOne, vlen);
  special_args = __riscv_vmor (abs_x_gt1, abs_x_eq1, vlen);
  special_args = __riscv_vmor (special_args, x_nan, vlen);
  UINT nb_special_args = __riscv_vcpop (special_args, vlen);

  vy_special = vx;
  if (nb_special_args > 0)
    {
      VFLOAT invalid_result = __riscv_vfmerge (abs_x_cmp, fp_sNaN, abs_x_gt1,
                                               vlen);
      invalid_result = __riscv_vfadd (abs_x_gt1, invalid_result, invalid_result,
                                      vlen);
      vy_special = __riscv_vmerge (vy_special, invalid_result, abs_x_gt1, vlen);

      VFLOAT vy_nan = __riscv_vfadd (x_nan, vx, vx, vlen);
      vy_special = __riscv_vmerge (vy_special, vy_nan, x_nan, vlen);

      VFLOAT at_one = __riscv_vmerge (VFMV_VF (RVV_ASINPIF_HALF_F, vlen),
                                      VFMV_VF (-RVV_ASINPIF_HALF_F, vlen), x_neg,
                                      vlen);
      vy_special = __riscv_vmerge (vy_special, at_one, abs_x_eq1, vlen);

      vx = __riscv_vfmerge (vx, fp_posZero, special_args, vlen);
      abs_x = __riscv_vfabs (vx, vlen);
    }

  /* Do not evaluate even the sanitized finite core when all lanes have
     already been resolved.  In a mixed vector, the sanitized lanes are
     deliberately safe inputs and therefore cannot leak discarded-lane
     exceptions.  */
  VBOOL main_args = __riscv_vmnot (special_args, vlen);
  if (__riscv_vcpop (main_args, vlen) == 0)
    {
      RESTORE_FRM;
      return vy_special;
    }

  /* Check for |x| < 0.5 */
  VBOOL x_lt_half = __riscv_vmflt (abs_x, RVV_ASINPIF_HALF_F, vlen);

  /* For |x| < 0.5: x2 = x * x
     For |x| >= 0.5: x2 = (1 - |x|) / 2 */
  /* For |x| < 2^-63, x*x is below the binary32 normal range while the
     cubic correction is far below one output ulp.  Compute that discarded
     intermediate from zero instead of allowing it to raise an underflow;
     the linear term still uses the original x through x_hi.  */
  VBOOL tiny_small
      = __riscv_vmflt (abs_x, RVV_ASINPIF_TINY_SQUARE_LIMIT_F, vlen);
  VFLOAT x_small_safe = __riscv_vfmerge (vx, fp_posZero, tiny_small, vlen);
  VFLOAT x2_small = __riscv_vfmul (x_small_safe, x_small_safe, vlen);
  VFLOAT one_minus_abs = __riscv_vfsub (VFMV_VF (1.0f, vlen), abs_x, vlen);
  VFLOAT x2_large = __riscv_vfmul (one_minus_abs, RVV_ASINPIF_HALF_F, vlen);
  VFLOAT x2 = __riscv_vmerge (x2_large, x2_small, x_lt_half, vlen);

  /* x = |x| (small) or sqrt(x2) (large), in double-float. */
  VBOOL sqrt_quarantine = __riscv_vmor (x_lt_half, abs_x_eq1, vlen);
  VFLOAT sqrt_input = __riscv_vfmerge (x2, fp_posOne, sqrt_quarantine, vlen);
  VFLOAT sqrt_hi, sqrt_lo;
  RVV_ASINPIF_SOURCE_MANGLE (asinpif_source_dfsqrt)
    (sqrt_input, &sqrt_hi, &sqrt_lo, vlen);
  VFLOAT x_hi = __riscv_vmerge (sqrt_hi, abs_x, x_lt_half, vlen);
  VFLOAT x_lo = __riscv_vmerge (sqrt_lo, VFMV_VF (0.0f, vlen), x_lt_half,
                                vlen);

  /* For |x| == 1, force x = 0 to avoid cancellation. */
  x_hi = __riscv_vmerge (x_hi, VFMV_VF (0.0f, vlen), abs_x_eq1, vlen);
  x_lo = __riscv_vmerge (x_lo, VFMV_VF (0.0f, vlen), abs_x_eq1, vlen);

  /* Polynomial: u = poly(x2) * x2 * x. */
  VFLOAT u = VFMV_VF (RVV_ASINPIF_POLY4_F, vlen);
  u = __riscv_vfmadd (u, x2, VFMV_VF (RVV_ASINPIF_POLY3_F, vlen), vlen);
  u = __riscv_vfmadd (u, x2, VFMV_VF (RVV_ASINPIF_POLY2_F, vlen), vlen);
  u = __riscv_vfmadd (u, x2, VFMV_VF (RVV_ASINPIF_POLY1_F, vlen), vlen);
  u = __riscv_vfmadd (u, x2, VFMV_VF (RVV_ASINPIF_POLY0_F, vlen), vlen);
  u = __riscv_vfmul (u, __riscv_vfmul (x2, x_hi, vlen), vlen);

  /* y = pi/4 - x - u in double-float. */
  VFLOAT t_hi, t_lo;
  RVV_ASINPIF_SOURCE_MANGLE (asinpif_source_dfsub_pair_pair)
    (VFMV_VF (RVV_ASINPIF_PI4_HI_F, vlen),
     VFMV_VF (RVV_ASINPIF_PI4_LO_F, vlen),
     x_hi, x_lo, &t_hi, &t_lo, vlen);
  VFLOAT y_hi, y_lo;
  RVV_ASINPIF_SOURCE_MANGLE (asinpif_source_dfsub_pair_scalar)
    (t_hi, t_lo, u, &y_hi, &y_lo, vlen);

  /* Select the reconstruction for small/large ranges. */
  VFLOAT small_result = __riscv_vfadd (u, x_hi, vlen);
  VFLOAT large_result = __riscv_vfmul
    (__riscv_vfadd (y_hi, y_lo, vlen), RVV_ASINPIF_TWO_F, vlen);
  VFLOAT asin_value = __riscv_vmerge (large_result, small_result,
                                      x_lt_half, vlen);

  /* SLEEF's vmulsign applies the input sign after the positive finite
     reconstruction.  Keep the original input because vx is sanitized for
     exceptional lanes.  */
  asin_value = __riscv_vfsgnjx (asin_value, x, vlen);

  /* Target-side normalized-pi composition: collapse the source asin result
     first, then apply the split 1/pi product exactly as the source-derived
     _original comparator does.  */
  VFLOAT prod_hi, prod_lo;
  RVV_SLEEF_ORIG_DFMUL_VF_VF
    (asin_value, VFMV_VF (RVV_ASINPIF_ONE_OVER_PI_HI_F, vlen), prod_hi, prod_lo,
     vlen);
  prod_lo = __riscv_vfadd
    (prod_lo, __riscv_vfmul (asin_value, RVV_ASINPIF_ONE_OVER_PI_LO_F, vlen),
     vlen);
  vy = __riscv_vfadd (prod_hi, prod_lo, vlen);

  /* Repair signed zero after the split product collapse. */
  vy = __riscv_vfsgnj (vy, x, vlen);

  vy = __riscv_vmerge (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
