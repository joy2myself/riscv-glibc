/* Single-precision vector asin function.
   asin(x) is defined for |x| <= 1.

   Based on SLEEF xasinf implementation.

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

/* ============================================================================
   Shared float-float helpers live in v_math.h.  */

/* Pinned SLEEF dfsqrt_vf2_vf operation graph.  The caller supplies a safe
   argument for lanes whose selected result is the small-path |x| value.  */
static inline void
RVV_ASINF_SOURCE_MANGLE (asinf_source_dfsqrt)
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

/* Pinned SLEEF dfsub_vf2_vf2_vf2 and dfsub_vf2_vf operation graphs.  */
static inline void
RVV_ASINF_SOURCE_MANGLE (asinf_source_dfsub_pair_pair)
  (VFLOAT x_hi, VFLOAT x_lo, VFLOAT y_hi, VFLOAT y_lo,
   VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  RVV_SLEEF_ORIG_DFSUB_VF2_VF2
    (x_hi, x_lo, y_hi, y_lo, *r_hi, *r_lo, vlen);
}

static inline void
RVV_ASINF_SOURCE_MANGLE (asinf_source_dfsub_pair_scalar)
  (VFLOAT x_hi, VFLOAT x_lo, VFLOAT y,
   VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT s = VFSUB_VV (x_hi, y, vlen);
  VFLOAT t = VFSUB_VV (x_hi, s, vlen);
  t = VFSUB_VV (t, y, vlen);
  *r_hi = s;
  *r_lo = VFADD_VV (t, x_lo, vlen);
}

/* asin(x) for |x| <= 1:
   For |x| < 0.5: asin(x) ~= x + x^3 * poly(x^2)
   For |x| >= 0.5: asin(x) = pi/2 - 2*asin(sqrt((1-|x|)/2))
   This matches SLEEF xasinf_u1 with double-float subtraction. */

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT vx, vy, vy_special;
  VBOOL x_nan, x_subnormal, special_args;

  SET_ROUNDTONEAREST;
  vx = x;

  VUINT vclass = __riscv_vfclass (vx, vlen);
  IDENTIFY_F (vclass, class_NaN, x_nan, vlen);
  IDENTIFY_F (vclass, class_Subnormal, x_subnormal, vlen);

  /* Compute |x| and select the reduction path. */
  VFLOAT abs_x = __riscv_vfabs (vx, vlen);
  VFLOAT abs_x_cmp = __riscv_vfmerge (abs_x, fp_posZero, x_nan, vlen);
  VBOOL abs_x_gt1 = __riscv_vmfgt (abs_x_cmp, fp_posOne, vlen);
  special_args = __riscv_vmor (abs_x_gt1, x_nan, vlen);
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

      vx = __riscv_vfmerge (vx, fp_posZero, special_args, vlen);
      abs_x = __riscv_vfabs (vx, vlen);
    }

  /* Do not evaluate the finite core when every lane is already resolved.
     This prevents discarded sanitized lanes from changing the caller's
     exception flags.  */
  VBOOL main_args = __riscv_vmnot (special_args, vlen);
  if (__riscv_vcpop (main_args, vlen) == 0)
    {
      RESTORE_FRM;
      return vy_special;
    }

  VBOOL x_lt_half = __riscv_vmflt (abs_x, 0.5f, vlen);

  /* x2 = x*x for small |x|; otherwise x2 = (1 - |x|) / 2.  The scalar
     asin path forms the small square in a wider domain.  Quarantine lanes
     whose binary32 square would underflow so the clean vector path does not
     add FE_UNDERFLOW that scalar RNE does not report; retain vx for the
     source-shaped small result and sign restoration.  */
  VBOOL tiny_square = __riscv_vmflt (abs_x, 0x1.0p-63f, vlen);
  VFLOAT square_arg = __riscv_vfmerge (vx, fp_posZero, tiny_square, vlen);
  VFLOAT x2_small = __riscv_vfmul (square_arg, square_arg, vlen);
  /* Subnormal scalar asin inputs do report the underflow from their square.
     Re-issue that event only for subnormal lanes, while the benign zero
     operand keeps normal tiny lanes from inheriting it.  The result is zero
     and is merged solely to keep the operation observable.  */
  VBOOL not_subnormal = __riscv_vmnot (x_subnormal, vlen);
  VFLOAT subnormal_square_arg
    = __riscv_vfmerge (vx, fp_posZero, not_subnormal, vlen);
  VFLOAT subnormal_square
    = __riscv_vfmul (subnormal_square_arg, subnormal_square_arg, vlen);
  x2_small = __riscv_vmerge (x2_small, subnormal_square, x_subnormal, vlen);
  VFLOAT one_minus_abs = __riscv_vfsub (VFMV_VF (1.0f, vlen), abs_x, vlen);
  VFLOAT x2_large = __riscv_vfmul (one_minus_abs, 0.5f, vlen);
  VFLOAT x2 = __riscv_vmerge (x2_large, x2_small, x_lt_half, vlen);

  /* Endpoint lanes also discard the sqrt result; x2 is exactly zero there.
     Quarantine them along with the small branch so the source-shaped
     reciprocal correction cannot raise FE_DIVBYZERO.  */
  VBOOL abs_x_eq1 = __riscv_vmfeq (abs_x, 1.0f, vlen);
  VBOOL sqrt_quarantine = __riscv_vmor (x_lt_half, abs_x_eq1, vlen);

  /* x = |x| (small) or sqrt(x2) (large), in double-float. */
  VFLOAT sqrt_hi, sqrt_lo;
  VFLOAT sqrt_arg = __riscv_vfmerge (x2, fp_posOne, sqrt_quarantine, vlen);
  RVV_ASINF_SOURCE_MANGLE (asinf_source_dfsqrt)
    (sqrt_arg, &sqrt_hi, &sqrt_lo, vlen);
  VFLOAT x_hi = __riscv_vmerge (sqrt_hi, abs_x, x_lt_half, vlen);
  VFLOAT x_lo = __riscv_vmerge (sqrt_lo, VFMV_VF (0.0f, vlen), x_lt_half,
                                vlen);

  /* For |x| == 1, force x = 0 to avoid cancellation. */
  x_hi = __riscv_vmerge (x_hi, VFMV_VF (0.0f, vlen), abs_x_eq1, vlen);
  x_lo = __riscv_vmerge (x_lo, VFMV_VF (0.0f, vlen), abs_x_eq1, vlen);

  /* Polynomial: u = poly(x2) * x2 * x. */
  VFLOAT u = VFMV_VF (RVV_ASINF_C4, vlen);
  u = __riscv_vfmadd (u, x2, VFMV_VF (RVV_ASINF_C3, vlen), vlen);
  u = __riscv_vfmadd (u, x2, VFMV_VF (RVV_ASINF_C2, vlen), vlen);
  u = __riscv_vfmadd (u, x2, VFMV_VF (RVV_ASINF_C1, vlen), vlen);
  u = __riscv_vfmadd (u, x2, VFMV_VF (RVV_ASINF_C0, vlen), vlen);
  u = __riscv_vfmul (u, __riscv_vfmul (x2, x_hi, vlen), vlen);

  /* y = pi/4 - x - u in double-float. */
  VFLOAT pi4_hi = VFMV_VF (RVV_ASINF_PIBY4_HI, vlen);
  VFLOAT pi4_lo = VFMV_VF (RVV_ASINF_PIBY4_LO, vlen);
  VFLOAT y1_hi, y1_lo;
  RVV_ASINF_SOURCE_MANGLE (asinf_source_dfsub_pair_pair)
    (pi4_hi, pi4_lo, x_hi, x_lo, &y1_hi, &y1_lo, vlen);
  VFLOAT y_hi, y_lo;
  RVV_ASINF_SOURCE_MANGLE (asinf_source_dfsub_pair_scalar)
    (y1_hi, y1_lo, u, &y_hi, &y_lo, vlen);

  /* Select the reconstruction for small/large ranges. */
  VFLOAT small = __riscv_vfadd (u, x_hi, vlen);
  VFLOAT large = __riscv_vfmul (__riscv_vfadd (y_hi, y_lo, vlen), 2.0f, vlen);
  vy = __riscv_vmerge (large, small, x_lt_half, vlen);

  /* Restore sign. */
  vy = __riscv_vfsgnjx (vy, vx, vlen);
  if (nb_special_args > 0)
    vy = __riscv_vmerge (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
