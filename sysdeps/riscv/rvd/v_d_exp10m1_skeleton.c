/* Double-precision vector exp10m1 function.
   exp10m1(x) computes 10^x - 1

   Ported from glibc aarch64 Advanced SIMD implementation.

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

static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (exp10m1_expm1_kernel, LMUL) (VFLOAT vx, size_t vlen)
{
  /* Argument reduction.  */
  VFLOAT flt_n = __riscv_vfmul (vx, RVV_LN2_INV_D, vlen);
  VINT n = __riscv_vfcvt_x (flt_n, vlen);
  flt_n = __riscv_vfcvt_f (n, vlen);

  VFLOAT r_tmp = __riscv_vfnmsac (vx, RVV_LN2_HI_D, flt_n, vlen);
  VFLOAT r = __riscv_vfnmsub (flt_n, RVV_LN2_LO_D, r_tmp, vlen);
  VFLOAT r_lo = __riscv_vfsub (r_tmp, r, vlen);
  r_lo = __riscv_vfnmsac (r_lo, RVV_LN2_LO_D, flt_n, vlen);
  /* r is the reduced argument in working precision; but r + r_lo is extra
     precise exp(r+r_lo) is 1 + (r+r_lo) + (r+r_lo)^2/2 + r^3 * polynomial(r)
     1 + r + r^2/2 + (r_lo + r * r_lo) + r^3 * polynomial(r).  */
  r_lo = __riscv_vfmacc (r_lo, r, r_lo, vlen);

  /* Compute P_head + P_tail as r + r^2/2 accurately.  */
  VFLOAT r_prime = __riscv_vfmul (r, RVV_HALF_D, vlen); /* 1/2 */
  VFLOAT P_head = __riscv_vfmadd (r, r_prime, r, vlen);
  VFLOAT P_tail = __riscv_vfsub (r, P_head, vlen);
  P_tail = __riscv_vfmacc (P_tail, r, r_prime, vlen);

  /* Polynomial computation, degree 11.  */
  VFLOAT poly_right = PSTEP (
      RVV_EXPM1_POLY_RIGHT0_D, r,
      PSTEP (RVV_EXPM1_POLY_RIGHT1_D, r,
             PSTEP (RVV_EXPM1_POLY_RIGHT2_D, RVV_EXPM1_POLY_RIGHT3_D, r, vlen),
             vlen),
      vlen);

  VFLOAT poly_mid = PSTEP (
      RVV_EXPM1_POLY_MID0_D, r,
      PSTEP (RVV_EXPM1_POLY_MID1_D, RVV_EXPM1_POLY_MID2_D, r, vlen), vlen);

  VFLOAT poly_left
      = PSTEP (RVV_EXPM1_POLY_LEFT0_D, r,
               PSTEP (RVV_EXPM1_POLY_LEFT1_D, RVV_EXPM1_POLY_LEFT2_D, r, vlen),
               vlen);

  VFLOAT r_sq = __riscv_vfmul (r, r, vlen);
  VFLOAT r_cube = __riscv_vfmul (r_sq, r, vlen);

  VFLOAT poly = __riscv_vfmadd (poly_right, r_cube, poly_mid, vlen);
  poly = __riscv_vfmadd (poly, r_cube, poly_left, vlen);

  /* Compute 1 - 2^(-n) accurately.  */
  VFLOAT One = VFMV_VF (fp_posOne, vlen);
  VINT n_clip = __riscv_vmin (n, RVV_EXP10M1_N_CLIP_MAX, vlen);
  VBOOL n_le53 = __riscv_vmsle (n_clip, RVV_EXP10M1_N_LE53, vlen);
  VBOOL n_le_neg54 = __riscv_vmsle (n_clip, RVV_EXP10M1_N_LE_NEG54, vlen);
  VINT I_tail = __riscv_vrsub (n_clip, -(EXP_BIAS + 2), vlen);
  /* u = -2^(-n_clip).  */
  VFLOAT u = I_AS_F (__riscv_vsll (I_tail, MAN_LEN, vlen));
  I_tail = __riscv_vmerge (I_tail, 0, n_le53, vlen);
  I_tail = __riscv_vmerge (I_tail, EXP_BIAS, n_le_neg54, vlen);
  VFLOAT a = I_AS_F (__riscv_vsll (I_tail, MAN_LEN, vlen));
  VFLOAT A = __riscv_vfadd (One, u, vlen);

  /* Compute A + a + P_head + P_tail + r_lo + r^3 * poly.  */
  P_tail = __riscv_vfadd (P_tail, a, vlen);
  P_tail = __riscv_vfadd (P_tail, r_lo, vlen);
  poly = __riscv_vfmadd (poly, r_cube, P_tail, vlen);
  P_head = __riscv_vfadd (P_head, poly, vlen);

  VFLOAT vy = __riscv_vfadd (P_head, A, vlen);
  /* vy is now exp(r) - 1 accurately.  */

  /* Need to compute 2^n * exp(r).  */
  FAST_LDEXP (vy, n, vlen);
  return vy;
}

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  const struct v_exp10m1_data *data = ptr_barrier (&__v_exp10m1_data);
  VFLOAT vx = x;
  VFLOAT x_orig = vx;
  VFLOAT vy, vy_special;
  VBOOL special_args;

  SET_ROUNDTONEAREST;

  /* Handle NaN/Inf.  */
  VUINT vclass = __riscv_vfclass (vx, vlen);
  IDENTIFY (vclass, class_NaN | class_Inf, special_args, vlen);
  vy_special = VFMV_VF (fp_posZero, vlen);
  UINT nb_special_args = __riscv_vcpop (special_args, vlen);
  if (nb_special_args > 0)
    {
      /* Substitute -Inf with -1.  */
      VBOOL id_mask;
      IDENTIFY (vclass, class_negInf, id_mask, vlen);
      vx = __riscv_vfmerge (vx, fp_negOne, id_mask, vlen);
      vy_special = __riscv_vfmul (special_args, vx, fp_posOne, vlen);
      vx = __riscv_vfmerge (vx, fp_posZero, special_args, vlen);
    }

  /* Clip to avoid overflow/underflow in the reduction path.  */
  vx = FCLIP (vx, RVV_EXP10_X_MIN_D, RVV_EXP10_X_MAX_D, vlen);

  VFLOAT abs_x = __riscv_vfsgnj (vx, fp_posOne, vlen);
  VBOOL tiny_input = __riscv_vmflt (abs_x, RVV_EXP10M1_TINY_BOUND_D, vlen);

  /* n = round(x/(log10(2)/N)).  */
  VFLOAT shift = VFMV_VF (RVV_EXP10M1_ROUND_SHIFT_D, vlen);
  VFLOAT z = __riscv_vfmadd (vx, RVV_EXP10M1_LOG10_2_D, shift, vlen);
  VUINT u = F_AS_U (z);
  VFLOAT n_flt = __riscv_vfsub (z, shift, vlen);

  /* r = x - n*log10(2)/N, with extra precision.  */
  VFLOAT r = __riscv_vfnmsac (vx, RVV_EXP10M1_LOG10_2_BY_N_HI_D, n_flt,
                              vlen);
  r = __riscv_vfnmsac (r, RVV_EXP10M1_LOG10_2_BY_N_LO_D, n_flt, vlen);

  /* p = 10^r - 1 = expm1(r * ln(10)).
     r is tiny due to N=128 reduction; use a short Taylor series.  */
  VFLOAT t = __riscv_vfmul (r, RVV_LN10_HI_D, vlen);
  VFLOAT t2 = __riscv_vfmul (t, t, vlen);
  VFLOAT poly = VFMV_VF (RVV_EXP10M1_OUTER_POLY0_D, vlen); /* 1/720 */
  poly = PSTEP (RVV_EXP10M1_OUTER_POLY1_D, t, poly, vlen);  /* 1/120 */
  poly = PSTEP (RVV_EXP10M1_OUTER_POLY2_D, t, poly, vlen);  /* 1/24 */
  poly = PSTEP (RVV_EXP10M1_OUTER_POLY3_D, t, poly, vlen);  /* 1/6 */
  poly = PSTEP (RVV_EXP10M1_OUTER_POLY4_D, t, poly, vlen);  /* 1/2 */
  VFLOAT p = __riscv_vfmadd (poly, t2, t, vlen);

  /* scale = 2^(n/N), using the same bit-manipulation trick as aarch64.  */
  VUINT ind = __riscv_vand (u, RVV_EXP10M1_TABLE_INDEX_MASK, vlen);
  ind = __riscv_vsll (ind, 3, vlen);
  VUINT scale_bits = __riscv_vluxei64 (data->exp_data, ind, vlen);
  VUINT tail_bits = __riscv_vluxei64 (data->exp_tail, ind, vlen);
  VUINT e = __riscv_vsll (u, RVV_EXP10M1_SCALE_SHIFT, vlen);
  VFLOAT scale = U_AS_F (__riscv_vadd (scale_bits, e, vlen));
  VFLOAT tail = U_AS_F (tail_bits);

  /* scalem1 = scale - 1, use a lookup table for small |x|.  */
  VFLOAT scalem1 = __riscv_vfsub (scale, fp_posOne, vlen);
  VBOOL is_small = __riscv_vmflt (abs_x, RVV_EXP10M1_TABLE_BOUND_D, vlen);
  if (__riscv_vcpop (is_small, vlen) > 0)
    {
      VUINT base_idx = __riscv_vand (u, RVV_EXP10M1_SMALL_INDEX_MASK, vlen);
      VBOOL x_lt_rnd2zero = __riscv_vmflt (vx, RVV_EXP10M1_ROUND2ZERO_D, vlen);
      VUINT off = __riscv_vmerge (VMVU_VX (0, vlen),
                                  RVV_EXP10M1_SMALL_INDEX_OFFSET,
                                  x_lt_rnd2zero, vlen);
      VUINT idx = __riscv_vadd (base_idx, off, vlen);
      idx = __riscv_vsll (idx, 3, vlen);
      VUINT sm1bits = __riscv_vluxei64 (data->scalem1, idx, vlen);
      VFLOAT lookup = U_AS_F (sm1bits);
      scalem1 = __riscv_vmerge (scalem1, lookup, is_small, vlen);
    }

  /* y = (scale - 1) + scale * p.  */
  vy = __riscv_vfmadd (scale, p, scalem1, vlen);
  /* Apply tail correction: for non-small x, scalem1 is computed from the head
     only and needs the missing scale_tail term; for small x, scalem1 already
     comes from a high-precision lookup table.  */
  VFLOAT scale_tail = __riscv_vfmul (scale, tail, vlen);
  vy = __riscv_vfmadd (scale_tail, p, vy, vlen);
  VFLOAT vy_tail_added = __riscv_vfadd (vy, scale_tail, vlen);
  vy = __riscv_vmerge (vy_tail_added, vy, is_small, vlen);

  /* Special path for very large |x| to avoid overflow in intermediate
     computations near the limit.  */
  VBOOL use_special = __riscv_vmfge (abs_x, RVV_EXP10M1_SPECIAL_BOUND_D, vlen);
  if (__riscv_vcpop (use_special, vlen) > 0)
    {
      VBOOL not_special = __riscv_vmnot (use_special, vlen);
      VFLOAT n_sc = VFMERGE_VFM (n_flt, fp_posZero, not_special, vlen);
      VFLOAT scale_sc = VFMERGE_VFM (scale, fp_posOne, not_special, vlen);
      VFLOAT p_sc = VFMERGE_VFM (p, fp_posZero, not_special, vlen);

      VBOOL n_neg = __riscv_vmflt (n_sc, fp_posZero, vlen);
      VUINT b = __riscv_vmerge (VMVU_VX (0, vlen),
                                RVV_EXP10M1_SPECIAL_OFFSET, n_neg, vlen);

      VUINT s1_bits = __riscv_vsub (VMVU_VX (RVV_EXP10M1_SPECIAL_BIAS1, vlen),
                                    b, vlen);
      VFLOAT s1 = U_AS_F (s1_bits);

      VUINT s2_bits = __riscv_vsub (F_AS_U (scale_sc),
                                    RVV_EXP10M1_SPECIAL_BIAS2, vlen);
      s2_bits = __riscv_vadd (s2_bits, b, vlen);
      VFLOAT s2 = U_AS_F (s2_bits);

      VFLOAT abs_n = __riscv_vfsgnj (n_sc, fp_posOne, vlen);
      VBOOL n_overflow = __riscv_vmfgt (abs_n, RVV_EXP10M1_SCALE_BOUND_D, vlen);

      VFLOAT s1_r0 = VFMERGE_VFM (s1, fp_posOne, n_overflow, vlen);
      VFLOAT s2_r0 = VFMERGE_VFM (s2, fp_posOne, n_overflow, vlen);
      VFLOAT p_r0 = VFMERGE_VFM (p_sc, fp_posZero, n_overflow, vlen);
      VFLOAT r0 = __riscv_vfmadd (p_r0, s2_r0, s2_r0, vlen);
      r0 = __riscv_vfmul (r0, s1_r0, vlen);

      VBOOL not_n_overflow = __riscv_vmnot (n_overflow, vlen);
      VFLOAT s1_r1 = VFMERGE_VFM (s1, fp_posOne, not_n_overflow, vlen);
      VFLOAT r1 = __riscv_vfmul (s1_r1, s1_r1, vlen);
      VFLOAT r_sel = __riscv_vmerge (r0, r1, n_overflow, vlen);
      VFLOAT y_special = __riscv_vfsub (r_sel, fp_posOne, vlen);

      vy = __riscv_vmerge (vy, y_special, use_special, vlen);
    }

  /* For small |x|, use expm1(x*ln10) to avoid cancellation in 10^x - 1.  */
  VBOOL use_expm1
    = __riscv_vmfle (abs_x, RVV_EXP10M1_USE_EXPM1_BOUND_D, vlen);
  if (__riscv_vcpop (use_expm1, vlen) > 0)
    {
      /* Compute x*ln10 in extra precision.  */
      VFLOAT t_hi = __riscv_vfmul (vx, RVV_LN10_HI_D, vlen);
      VFLOAT t_err = __riscv_vfmsub (vx, RVV_LN10_HI_D, t_hi, vlen);
      VFLOAT t_lo = __riscv_vfmadd (vx, RVV_LN10_LO_D, t_err, vlen);
      VFLOAT t_expm1 = __riscv_vfadd (t_hi, t_lo, vlen);
      VFLOAT y_expm1 = RVV_LIBMVEC_IMPL_NAME (exp10m1_expm1_kernel, LMUL) (
          t_expm1, vlen);
      vy = __riscv_vmerge (vy, y_expm1, use_expm1, vlen);
    }

  if (__riscv_vcpop (tiny_input, vlen) > 0)
    {
      /* For subnormals and signed zero, exp10m1(x) rounds to x * ln(10).  */
      VFLOAT tiny_y = __riscv_vfmul (vx, RVV_LN10_HI_D, vlen);
      tiny_y = __riscv_vfmadd (vx, RVV_LN10_LO_D, tiny_y, vlen);
      vy = __riscv_vmerge (vy, tiny_y, tiny_input, vlen);
    }

  VBOOL negzero = __riscv_vmseq (F_AS_U (x_orig),
                                 VMVU_VX (RVV_EXP10M1_NEGZERO_BITS, vlen),
                                 vlen);
  vy = __riscv_vfmerge (vy, fp_negZero, negzero, vlen);

  vy = __riscv_vmerge (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
