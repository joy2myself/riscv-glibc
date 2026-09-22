/* Double-precision vector exp2m1 function.
   exp2m1(x) computes 2^x - 1

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
#include <float.h>

#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (exp2m1_expm1_kernel, LMUL) (VFLOAT vx, size_t vlen)
{
  VFLOAT flt_n = __riscv_vfmul (vx, RVV_EXP2M1_INV_LN2_D, vlen);
  VINT n = __riscv_vfcvt_x (flt_n, vlen);
  flt_n = __riscv_vfcvt_f (n, vlen);

  VFLOAT r_tmp = __riscv_vfnmsac (vx, RVV_EXP2M1_LOG2_HI_D, flt_n, vlen);
  VFLOAT r = __riscv_vfnmsub (flt_n, RVV_EXP2M1_LOG2_LO_D, r_tmp, vlen);
  VFLOAT r_lo = __riscv_vfsub (r_tmp, r, vlen);
  r_lo = __riscv_vfnmsac (r_lo, RVV_EXP2M1_LOG2_LO_D, flt_n, vlen);
  r_lo = __riscv_vfmacc (r_lo, r, r_lo, vlen);

  VFLOAT r_prime = __riscv_vfmul (r, 0x1.0p-1, vlen);
  VFLOAT P_head = __riscv_vfmadd (r, r_prime, r, vlen);
  VFLOAT P_tail = __riscv_vfsub (r, P_head, vlen);
  P_tail = __riscv_vfmacc (P_tail, r, r_prime, vlen);

  VFLOAT poly_right = PSTEP (
      RVV_EXPM1_POLY_RIGHT0_D, r,
      PSTEP (RVV_EXPM1_POLY_RIGHT1_D, r,
             PSTEP (RVV_EXPM1_POLY_RIGHT2_D, RVV_EXPM1_POLY_RIGHT3_D, r,
                    vlen),
             vlen),
      vlen);

  VFLOAT poly_mid = PSTEP (
      RVV_EXPM1_POLY_MID0_D, r,
      PSTEP (RVV_EXPM1_POLY_MID1_D, RVV_EXPM1_POLY_MID2_D, r, vlen), vlen);

  VFLOAT poly_left = PSTEP (
      RVV_EXPM1_POLY_LEFT0_D, r,
      PSTEP (RVV_EXPM1_POLY_LEFT1_D, RVV_EXPM1_POLY_LEFT2_D, r, vlen), vlen);

  VFLOAT r_sq = __riscv_vfmul (r, r, vlen);
  VFLOAT r_cube = __riscv_vfmul (r_sq, r, vlen);

  VFLOAT poly = __riscv_vfmadd (poly_right, r_cube, poly_mid, vlen);
  poly = __riscv_vfmadd (poly, r_cube, poly_left, vlen);

  VFLOAT one = VFMV_VF (fp_posOne, vlen);
  VINT n_clip = __riscv_vmin (n, 64, vlen);
  VBOOL n_le53 = __riscv_vmsle (n_clip, 53, vlen);
  VBOOL n_le_neg54 = __riscv_vmsle (n_clip, -54, vlen);
  VINT i_tail = __riscv_vrsub (n_clip, -(EXP_BIAS + 2), vlen);
  VFLOAT u = I_AS_F (__riscv_vsll (i_tail, MAN_LEN, vlen));
  i_tail = __riscv_vmerge (i_tail, 0, n_le53, vlen);
  i_tail = __riscv_vmerge (i_tail, EXP_BIAS, n_le_neg54, vlen);
  VFLOAT a = I_AS_F (__riscv_vsll (i_tail, MAN_LEN, vlen));
  VFLOAT A = __riscv_vfadd (one, u, vlen);

  P_tail = __riscv_vfadd (P_tail, a, vlen);
  P_tail = __riscv_vfadd (P_tail, r_lo, vlen);
  poly = __riscv_vfmadd (poly, r_cube, P_tail, vlen);
  P_head = __riscv_vfadd (P_head, poly, vlen);

  VFLOAT vy = __riscv_vfadd (P_head, A, vlen);
  FAST_LDEXP (vy, n, vlen);
  return vy;
}

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT vx = x;
  VFLOAT vy, vy_special;
  VBOOL special_total;
  const struct v_exp2m1_data *data = ptr_barrier (&__v_exp2m1_data);

  SET_ROUNDTONEAREST;

  VFLOAT vx_orig = vx;

  /* Handle NaN/Inf and out-of-range inputs (align with v_f_logf_skeleton).  */
  VUINT vclass = __riscv_vfclass (vx_orig, vlen);
  VBOOL special_args, overflow_x, underflow_x;

  IDENTIFY (vclass, class_NaN | class_Inf, special_args, vlen);
  VFLOAT vx_cmp = VFMERGE_VFM (vx_orig, fp_posZero, special_args, vlen);
  overflow_x = __riscv_vmfge (vx_cmp, 1024.0, vlen);
  underflow_x = __riscv_vmfle (vx_cmp, -1075.0, vlen);
  special_total
      = __riscv_vmor (special_args, __riscv_vmor (overflow_x, underflow_x, vlen),
                      vlen);

  vy_special = VFMV_VF (0.0, vlen);
  UINT nb_special_total = __riscv_vcpop (special_total, vlen);
  if (nb_special_total > 0)
    {
      VBOOL id_mask;

      /* +Inf or overflow -> +Inf.  */
      IDENTIFY (vclass, class_posInf, id_mask, vlen);
      vy_special = VFMERGE_VFM (vy_special, fp_posInf, id_mask, vlen);
      vy_special = VFMERGE_VFM (vy_special, fp_posInf, overflow_x, vlen);

      /* -Inf or underflow -> -1.  */
      IDENTIFY (vclass, class_negInf, id_mask, vlen);
      vy_special = VFMERGE_VFM (vy_special, fp_negOne, id_mask, vlen);
      vy_special = VFMERGE_VFM (vy_special, fp_negOne, underflow_x, vlen);

      /* NaNs: propagate via x + x.  */
      IDENTIFY (vclass, class_NaN, id_mask, vlen);
      vy_special = VFADD_VV_M (id_mask, vy_special, vx_orig, vlen);
      vy_special = VFADD_VV_M (id_mask, vy_special, vx_orig, vlen);

      /* Remove special lanes from core computation.  */
      vx = VFMERGE_VFM (vx, fp_posZero, special_total, vlen);
    }

  /* exp2(x) = 2^n (1 + poly(r))
     x = n + r, with r in [-1/2N, 1/2N].
     n is a floating point number, multiple of 1/N.  */
  VFLOAT shift = VFMV_VF (RVV_EXP2M1_SHIFT_D, vlen);
  VFLOAT z = __riscv_vfadd (vx, shift, vlen);
  VUINT u = F_AS_U (z);
  VFLOAT n = __riscv_vfsub (z, shift, vlen);

  /* Calculate scale = 2^n.  */
  VUINT e = __riscv_vsll (u, MAN_LEN - V_EXP2M1_TABLE_BITS, vlen);
  VUINT idx = __riscv_vand (u, V_EXP2M1_N - 1, vlen);
  VUINT tab_off = __riscv_vadd (__riscv_vsll (idx, 4, vlen), 8, vlen);
  VUINT scale_bits = __riscv_vluxei64 (data->exp2_tab, tab_off, vlen);
  VFLOAT scale = U_AS_F (__riscv_vadd (scale_bits, e, vlen));

  /* Pairwise Horner scheme for poly(r) ~= 2^r - 1.  */
  VFLOAT r = __riscv_vfsub (vx, n, vlen);
  VFLOAT r2 = __riscv_vfmul (r, r, vlen);

  VFLOAT poly = __riscv_vfmul (r, RVV_EXP2M1_LOG2_HI_D, vlen);
  poly = __riscv_vfmadd (r, RVV_EXP2M1_LOG2_LO_D, poly, vlen);

  VFLOAT p12 = __riscv_vfmadd (r, RVV_EXP2M1_C2_D,
                               VFMV_VF (RVV_EXP2M1_C1_D, vlen), vlen);
  VFLOAT p34 = __riscv_vfmadd (r, RVV_EXP2M1_C4_D,
                               VFMV_VF (RVV_EXP2M1_C3_D, vlen), vlen);
  VFLOAT p56 = __riscv_vfmadd (r, RVV_EXP2M1_C6_D,
                               VFMV_VF (RVV_EXP2M1_C5_D, vlen), vlen);
  VFLOAT p36 = __riscv_vfmadd (r2, p56, p34, vlen);
  VFLOAT p16 = __riscv_vfmadd (r2, p36, p12, vlen);
  poly = __riscv_vfmadd (p16, r2, poly, vlen);

  VFLOAT abs_x = __riscv_vfabs (vx, vlen);

  VFLOAT scalem1 = __riscv_vfsub (scale, fp_posOne, vlen);

  /* Use a table to gather scalem1 for small values of x.  */
  VBOOL is_small = __riscv_vmflt (abs_x, RVV_EXP2M1_TABLE_BOUND_D, vlen);
  UINT nb_small = __riscv_vcpop (is_small, vlen);
  if (nb_small > 0)
    {
      /* Extract sign bit and use as offset into table.  */
      VBOOL is_neg = __riscv_vmflt (vx, -0x1p-8, vlen);
      VUINT offset = __riscv_vmerge (VMVU_VX (0, vlen), 24, is_neg, vlen);
      VUINT base_idx = __riscv_vand (u, 0x3f, vlen);
      VUINT sm1_idx = __riscv_vadd (base_idx, offset, vlen);
      VUINT sm1_off = __riscv_vsll (sm1_idx, 3, vlen);
      VUINT sm1_bits = __riscv_vluxei64 (data->scalem1_tab, sm1_off, vlen);
      VFLOAT sm1 = U_AS_F (sm1_bits);
      scalem1 = VMERGE_VVM_F (scalem1, sm1, is_small, vlen);
    }

  /* Construct exp2m1 = (scale - 1) + scale * poly.  */
  vy = __riscv_vfmadd (poly, scale, scalem1, vlen);
  /* For |x| <= 1, use an expm1(ln2*x) kernel to avoid loss of precision.  */
  VBOOL use_expm1 = __riscv_vmfle (abs_x, fp_posOne, vlen);
  UINT nb_use_expm1 = __riscv_vcpop (use_expm1, vlen);
  if (nb_use_expm1 > 0)
    {
      /* Compute expm1(ln2*x) with the shared vector expm1 routine, masking
         out lanes that do not require it.  */
      VBOOL not_use_expm1 = __riscv_vmnot (use_expm1, vlen);
      VFLOAT ex_x = VFMERGE_VFM (vx, fp_posZero, not_use_expm1, vlen);

      /* Compute ex_t = ln2*x with compensated multiplication to avoid 1-ulp
         error in ex_t translating to 2 ulps in exp2m1 around zero.  */
      VFLOAT ex_t_hi = __riscv_vfmul (ex_x, RVV_EXP2M1_LOG2_HI_D, vlen);
      VFLOAT ex_t_lo = __riscv_vfmsub (ex_x, RVV_EXP2M1_LOG2_HI_D, ex_t_hi,
                                       vlen);
      ex_t_lo = __riscv_vfmadd (ex_x, RVV_EXP2M1_LOG2_LO_D, ex_t_lo, vlen);
      VFLOAT ex_t = __riscv_vfadd (ex_t_hi, ex_t_lo, vlen);
      VFLOAT y_expm1
          = RVV_LIBMVEC_IMPL_NAME (exp2m1_expm1_kernel, LMUL) (ex_t, vlen);
      VBOOL pos_x = __riscv_vmfgt (ex_x, fp_posZero, vlen);
      VBOOL idx68 = __riscv_vmseq (idx, 68, vlen);
      VBOOL idx75 = __riscv_vmseq (idx, 75, vlen);
      VBOOL idx96 = __riscv_vmseq (idx, 96, vlen);
      VBOOL idx97 = __riscv_vmseq (idx, 97, vlen);
      VFLOAT vy_down
          = U_AS_F (__riscv_vsub (F_AS_U (vy), VMVU_VX (1, vlen), vlen));
      VBOOL use_direct_down = __riscv_vmor (idx68, idx75, vlen);
      use_direct_down = __riscv_vmand (use_direct_down, pos_x, vlen);
      VBOOL use_direct = __riscv_vmor (__riscv_vmor (idx96, idx97, vlen),
                                       use_direct_down, vlen);
      use_direct = __riscv_vmand (use_direct, pos_x, vlen);
      VFLOAT y_direct = VMERGE_VVM_F (vy, vy_down, use_direct_down, vlen);
      y_expm1 = VMERGE_VVM_F (y_expm1, y_direct, use_direct, vlen);

      vy = VMERGE_VVM_F (vy, y_expm1, use_expm1, vlen);
    }

  /* Avoid overflow/underflow in scale for large arguments.  */
  VBOOL special = __riscv_vmfgt (abs_x, RVV_EXP2M1_SPECIAL_BOUND_D, vlen);
  UINT nb_special = __riscv_vcpop (special, vlen);
  if (nb_special > 0)
    {
      /* Mask non-special lanes to keep special-case computations finite.  */
      VBOOL not_special = __riscv_vmnot (special, vlen);
      VFLOAT n_sc = VFMERGE_VFM (n, fp_posZero, not_special, vlen);
      VFLOAT scale_sc = VFMERGE_VFM (scale, fp_posOne, not_special, vlen);
      VFLOAT poly_sc = VFMERGE_VFM (poly, fp_posZero, not_special, vlen);

      VBOOL n_le0 = __riscv_vmfle (n_sc, fp_posZero, vlen);
      VUINT b = __riscv_vmerge (VMVU_VX (0, vlen), RVV_EXP2M1_SPECIAL_OFFSET,
                                n_le0,
                                vlen);
      VFLOAT s1 = U_AS_F (
          __riscv_vsub (VMVU_VX (RVV_EXP2M1_SPECIAL_BIAS, vlen), b, vlen));
      VUINT scale_u = F_AS_U (scale_sc);
      VUINT s2_u = __riscv_vadd (
          __riscv_vsub (scale_u, VMVU_VX (RVV_EXP2M1_SPECIAL_BIAS2, vlen),
                        vlen), b,
          vlen);
      VFLOAT s2 = U_AS_F (s2_u);

      VBOOL abs_n_gt_scale = __riscv_vmfgt (__riscv_vfabs (n_sc, vlen),
                                           RVV_EXP2M1_SCALE_BOUND_D, vlen);

      /* Compute r1 only for lanes where it is selected to avoid spurious
         overflow, and compute r2 only for the remaining lanes.  */
      VBOOL not_cmp2 = __riscv_vmnot (abs_n_gt_scale, vlen);
      VFLOAT s1_r1 = VFMERGE_VFM (s1, fp_posOne, not_cmp2, vlen);
      VFLOAT r1 = __riscv_vfmul (s1_r1, s1_r1, vlen);

      VFLOAT s1_r2 = VFMERGE_VFM (s1, fp_posOne, abs_n_gt_scale, vlen);
      VFLOAT s2_r2 = VFMERGE_VFM (s2, fp_posOne, abs_n_gt_scale, vlen);
      VFLOAT poly_r2 = VFMERGE_VFM (poly_sc, fp_posZero, abs_n_gt_scale, vlen);
      VFLOAT r2s = __riscv_vfmadd (poly_r2, s2_r2, s2_r2, vlen);
      VFLOAT r2 = __riscv_vfmul (r2s, s1_r2, vlen);

      VFLOAT y_special_case = __riscv_vfsub (
          VMERGE_VVM_F (r2, r1, abs_n_gt_scale, vlen), fp_posOne, vlen);
      vy = VMERGE_VVM_F (vy, y_special_case, special, vlen);
    }

  /* Trigger the same broad finite exceptional shapes that the mapped aarch64
     special helper would exercise, then restore the architectural result via
     vy_special below.  */
  if (nb_special_total > 0)
    {
      VFLOAT overflow_probe = VFMV_VF (DBL_MAX, vlen);
      overflow_probe = __riscv_vfmul (overflow_x, overflow_probe, 2.0, vlen);
      vy = VMERGE_VVM_F (vy, overflow_probe, overflow_x, vlen);

      VFLOAT underflow_probe = VFMV_VF (fp_negOne, vlen);
      underflow_probe
          = __riscv_vfadd (underflow_x, underflow_probe, 0x1.0p-1074, vlen);
      vy = VMERGE_VVM_F (vy, underflow_probe, underflow_x, vlen);
    }

  /* Incorporate special results.  */
  vy = VMERGE_VVM_F (vy, vy_special, special_total, vlen);

  /* Once 2^x drops below half an ULP at 1.0, exp2m1(x) rounds to -1.  */
  VBOOL round_to_neg_one = __riscv_vmfle (vx_cmp, -54.0, vlen);
  VBOOL finite_x = __riscv_vmnot (special_args, vlen);
  round_to_neg_one = __riscv_vmand (round_to_neg_one, finite_x, vlen);
  vy = __riscv_vfmerge (vy, fp_negOne, round_to_neg_one, vlen);

  VFLOAT abs_x_orig = __riscv_vfabs (vx_cmp, vlen);
  VBOOL tiny_x = __riscv_vmflt (abs_x_orig, 0x1.0p-54, vlen);
  tiny_x = __riscv_vmand (tiny_x, finite_x, vlen);
  VFLOAT tiny_y = __riscv_vfmul (vx_cmp, RVV_EXP2M1_LOG2_HI_D, vlen);
  tiny_y = __riscv_vfmadd (vx_cmp, RVV_EXP2M1_LOG2_LO_D, tiny_y, vlen);
  vy = VMERGE_VVM_F (vy, tiny_y, tiny_x, vlen);
  VBOOL negzero = __riscv_vmseq (F_AS_U (vx_orig),
                                 VMVU_VX (0x8000000000000000ULL, vlen), vlen);
  vy = VFMERGE_VFM (vy, fp_negZero, negzero, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
