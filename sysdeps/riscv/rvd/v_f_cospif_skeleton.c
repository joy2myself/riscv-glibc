/* Single-precision vector cospi function.
   cospif(x) = cos(pi * x)

   The normal path follows SLEEF xcospif_u05/cospifk trig-pi reduction and
   polynomial reconstruction.  The special-case handling is adapted for glibc
   libmvec.

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

static inline void
RVV_COSPIF_SOURCE_DFMUL_PAIR_PAIR
  (VFLOAT a_hi, VFLOAT a_lo, VFLOAT b_hi, VFLOAT b_lo,
   VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT s = VFMUL_VV (a_hi, b_hi, vlen);
  VFLOAT residual = VFMSUB_VV (a_hi, b_hi, s, vlen);
  VFLOAT cross = VFMADD_VV (a_lo, b_hi, residual, vlen);
  *r_hi = s;
  *r_lo = VFMADD_VV (a_hi, b_lo, cross, vlen);
}

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT vx_orig = x;
  VFLOAT vx, vy, vy_special;
  VBOOL special_args;

  SET_ROUNDTONEAREST;

  vx = vx_orig;

  /* Special cases: NaN/Inf, tiny |x| (cospif rounds to 1), and |x| >= 2^24
     (float can represent only even integers, so cospif(x) = 1).  */
  VUINT vclass = __riscv_vfclass (vx, vlen);
  VBOOL NaN_Inf;
  IDENTIFY_F (vclass, class_NaN | class_Inf, NaN_Inf, vlen);

  VFLOAT ax = __riscv_vfabs (vx, vlen);
  VUINT expo_x = __riscv_vsrl (F_AS_U_F (ax), MAN_LEN_F, vlen);
  expo_x = __riscv_vand (expo_x, 0xff, vlen);
  VBOOL tiny_x = __riscv_vmsleu (expo_x, RVV_COSPIF_TINY_EXPONENT, vlen);
  VBOOL x_large = __riscv_vmsgeu (expo_x, RVV_COSPIF_LARGE_EXPONENT, vlen);
  special_args = __riscv_vmor (NaN_Inf, __riscv_vmor (tiny_x, x_large, vlen),
                               vlen);

  UINT nb_special_args = __riscv_vcpop (special_args, vlen);
  if (nb_special_args > 0)
    {
      /* Default special-case result: 1.0.  */
      vy_special = VFMV_VF_F (fp_posOne_f, vlen);

      /* For tiny |x|, use 1 + x so nonzero x raises inexact.  */
      VFLOAT vy_tiny = __riscv_vfadd (tiny_x, vx_orig, fp_posOne_f, vlen);
      vy_special = VMERGE_VVM_F (vy_special, vy_tiny, tiny_x, vlen);

      /* Substitute Inf with sNaN and propagate NaNs via x+x.  */
      VBOOL inf_mask;
      IDENTIFY_F (vclass, class_Inf, inf_mask, vlen);
      VFLOAT vy_NaN_Inf = __riscv_vfmerge (vx_orig, fp_sNaN_f, inf_mask, vlen);
      vy_NaN_Inf = __riscv_vfadd (NaN_Inf, vy_NaN_Inf, vy_NaN_Inf, vlen);
      vy_special = VMERGE_VVM_F (vy_special, vy_NaN_Inf, NaN_Inf, vlen);

      /* Mask special lanes out from the main path to avoid spurious
         exceptions (e.g. NaN->int conversion).  */
      vx = VFMERGE_VFM (vx, fp_posZero_f, special_args, vlen);
      ax = VFMERGE_VFM (ax, fp_posZero_f, special_args, vlen);
    }

  /* SLEEF cospifk reduction: u = 4*x, q = even truncation of u, s = u - q.
     The current EH/special-lane handling above is preserved.  */
  VFLOAT u4 = __riscv_vfmul (ax, 4.0f, vlen);
  VINT q = __riscv_vfcvt_rtz_x (u4, vlen);
  q = __riscv_vand (__riscv_vadd (q, 1, vlen), ~1, vlen);
  VFLOAT q_flt = __riscv_vfcvt_f (q, vlen);
  VFLOAT t_arg = __riscv_vfsub (u4, q_flt, vlen);
  VFLOAT s = __riscv_vfmul (t_arg, t_arg, vlen);
  VFLOAT s2_hi, s2_lo;
  DFMUL_VF2_VF_VF (t_arg, t_arg, s2_hi, s2_lo, vlen);

  VBOOL use_cos_poly
      = __riscv_vmseq (__riscv_vand (I_AS_U_F (q), 2, vlen), 0, vlen);

  VFLOAT p = VFMV_VF_F (RVV_COSPIF_POLY_COS0_F, vlen);
  p = VFMERGE_VFM (p, RVV_COSPIF_POLY_SIN0_F, use_cos_poly, vlen);
  VFLOAT c = VFMV_VF_F (RVV_COSPIF_POLY_COS1_F, vlen);
  c = VFMERGE_VFM (c, RVV_COSPIF_POLY_SIN1_F, use_cos_poly, vlen);
  p = __riscv_vfmadd (p, s, c, vlen);
  c = VFMV_VF_F (RVV_COSPIF_POLY_COS2_F, vlen);
  c = VFMERGE_VFM (c, RVV_COSPIF_POLY_SIN2_F, use_cos_poly, vlen);
  p = __riscv_vfmadd (p, s, c, vlen);

  VFLOAT x_hi, x_lo;
  VFLOAT c_hi = VFMV_VF_F (RVV_COSPIF_DF_FIRST_COS_HI_F, vlen);
  VFLOAT c_lo = VFMV_VF_F (RVV_COSPIF_DF_FIRST_COS_LO_F, vlen);
  c_hi = VFMERGE_VFM
    (c_hi, RVV_COSPIF_DF_FIRST_SIN_HI_F, use_cos_poly, vlen);
  c_lo = VFMERGE_VFM
    (c_lo, RVV_COSPIF_DF_FIRST_SIN_LO_F, use_cos_poly, vlen);
  DFADD2_VF2_VF_VF2 (__riscv_vfmul (p, s, vlen), c_hi, c_lo, x_hi, x_lo,
                     vlen);

  VFLOAT tmp_hi, tmp_lo;
  RVV_COSPIF_SOURCE_DFMUL_PAIR_PAIR
    (s2_hi, s2_lo, x_hi, x_lo, &tmp_hi, &tmp_lo, vlen);
  c_hi = VFMV_VF_F (RVV_COSPIF_DF_SECOND_COS_HI_F, vlen);
  c_lo = VFMV_VF_F (RVV_COSPIF_DF_SECOND_COS_LO_F, vlen);
  c_hi = VFMERGE_VFM
    (c_hi, RVV_COSPIF_DF_SECOND_SIN_HI_F, use_cos_poly, vlen);
  c_lo = VFMERGE_VFM
    (c_lo, RVV_COSPIF_DF_SECOND_SIN_LO_F, use_cos_poly, vlen);
  DFADD2_VF2_VF2_VF2 (tmp_hi, tmp_lo, c_hi, c_lo, x_hi, x_lo, vlen);

  VFLOAT mul_hi = t_arg;
  VFLOAT mul_lo = VFMV_VF_F (fp_posZero_f, vlen);
  mul_hi = VMERGE_VVM_F (mul_hi, s2_hi, use_cos_poly, vlen);
  mul_lo = VMERGE_VVM_F (mul_lo, s2_lo, use_cos_poly, vlen);
  RVV_COSPIF_SOURCE_DFMUL_PAIR_PAIR
    (x_hi, x_lo, mul_hi, mul_lo, &tmp_hi, &tmp_lo, vlen);
  x_hi = tmp_hi;
  x_lo = tmp_lo;

  VFLOAT add_hi, add_lo;
  DFADD2_VF2_VF (x_hi, x_lo, VFMV_VF_F (fp_posOne_f, vlen), add_hi, add_lo,
                 vlen);
  x_hi = VMERGE_VVM_F (x_hi, add_hi, use_cos_poly, vlen);
  x_lo = VMERGE_VVM_F (x_lo, add_lo, use_cos_poly, vlen);
  vy = __riscv_vfadd (x_hi, x_lo, vlen);

  /* SLEEF cospifk sign: flip when ((q + 2) & 4) is set.  */
  VUINT sign_bits
      = __riscv_vsll (__riscv_vand (I_AS_U_F (__riscv_vadd (q, 2, vlen)), 4,
                                    vlen),
                      BIT_WIDTH_F - 3, vlen);
  vy = __riscv_vfsgnjx (vy, U_AS_F_F (sign_bits), vlen);

  /* Preserve exact half-integer zeros expected by the glibc tests.  */
  VBOOL exact_zero = __riscv_vmand (
      __riscv_vmseq (F_AS_U_F (t_arg), 0, vlen),
      __riscv_vmnot (use_cos_poly, vlen), vlen);
  vy = __riscv_vmerge (vy, VFMV_VF_F (fp_posZero_f, vlen), exact_zero, vlen);

  if (nb_special_args > 0)
    vy = VMERGE_VVM_F (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
