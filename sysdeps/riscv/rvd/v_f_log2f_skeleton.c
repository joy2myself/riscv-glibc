/* Single-precision vector log2f function.
   log2f(x) computes the base-2 logarithm of x

   Based on logf implementation with log2(e) scaling.

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
  VFLOAT vx, vy, vy_special;
  VBOOL special_args;
  VINT n_adjust, e_adjust;

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
      n_adjust = VMERGE_VXM_I
        (n_adjust, RVV_LOG2F_SUBNORMAL_EXP_ADJUST, id_mask, vlen);
      VFLOAT vx_normalized = VFMUL_VF_M
        (id_mask, vx, RVV_LOG2F_SCALE_F, vlen);
      vx = VMERGE_VVM_F(vx, vx_normalized, id_mask, vlen);
    }

  /* Scale large finite inputs to avoid spurious overflow in vx * (4/3). */
  VBOOL adjust_large = __riscv_vmfge(vx, RVV_LOG2F_CLEAN_BOUND_F, vlen);
  e_adjust = VMVI_VX(0, vlen);
  e_adjust = __riscv_vmerge
    (e_adjust, RVV_LOG2F_CLEAN_LARGE_EXP_ADJUST, adjust_large, vlen);
  vx = __riscv_vfmul(
      vx, __riscv_vfmerge
        (VFMV_VF(1.0f, vlen), RVV_LOG2F_CLEAN_SCALE_DOWN_F, adjust_large,
         vlen),
      vlen);

  /* Range reduction: multiply by 4/3 */
  VFLOAT vx_adj = __riscv_vfmul(vx, RVV_LOG2F_INV_075_F, vlen);

  /* Extract e = vilogb2k((4/3) * vx), then restore the exponent removed by
     the clean large-input adapter.  */
  VINT e_core = RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_log2f_vilogb2k_m, LMUL) (vx_adj, vlen);
  VINT e = __riscv_vadd(e_core, e_adjust, vlen);

  /* Compute m = vldexp3(vx, -e_core).  */
  VINT neg_e = __riscv_vneg(e_core, vlen);
  VFLOAT m = RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_log2f_vldexp3_m, LMUL)
    (vx, neg_e, vlen);
  e = __riscv_vsub(e, n_adjust, vlen);

  /* Float-float division: retain SLEEF's scalar-first dfadd2 and FMA
     reciprocal-residual association.  The EH envelope above/below is the
     clean adaptation, not a reason to change this finite graph.  */
  VFLOAT num_hi, num_lo;
  VFLOAT den_hi, den_lo;
  RVV_SLEEF_ORIG_DFADD2_VF_VF (VFMV_VF (-1.0f, vlen), m,
                               num_hi, num_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF_VF (VFMV_VF (1.0f, vlen), m,
                               den_hi, den_lo, vlen);

  VFLOAT x_hi, x_lo;
  RVV_SLEEF_ORIG_DFDIV_VF2_VF2 (num_hi, num_lo, den_hi, den_lo,
                                x_hi, x_lo, vlen);

  /* x2 = x_hi^2, poly = C2 + x2*(C1 + x2*C0), poly_tail = x2 * x_hi * poly */
  VFLOAT x2 = __riscv_vfmul(x_hi, x_hi, vlen);

  VFLOAT poly = VFMV_VF (RVV_LOG2F_POLY0_F, vlen);
  poly = __riscv_vfmacc
    (VFMV_VF (RVV_LOG2F_POLY1_F, vlen), poly, x2, vlen);
  poly = __riscv_vfmacc
    (VFMV_VF (RVV_LOG2F_POLY2_F, vlen), poly, x2, vlen);

  VFLOAT poly_tail = __riscv_vfmul(__riscv_vfmul(x2, x_hi, vlen), poly, vlen);

  VFLOAT ef = __riscv_vfcvt_f(e, vlen);

  /* Compute log2(x) in double-float: s = e + x*(2/ln2) + poly_tail */
  VFLOAT x_scale_hi, x_scale_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (x_hi, x_lo, VFMV_VF (RVV_LOG2F_X_HI_F, vlen),
     VFMV_VF (RVV_LOG2F_X_LO_F, vlen), x_scale_hi, x_scale_lo, vlen);

  VFLOAT s_hi, s_lo;
  RVV_SLEEF_ORIG_DFADD2_VF_VF2 (ef, x_scale_hi, x_scale_lo,
                                s_hi, s_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF (s_hi, s_lo, poly_tail,
                                s_hi, s_lo, vlen);

  vy = __riscv_vfadd(s_hi, s_lo, vlen);

  vy = VMERGE_VVM_F(vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
