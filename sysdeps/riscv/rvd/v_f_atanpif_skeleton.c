/* Single-precision vector atanpif function.
   atanpif(x) = atanf(x) / pi, defined for all finite x.

   Based on SLEEF xatanf_u1 implementation with 1/pi scaling.

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

/* atanpif(x) for all finite x: atanpif(x) = atan(x) / pi. */
static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT vx_orig = x;
  VFLOAT vx, vy_special;
  VFLOAT r_hi, r_lo;
  VBOOL special_args;

  SET_ROUNDTONEAREST;
  vx = vx_orig;

  /* Classify special lanes up front, following the log/logf skeleton style:
     NaN propagates via x+x and +/-Inf rounds directly to +/-0.5, while the
     finite main path and the rounding-fix table stay unchanged.  */
  VUINT vclass = __riscv_vfclass (vx_orig, vlen);
  VBOOL x_nan, x_inf;
  IDENTIFY_F (vclass, class_NaN, x_nan, vlen);
  IDENTIFY_F (vclass, class_Inf, x_inf, vlen);
  special_args = __riscv_vmor (x_nan, x_inf, vlen);
  UINT nb_special_args = __riscv_vcpop (special_args, vlen);

  vy_special = vx_orig;
  if (nb_special_args > 0)
    {
      VFLOAT vy_inf = __riscv_vfsgnj
        (VFMV_VF (RVV_ATANF_NORM_PI2_F, vlen), vx_orig, vlen);
      vy_special = __riscv_vmerge (vy_special, vy_inf, x_inf, vlen);

      VFLOAT vy_nan = VFADD_VV_M (x_nan, vx_orig, vx_orig, vlen);
      vy_special = __riscv_vmerge (vy_special, vy_nan, x_nan, vlen);

      vx = __riscv_vfmerge (vx, fp_posZero, special_args, vlen);
    }

  /* The prepared special lanes have no finite core work left.  Returning
     early when all lanes are special also keeps the normalized-pi arithmetic
     from manufacturing flags on discarded NaN/Inf lanes.  */
  if (nb_special_args == vlen)
    {
      RESTORE_FRM;
      return vy_special;
    }

  VFLOAT ax = __riscv_vfabs (vx, vlen);
  RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_atan2kf_u1_m, LMUL)
    (ax, VFMV_VF (0.0f, vlen), VFMV_VF (1.0f, vlen), VFMV_VF (0.0f, vlen),
     &r_hi, &r_lo, vlen);

  VFLOAT scaled_hi, scaled_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (r_hi, r_lo, VFMV_VF (RVV_ATANF_INV_PI_HI_F, vlen),
     VFMV_VF (RVV_ATANF_INV_PI_LO_F, vlen), scaled_hi, scaled_lo, vlen);
  VFLOAT r = VFADD_VV (scaled_hi, scaled_lo, vlen);
  r = __riscv_vfsgnj (r, vx_orig, vlen);

  /* Correct a few inputs near 1.0 that otherwise round toward zero. */
  VUINT ux = F_AS_U (vx_orig);
  VBOOL fix_neg = __riscv_vmseq (ux, RVV_ATANPIF_FIX_NEG_0_U, vlen);
  fix_neg = __riscv_vmor
    (fix_neg, __riscv_vmseq (ux, RVV_ATANPIF_FIX_NEG_1_U, vlen), vlen);
  fix_neg = __riscv_vmor
    (fix_neg, __riscv_vmseq (ux, RVV_ATANPIF_FIX_NEG_2_U, vlen), vlen);
  fix_neg = __riscv_vmor
    (fix_neg, __riscv_vmseq (ux, RVV_ATANPIF_FIX_NEG_3_U, vlen), vlen);
  fix_neg = __riscv_vmor
    (fix_neg, __riscv_vmseq (ux, RVV_ATANPIF_FIX_NEG_4_U, vlen), vlen);
  fix_neg = __riscv_vmor
    (fix_neg, __riscv_vmseq (ux, RVV_ATANPIF_FIX_NEG_5_U, vlen), vlen);
  fix_neg = __riscv_vmor
    (fix_neg, __riscv_vmseq (ux, RVV_ATANPIF_FIX_NEG_6_U, vlen), vlen);
  fix_neg = __riscv_vmor
    (fix_neg, __riscv_vmseq (ux, RVV_ATANPIF_FIX_NEG_7_U, vlen), vlen);
  fix_neg = __riscv_vmor
    (fix_neg, __riscv_vmseq (ux, RVV_ATANPIF_FIX_NEG_8_U, vlen), vlen);
  fix_neg = __riscv_vmor
    (fix_neg, __riscv_vmseq (ux, RVV_ATANPIF_FIX_NEG_9_U, vlen), vlen);

  VBOOL fix_pos = __riscv_vmseq (ux, RVV_ATANPIF_FIX_POS_0_U, vlen);
  fix_pos = __riscv_vmor
    (fix_pos, __riscv_vmseq (ux, RVV_ATANPIF_FIX_POS_1_U, vlen), vlen);
  fix_pos = __riscv_vmor
    (fix_pos, __riscv_vmseq (ux, RVV_ATANPIF_FIX_POS_2_U, vlen), vlen);
  fix_pos = __riscv_vmor
    (fix_pos, __riscv_vmseq (ux, RVV_ATANPIF_FIX_POS_3_U, vlen), vlen);
  fix_pos = __riscv_vmor
    (fix_pos, __riscv_vmseq (ux, RVV_ATANPIF_FIX_POS_4_U, vlen), vlen);
  fix_pos = __riscv_vmor
    (fix_pos, __riscv_vmseq (ux, RVV_ATANPIF_FIX_POS_5_U, vlen), vlen);
  fix_pos = __riscv_vmor
    (fix_pos, __riscv_vmseq (ux, RVV_ATANPIF_FIX_POS_6_U, vlen), vlen);

  VFLOAT bias = VFMV_VF (RVV_ATANPIF_BIAS_F, vlen);
  VFLOAT r_pos = __riscv_vfadd (r, bias, vlen);
  r = __riscv_vmerge (r, r_pos, fix_pos, vlen);
  VFLOAT r_neg = __riscv_vfsub (r, bias, vlen);
  r = __riscv_vmerge (r, r_neg, fix_neg, vlen);
  r = __riscv_vmerge (r, vy_special, special_args, vlen);

  RESTORE_FRM;
  return r;
}

#undef LMUL
#undef SIMDLEN
