/* Single-precision vector atan function.
   atan(x) is defined for all finite x.

   Based on SLEEF xatanf_u1 implementation.

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

/* atan(x) for all finite x:
   For |x| <= 1: atan(x) ~= x + x^3 * poly(x^2)
   For |x| > 1: atan(x) = sign(x) * pi/2 - atan(1/|x|)  */

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT vx, vy_special;
  VFLOAT r_hi, r_lo;
  VBOOL special_args;

  SET_ROUNDTONEAREST;
  vx = x;

  /* Classify special lanes up front, following the log/logf skeleton style:
     NaN propagates via x+x and +/-Inf rounds directly to +/-pi/2, while the
     main path stays unchanged for finite arguments.  */
  VUINT vclass = __riscv_vfclass (vx, vlen);
  VBOOL x_nan, x_inf;
  IDENTIFY_F (vclass, class_NaN, x_nan, vlen);
  IDENTIFY_F (vclass, class_Inf, x_inf, vlen);
  special_args = __riscv_vmor (x_nan, x_inf, vlen);
  UINT nb_special_args = __riscv_vcpop (special_args, vlen);

  vy_special = vx;
  if (nb_special_args > 0)
    {
      VFLOAT vy_inf = __riscv_vfsgnj (VFMV_VF (RVV_ATANF_PI2_F, vlen), vx,
                                      vlen);
      vy_special = __riscv_vmerge (vy_special, vy_inf, x_inf, vlen);

      VFLOAT vy_nan = __riscv_vfadd (x_nan, vx, vx, vlen);
      vy_special = __riscv_vmerge (vy_special, vy_nan, x_nan, vlen);

      vx = __riscv_vfmerge (vx, fp_posZero, special_args, vlen);
    }

  VFLOAT ax = __riscv_vfabs (vx, vlen);
  RVV_LIBMVEC_IMPL_NAME (atan2kf_u1, LMUL)
    (ax, VFMV_VF (0.0f, vlen), VFMV_VF (1.0f, vlen), VFMV_VF (0.0f, vlen),
     &r_hi, &r_lo, vlen);

  VFLOAT r = __riscv_vfadd (r_hi, r_lo, vlen);
  r = __riscv_vfsgnj (r, vx, vlen);
  r = __riscv_vmerge (r, vy_special, special_args, vlen);

  RESTORE_FRM;
  return r;
}

#undef LMUL
#undef SIMDLEN
