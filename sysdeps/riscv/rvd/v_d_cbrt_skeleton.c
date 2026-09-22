/* Double-precision vector cbrt function.
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
  VFLOAT vx_orig, vx, vy, vy_special;
  VBOOL special_args;
  VINT n_adjust;

  SET_ROUNDTONEAREST;
  vx_orig = x;
  vy_special = vx_orig;
  n_adjust = VMVI_VX (0, vlen);

  RVV_LIBMVEC_IMPL_NAME (rvv_cbrt_exception_handling, LMUL)
    (&vx_orig, &special_args, &vy_special, &n_adjust, vlen);
  VBOOL active_args = __riscv_vmnot (special_args, vlen);
  vx = __riscv_vfsgnj (vx_orig, fp_posOne, vlen);

  VINT n = __riscv_vsra (F_AS_I (vx), MAN_LEN, vlen);
  n = __riscv_vsub (n, EXP_BIAS, vlen);
  vx = I_AS_F (__riscv_vsub (F_AS_I (vx), __riscv_vsll (n, MAN_LEN, vlen),
                             vlen));

  VFLOAT poly_left = RVV_LIBMVEC_IMPL_NAME (rvv_cbrt_pstep_masked, LMUL) (
      active_args, 0x1.c7feaf5d6cc3bp+0, vx,
      RVV_LIBMVEC_IMPL_NAME (rvv_cbrt_pstep_masked, LMUL) (
          active_args, -0x1.910e22c54a1eap+0, vx,
          VFMV_VF (0x1.3e9d3512b6a5ap+0, vlen), vlen),
      vlen);

  VFLOAT xcube = VFMUL_VV_M (active_args, vx, vx, vlen);
  xcube = VFMUL_VV_M (active_args, xcube, vx, vlen);

  VFLOAT poly_right = RVV_LIBMVEC_IMPL_NAME (rvv_cbrt_pstep_masked, LMUL) (
      active_args, -0x1.3261c716ecf2dp-1, vx,
      RVV_LIBMVEC_IMPL_NAME (rvv_cbrt_pstep_masked, LMUL) (
          active_args, 0x1.3ffc61ff0985dp-3, vx,
          VFMV_VF (-0x1.173278cb4b00fp-6, vlen), vlen),
      vlen);

  VFLOAT z = VFMADD_VV_M (active_args, poly_right, xcube, poly_left, vlen);

  VINT l = __riscv_vadd (n, 1023, vlen);
  VINT m = __riscv_vsra (__riscv_vmul (l, 1366, vlen), 12, vlen);
  VINT j = __riscv_vadd (__riscv_vsll (m, 1, vlen), m, vlen);
  j = __riscv_vsub (l, j, vlen);
  m = __riscv_vsub (m, 341, vlen);

  VINT r = VMVI_VX (RVV_ONE_Q62, vlen);
  VBOOL j_is_1 = __riscv_vmseq (j, 1, vlen);
  VBOOL j_is_2 = __riscv_vmseq (j, 2, vlen);
  r = __riscv_vmerge (r, RVV_CBRT_2_Q62, j_is_1, vlen);
  r = __riscv_vmerge (r, RVV_CBRT_4_Q62, j_is_2, vlen);

  VFLOAT a = VFMUL_VV_M (active_args, z, z, vlen);
  VFLOAT b = VFMUL_VV_M (active_args, vx, z, vlen);
  b = VFMUL_VV_M (active_args, b, a, vlen);
  VFLOAT c = VFMV_VF (RVV_ONE_OVER_3_D, vlen);
  VFLOAT delta = VFNMSUB_VV_M (active_args, b,
                               VFMV_VF (RVV_ONE_OVER_3_D, vlen), c, vlen);
  z = VFMACC_VV_M (active_args, z, delta, z, vlen);

  VFLOAT z_scaled = VFMUL_VV_M (active_args, z, VFMV_VF (0x1.0p62, vlen),
                                vlen);
  z_scaled = __riscv_vfmerge (z_scaled, fp_posZero, special_args, vlen);
  VFLOAT x_scaled = VFMUL_VV_M (active_args, vx,
                                VFMV_VF (0x1.0p62, vlen), vlen);
  x_scaled = __riscv_vfmerge (x_scaled, fp_posZero, special_args, vlen);
  VINT z_q62 = __riscv_vfcvt_x (z_scaled, vlen);
  VINT x_q62 = __riscv_vfcvt_x (x_scaled, vlen);
  VINT a_q = __riscv_vsll (__riscv_vsmul (z_q62, z_q62, 1, vlen), 1, vlen);
  VINT b_q = __riscv_vsll (__riscv_vsmul (x_q62, z_q62, 1, vlen), 1, vlen);
  b_q = __riscv_vsll (__riscv_vsmul (a_q, b_q, 1, vlen), 1, vlen);
  b_q = __riscv_vsll (__riscv_vsmul (b_q, RVV_THIRD_Q62, 1, vlen), 1,
                      vlen);
  VINT delta_q = __riscv_vrsub (b_q, RVV_THIRD_Q62, vlen);
  a_q = __riscv_vsll (__riscv_vsmul (delta_q, z_q62, 1, vlen), 1, vlen);
  z_q62 = __riscv_vadd (z_q62, a_q, vlen);

  z_q62 = __riscv_vsll (__riscv_vsmul (z_q62, z_q62, 1, vlen), 1, vlen);
  z_q62 = __riscv_vsll (__riscv_vsmul (z_q62, x_q62, 1, vlen), 1, vlen);
  r = __riscv_vsmul (r, z_q62, 1, vlen);

  m = __riscv_vadd (m, EXP_BIAS - 61, vlen);
  m = __riscv_vadd (m, n_adjust, vlen);
  VFLOAT scale = I_AS_F (__riscv_vsll (m, MAN_LEN, vlen));
  r = __riscv_vmerge (r, 0, special_args, vlen);
  vy = __riscv_vfcvt_f (r, vlen);
  vy = VFMUL_VV_M (active_args, vy, scale, vlen);
  vy = __riscv_vfsgnj (vy, vx_orig, vlen);
  vy = __riscv_vmerge (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
