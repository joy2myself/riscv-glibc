/* Double-precision vector erfc function.
   Copyright (C) 2026 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public License as
   published by the Free Software Foundation; either version 2.1 of the
   License, or (at your option) any later version.

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
  VFLOAT vx = x;
  VFLOAT vx_orig = vx;
  VFLOAT vy, vy_special = VFMV_VF (fp_posZero, vlen);
  VBOOL special_args;

  SET_ROUNDTONEAREST;

  RVV_ERFC_EXCEPTION_HANDLING (&vx, &special_args, &vy_special, vlen);

  vx = __riscv_vfsgnj (vx, fp_posOne, vlen);
  vx = __riscv_vfmin (vx, RVV_ERFC_MAX_CLIP_D, vlen);
  vx = __riscv_vfmax (vx, RVV_ERFC_MIN_CLIP_D, vlen);

  VINT r = RVV_ERFC_X_TRANSFORM (vx, vlen);

  VINT n, a;
  vy = vx;
  RVV_ERFC_EXP_NEG_AB (vx, vy, &n, &a, vlen);

  VINT p_right = PSTEP_I (
      RVV_ERFC_P16_I, r,
      PSTEP_I (
          RVV_ERFC_P17_I, r,
          PSTEP_I (RVV_ERFC_P18_I, r,
                   PSTEP_I (RVV_ERFC_P19_I, r,
                            PSTEP_I (RVV_ERFC_P20_I, r,
                                     PSTEP_I (RVV_ERFC_P21_I, r,
                                              PSTEP_I (RVV_ERFC_P22_I,
                                                       RVV_ERFC_P23_I,
                                                       r, vlen),
                                              vlen),
                                     vlen),
                            vlen),
                   vlen),
          vlen),
      vlen);

  VINT r_8 = __riscv_vsmul (r, r, 1, vlen);
  r_8 = __riscv_vsmul (r_8, r_8, 1, vlen);
  r_8 = __riscv_vsmul (r_8, r_8, 1, vlen);

  VINT p_mid = PSTEP_I (
      RVV_ERFC_P8_I, r,
      PSTEP_I (
          RVV_ERFC_P9_I, r,
          PSTEP_I (RVV_ERFC_P10_I, r,
                   PSTEP_I (RVV_ERFC_P11_I, r,
                            PSTEP_I (RVV_ERFC_P12_I, r,
                                     PSTEP_I (RVV_ERFC_P13_I, r,
                                              PSTEP_I (RVV_ERFC_P14_I,
                                                       RVV_ERFC_P15_I,
                                                       r, vlen),
                                              vlen),
                                     vlen),
                            vlen),
                   vlen),
          vlen),
      vlen);

  p_right = __riscv_vsmul (r_8, p_right, 1, vlen);
  p_right = __riscv_vadd (p_right, p_mid, vlen);
  p_right = __riscv_vsmul (r_8, p_right, 1, vlen);

  VINT p_left = PSTEP_I (
      RVV_ERFC_P0_I, r,
      PSTEP_I (
          RVV_ERFC_P1_I, r,
          PSTEP_I (RVV_ERFC_P2_I, r,
                   PSTEP_I (RVV_ERFC_P3_I, r,
                            PSTEP_I (RVV_ERFC_P4_I, r,
                                     PSTEP_I (RVV_ERFC_P5_I, r,
                                              PSTEP_I (RVV_ERFC_P6_I,
                                                       RVV_ERFC_P7_I,
                                                       r, vlen),
                                              vlen),
                                     vlen),
                            vlen),
                   vlen),
          vlen),
      vlen);

  VINT p = __riscv_vadd (p_left, p_right, vlen);

  VINT m, b;
  RVV_ERFC_RECIP_SCALE (vx, &b, &m, vlen);

  p = __riscv_vsmul (p, a, 1, vlen);
  p = __riscv_vsmul (p, b, 1, vlen);
  n = __riscv_vsub (n, m, vlen);
  n = __riscv_vadd (n, 2, vlen);

  VUINT ell = I_AS_U (__riscv_vrsub (n, -60, vlen));
  ell = __riscv_vminu (ell, 63, vlen);
  VINT pp = __riscv_vsra (p, ell, vlen);
  VINT q = VMVI_VX (1, vlen);
  q = __riscv_vsll (q, 61, vlen);
  q = __riscv_vsub (q, pp, vlen);
  VFLOAT vz = __riscv_vfcvt_f (q, vlen);
  vz = __riscv_vfmul (vz, RVV_Q60_SCALE_D, vlen);

  vy = __riscv_vfcvt_f (p, vlen);
  FAST_LDEXP (vy, n, vlen);

  VBOOL x_is_neg = __riscv_vmflt (vx_orig, fp_posZero, vlen);
  vy = __riscv_vmerge (vy, vz, x_is_neg, vlen);
  vy = __riscv_vmerge (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
