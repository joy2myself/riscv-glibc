/* Double-precision vector erfinv function.
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
  VFLOAT vx = x;
  VFLOAT vx_orig = vx;
  VFLOAT vy, vy_special = VFMV_VF (fp_posZero, vlen);
  VBOOL special_args;

  SET_ROUNDTONEAREST;

  RVV_LIBMVEC_IMPL_NAME (rvv_erfinv_exception_handling, LMUL)
    (&vx, &special_args, &vy_special, vlen);

  vx = __riscv_vfsgnj (vx, fp_posOne, vlen);
  VBOOL x_in_right = __riscv_vmfge (vx, RVV_ERFINV_RIGHT_THRESHOLD_D, vlen);

  VFLOAT w_hi_right = VFMV_VF (fp_posZero, vlen);
  VFLOAT w_lo_right = VFMV_VF (fp_posZero, vlen);
  VINT t_right = VMVI_VX (0, vlen);

  if (__riscv_vcpop (x_in_right, vlen) > 0)
    {
      VFLOAT vx_right
        = VFMERGE_VFM (vx, RVV_ERFINV_RIGHT_SAFE_X_D,
                       __riscv_vmnot (x_in_right, vlen), vlen);
      VFLOAT one_minus_x = __riscv_vfrsub (vx_right, fp_posOne, vlen);
      VINT n_adjust = VMVI_VX (RVV_ERFINV_N_ADJUST, vlen);
      VFLOAT y_hi, y_lo;
      RVV_LIBMVEC_IMPL_NAME (rvv_erfinv_neg_logx, LMUL)
        (one_minus_x, n_adjust, &y_hi, &y_lo, vlen);
      RVV_LIBMVEC_IMPL_NAME (rvv_erfinv_sqrtx, LMUL)
        (y_hi, y_lo, &w_hi_right, &w_lo_right, &t_right, RVV_Q63_SCALE_D,
         RVV_Q63_INV_SCALE_D, vlen);
    }

  VINT t = __riscv_vfcvt_x (__riscv_vfmul (vx, RVV_Q63_SCALE_D, vlen),
                            vlen);
  t = VMERGE_VVM_I (t, t_right, x_in_right, vlen);

  VFLOAT w_hi = vx;
  VFLOAT w_lo = VFMV_VF (fp_posZero, vlen);
  w_hi = VMERGE_VVM_F (w_hi, w_hi_right, x_in_right, vlen);
  w_lo = VMERGE_VVM_F (w_lo, w_lo_right, x_in_right, vlen);

  VINT p = VMVI_VX (RVV_ERFINV_P_LEFT10_I, vlen);
  p = VMERGE_VXM_I (p, RVV_ERFINV_P_RIGHT10_I, x_in_right, vlen);
  p = PSTEP_I_ab (x_in_right, RVV_ERFINV_P_RIGHT6_I,
                  RVV_ERFINV_P_LEFT6_I, t,
                  PSTEP_I_ab (x_in_right, RVV_ERFINV_P_RIGHT7_I,
                              RVV_ERFINV_P_LEFT7_I,
                              t,
                              PSTEP_I_ab (x_in_right, RVV_ERFINV_P_RIGHT8_I,
                                          RVV_ERFINV_P_LEFT8_I, t,
                                          PSTEP_I_ab (x_in_right,
                                                      RVV_ERFINV_P_RIGHT9_I,
                                                      RVV_ERFINV_P_LEFT9_I,
                                                      t, p,
                                                      vlen),
                                          vlen),
                              vlen),
                  vlen);

  VINT q = VMVI_VX (RVV_ERFINV_Q_LEFT10_I, vlen);
  q = VMERGE_VXM_I (q, RVV_ERFINV_Q_RIGHT10_I, x_in_right, vlen);
  q = PSTEP_I_ab (x_in_right, RVV_ERFINV_Q_RIGHT6_I,
                  RVV_ERFINV_Q_LEFT6_I, t,
                  PSTEP_I_ab (x_in_right, RVV_ERFINV_Q_RIGHT7_I,
                              RVV_ERFINV_Q_LEFT7_I,
                              t,
                              PSTEP_I_ab (x_in_right, RVV_ERFINV_Q_RIGHT8_I,
                                          RVV_ERFINV_Q_LEFT8_I, t,
                                          PSTEP_I_ab (x_in_right,
                                                      RVV_ERFINV_Q_RIGHT9_I,
                                                      RVV_ERFINV_Q_LEFT9_I,
                                                      t, q,
                                                      vlen),
                                          vlen),
                              vlen),
                  vlen);

  p = PSTEP_I_ab (x_in_right, RVV_ERFINV_P_RIGHT0_I,
                  RVV_ERFINV_P_LEFT0_I, t,
                  PSTEP_I_ab (x_in_right, RVV_ERFINV_P_RIGHT1_I,
                              RVV_ERFINV_P_LEFT1_I,
                              t,
                              PSTEP_I_ab (x_in_right, RVV_ERFINV_P_RIGHT2_I,
                                          RVV_ERFINV_P_LEFT2_I, t,
                                          PSTEP_I_ab (x_in_right,
                                                      RVV_ERFINV_P_RIGHT3_I,
                                                      RVV_ERFINV_P_LEFT3_I,
                                                      t,
                                                      PSTEP_I_ab (
                                                          x_in_right,
                                                          RVV_ERFINV_P_RIGHT4_I,
                                                          RVV_ERFINV_P_LEFT4_I,
                                                          t,
                                                          PSTEP_I_ab (
                                                              x_in_right,
                                                              RVV_ERFINV_P_RIGHT5_I,
                                                              RVV_ERFINV_P_LEFT5_I,
                                                              t, p, vlen),
                                                          vlen),
                                                      vlen),
                                          vlen),
                              vlen),
                  vlen);

  q = PSTEP_I_ab (x_in_right, RVV_ERFINV_Q_RIGHT0_I,
                  RVV_ERFINV_Q_LEFT0_I, t,
                  PSTEP_I_ab (x_in_right, RVV_ERFINV_Q_RIGHT1_I,
                              RVV_ERFINV_Q_LEFT1_I,
                              t,
                              PSTEP_I_ab (x_in_right, RVV_ERFINV_Q_RIGHT2_I,
                                          RVV_ERFINV_Q_LEFT2_I, t,
                                          PSTEP_I_ab (x_in_right,
                                                      RVV_ERFINV_Q_RIGHT3_I,
                                                      RVV_ERFINV_Q_LEFT3_I,
                                                      t,
                                                      PSTEP_I_ab (
                                                          x_in_right,
                                                          RVV_ERFINV_Q_RIGHT4_I,
                                                          RVV_ERFINV_Q_LEFT4_I,
                                                          t,
                                                          PSTEP_I_ab (
                                                              x_in_right,
                                                              RVV_ERFINV_Q_RIGHT5_I,
                                                              RVV_ERFINV_Q_LEFT5_I,
                                                              t, q, vlen),
                                                          vlen),
                                                      vlen),
                                          vlen),
                              vlen),
                  vlen);

  VFLOAT p_hi = __riscv_vfcvt_f (p, vlen);
  VFLOAT p_lo
    = __riscv_vfcvt_f (__riscv_vsub (p, __riscv_vfcvt_x (p_hi, vlen), vlen),
                       vlen);
  VFLOAT delta_p0 = VFMV_VF (RVV_ERFINV_DELTA_P0_LEFT_D, vlen);
  delta_p0 = VFMERGE_VFM (delta_p0, RVV_ERFINV_DELTA_P0_RIGHT_D,
                          x_in_right, vlen);
  p_lo = __riscv_vfadd (p_lo, delta_p0, vlen);

  VFLOAT q_hi = __riscv_vfcvt_f (q, vlen);
  VFLOAT q_lo
    = __riscv_vfcvt_f (__riscv_vsub (q, __riscv_vfcvt_x (q_hi, vlen), vlen),
                       vlen);
  VFLOAT delta_q0 = VFMV_VF (RVV_ERFINV_DELTA_Q0_LEFT_D, vlen);
  delta_q0 = VFMERGE_VFM (delta_q0, RVV_ERFINV_DELTA_Q0_RIGHT_D,
                          x_in_right, vlen);
  q_lo = __riscv_vfadd (q_lo, delta_q0, vlen);

  VFLOAT y_hi = __riscv_vfmul (w_hi, p_hi, vlen);
  VFLOAT y_lo = __riscv_vfmsub (w_hi, p_hi, y_hi, vlen);
  y_lo = __riscv_vfmacc (y_lo, w_hi, p_lo, vlen);
  y_lo = __riscv_vfmacc (y_lo, w_lo, p_hi, vlen);

  DIV_N2D2 (y_hi, y_lo, q_hi, q_lo, w_hi, vlen);
  vy = w_hi;

  vy = __riscv_vfsgnj (vy, vx_orig, vlen);
  vy = VMERGE_VVM_F (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
