/* Double-precision vector erfcinv function.
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
  VFLOAT vx_sign, vy, vy_special = VFMV_VF (fp_posZero, vlen);
  VBOOL special_args;

  SET_ROUNDTONEAREST;

  RVV_LIBMVEC_IMPL_NAME (rvv_erfcinv_exception_handling, LMUL)
    (&vx, &special_args, &vy_special, vlen);

  vx_sign = __riscv_vfrsub (vx, fp_posOne, vlen);
  VFLOAT two_minus_x = __riscv_vfadd (vx_sign, fp_posOne, vlen);
  VBOOL x_gt_1 = __riscv_vmflt (vx_sign, fp_posZero, vlen);
  vx = VMERGE_VVM_F (vx, two_minus_x, x_gt_1, vlen);

  VBOOL x_in_left = __riscv_vmfle (vx, RVV_ERFCINV_LEFT_THRESHOLD_D, vlen);

  VFLOAT w_hi_left = VFMV_VF (fp_posZero, vlen);
  VFLOAT w_lo_left = VFMV_VF (fp_posZero, vlen);
  VINT t_left = VMVI_VX (0, vlen);
  VINT t_tiny = VMVI_VX (0, vlen);
  VBOOL x_is_tiny = __riscv_vmsne (VMVU_VX (0, vlen), 0, vlen);

  if (__riscv_vcpop (x_in_left, vlen) > 0)
    {
      VFLOAT x_left = VFMV_VF (RVV_ERFCINV_LEFT_SEED_D, vlen);
      x_left = VMERGE_VVM_F (x_left, vx, x_in_left, vlen);
      x_is_tiny = __riscv_vmflt (x_left, RVV_ERFCINV_TINY_THRESHOLD_D, vlen);
      x_left = __riscv_vfmul (x_left, RVV_ERFCINV_SCALE_X_D, vlen);
      VFLOAT y_hi, y_lo;
      RVV_LIBMVEC_IMPL_NAME (rvv_erfcinv_neg_logx, LMUL)
        (x_left, RVV_ERFCINV_N_ADJUST, &y_hi, &y_lo, vlen);
      RVV_LIBMVEC_IMPL_NAME (rvv_erfcinv_sqrtx, LMUL)
        (y_hi, y_lo, &w_hi_left, &w_lo_left, &t_left, RVV_Q63_SCALE_D,
         RVV_Q63_INV_SCALE_D, vlen);
      if (__riscv_vcpop (x_is_tiny, vlen) > 0)
        {
          /* The 2^64 scale can overflow signed integer conversion for
             discarded non-tiny lanes.  Keep the source graph unchanged on
             tiny lanes while giving those lanes a defined, in-range input.  */
          VFLOAT y_hi_tiny
            = VFMV_VF (RVV_ERFCINV_TINY_QUARANTINE_Y_D, vlen);
          VFLOAT y_lo_tiny = VFMV_VF (fp_posZero, vlen);
          y_hi_tiny = VMERGE_VVM_F (y_hi_tiny, y_hi, x_is_tiny, vlen);
          y_lo_tiny = VMERGE_VVM_F (y_lo_tiny, y_lo, x_is_tiny, vlen);
          VFLOAT w_hi_dummy = VFMV_VF (fp_posZero, vlen);
          VFLOAT w_lo_dummy = VFMV_VF (fp_posZero, vlen);
          RVV_LIBMVEC_IMPL_NAME (rvv_erfcinv_sqrtx, LMUL)
            (y_hi_tiny, y_lo_tiny, &w_hi_dummy, &w_lo_dummy, &t_tiny,
             RVV_ERFCINV_TINY_T_SCALE_D,
             RVV_ERFCINV_TINY_T_SCALE_INV_D, vlen);
        }
    }

  VFLOAT w_hi = __riscv_vfrsub (vx, fp_posOne, vlen);
  VFLOAT w_lo = __riscv_vfrsub (w_hi, fp_posOne, vlen);
  w_lo = __riscv_vfsub (w_lo, vx, vlen);

  VINT t = __riscv_vfcvt_x (__riscv_vfmul (w_hi, RVV_Q63_SCALE_D, vlen),
                            vlen);
  VFLOAT delta_t = __riscv_vfmul (w_lo, RVV_Q63_SCALE_D, vlen);
  t = __riscv_vadd (t, __riscv_vfcvt_x (delta_t, vlen), vlen);

  t = VMERGE_VVM_I (t, t_left, x_in_left, vlen);
  w_hi = VMERGE_VVM_F (w_hi, w_hi_left, x_in_left, vlen);
  w_lo = VMERGE_VVM_F (w_lo, w_lo_left, x_in_left, vlen);

  VINT p = VMVI_VX (RVV_ERFCINV_P_RIGHT10_I, vlen);
  p = VMERGE_VXM_I (p, RVV_ERFCINV_P_LEFT10_I, x_in_left, vlen);
  p = PSTEP_I_ab (
    x_in_left, RVV_ERFCINV_P_LEFT6_I, RVV_ERFCINV_P_RIGHT6_I, t,
    PSTEP_I_ab (
      x_in_left, RVV_ERFCINV_P_LEFT7_I, RVV_ERFCINV_P_RIGHT7_I, t,
      PSTEP_I_ab (
        x_in_left, RVV_ERFCINV_P_LEFT8_I, RVV_ERFCINV_P_RIGHT8_I, t,
        PSTEP_I_ab (x_in_left, RVV_ERFCINV_P_LEFT9_I,
                    RVV_ERFCINV_P_RIGHT9_I, t, p,
                    vlen),
        vlen),
      vlen),
    vlen);

  VINT q = VMVI_VX (RVV_ERFCINV_Q_RIGHT10_I, vlen);
  q = VMERGE_VXM_I (q, RVV_ERFCINV_Q_LEFT10_I, x_in_left, vlen);
  q = PSTEP_I_ab (
    x_in_left, RVV_ERFCINV_Q_LEFT6_I, RVV_ERFCINV_Q_RIGHT6_I, t,
    PSTEP_I_ab (
      x_in_left, RVV_ERFCINV_Q_LEFT7_I, RVV_ERFCINV_Q_RIGHT7_I, t,
      PSTEP_I_ab (
        x_in_left, RVV_ERFCINV_Q_LEFT8_I, RVV_ERFCINV_Q_RIGHT8_I, t,
        PSTEP_I_ab (x_in_left, RVV_ERFCINV_Q_LEFT9_I,
                    RVV_ERFCINV_Q_RIGHT9_I, t, q,
                    vlen),
        vlen),
      vlen),
    vlen);

  p = PSTEP_I_ab (
    x_in_left, RVV_ERFCINV_P_LEFT0_I, RVV_ERFCINV_P_RIGHT0_I, t,
    PSTEP_I_ab (
      x_in_left, RVV_ERFCINV_P_LEFT1_I, RVV_ERFCINV_P_RIGHT1_I, t,
      PSTEP_I_ab (
        x_in_left, RVV_ERFCINV_P_LEFT2_I, RVV_ERFCINV_P_RIGHT2_I, t,
        PSTEP_I_ab (
          x_in_left, RVV_ERFCINV_P_LEFT3_I, RVV_ERFCINV_P_RIGHT3_I, t,
          PSTEP_I_ab (
            x_in_left, RVV_ERFCINV_P_LEFT4_I, RVV_ERFCINV_P_RIGHT4_I, t,
            PSTEP_I_ab (x_in_left, RVV_ERFCINV_P_LEFT5_I,
                        RVV_ERFCINV_P_RIGHT5_I, t, p,
                        vlen),
            vlen),
          vlen),
        vlen),
      vlen),
    vlen);

  q = PSTEP_I_ab (
    x_in_left, RVV_ERFCINV_Q_LEFT0_I, RVV_ERFCINV_Q_RIGHT0_I, t,
    PSTEP_I_ab (
      x_in_left, RVV_ERFCINV_Q_LEFT1_I, RVV_ERFCINV_Q_RIGHT1_I, t,
      PSTEP_I_ab (
        x_in_left, RVV_ERFCINV_Q_LEFT2_I, RVV_ERFCINV_Q_RIGHT2_I, t,
        PSTEP_I_ab (
          x_in_left, RVV_ERFCINV_Q_LEFT3_I, RVV_ERFCINV_Q_RIGHT3_I, t,
          PSTEP_I_ab (
            x_in_left, RVV_ERFCINV_Q_LEFT4_I, RVV_ERFCINV_Q_RIGHT4_I, t,
            PSTEP_I_ab (x_in_left, RVV_ERFCINV_Q_LEFT5_I,
                        RVV_ERFCINV_Q_RIGHT5_I, t, q,
                        vlen),
            vlen),
          vlen),
        vlen),
      vlen),
    vlen);

  VFLOAT p_hi = __riscv_vfcvt_f (p, vlen);
  VFLOAT p_lo
    = __riscv_vfcvt_f (__riscv_vsub (p, __riscv_vfcvt_x (p_hi, vlen), vlen),
                       vlen);
  VFLOAT delta_p0 = VFMV_VF (RVV_ERFCINV_DELTA_P0_RIGHT_D, vlen);
  delta_p0 = VFMERGE_VFM (delta_p0, RVV_ERFCINV_DELTA_P0_LEFT_D,
                          x_in_left, vlen);
  p_lo = __riscv_vfadd (p_lo, delta_p0, vlen);

  VFLOAT q_hi = __riscv_vfcvt_f (q, vlen);
  VFLOAT q_lo
    = __riscv_vfcvt_f (__riscv_vsub (q, __riscv_vfcvt_x (q_hi, vlen), vlen),
                       vlen);
  VFLOAT delta_q0 = VFMV_VF (RVV_ERFCINV_DELTA_Q0_RIGHT_D, vlen);
  delta_q0 = VFMERGE_VFM (delta_q0, RVV_ERFCINV_DELTA_Q0_LEFT_D,
                          x_in_left, vlen);
  q_lo = __riscv_vfadd (q_lo, delta_q0, vlen);

  if (__riscv_vcpop (x_is_tiny, vlen) > 0)
    {
      VFLOAT p_hi_tiny, p_lo_tiny, q_hi_tiny, q_lo_tiny;
      RVV_LIBMVEC_IMPL_NAME (rvv_erfcinv_pq_hilo_tiny, LMUL)
        (t_tiny, &p_hi_tiny, &p_lo_tiny, &q_hi_tiny, &q_lo_tiny,
         vlen);
      p_hi = VMERGE_VVM_F (p_hi, p_hi_tiny, x_is_tiny, vlen);
      p_lo = VMERGE_VVM_F (p_lo, p_lo_tiny, x_is_tiny, vlen);
      q_hi = VMERGE_VVM_F (q_hi, q_hi_tiny, x_is_tiny, vlen);
      q_lo = VMERGE_VVM_F (q_lo, q_lo_tiny, x_is_tiny, vlen);
    }

  VFLOAT y_hi = __riscv_vfmul (w_hi, p_hi, vlen);
  VFLOAT y_lo = __riscv_vfmsub (w_hi, p_hi, y_hi, vlen);
  y_lo = __riscv_vfmacc (y_lo, w_hi, p_lo, vlen);
  y_lo = __riscv_vfmacc (y_lo, w_lo, p_hi, vlen);

  DIV_N2D2 (y_hi, y_lo, q_hi, q_lo, w_hi, vlen);

  vy = w_hi;
  vy = __riscv_vfsgnj (vy, vx_sign, vlen);
  vy = VMERGE_VVM_F (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
