/* Single-precision vector atanh function.
   atanh(x) = (1/2) * log((1+x)/(1-x)) for |x| < 1.

   Based on SLEEF xatanhf implementation.

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

/* Shared constants and importer-qualified helper names live in v_math.h.  */

static inline VINT
RVV_LIBMVEC_ATANHF_VILOGBK (VFLOAT d, size_t vlen)
{
  VBOOL small = __riscv_vmflt (d, 5.421010862427522e-20f, vlen);
  VFLOAT scaled = VFMUL_VF (d, 1.8446744073709552e19f, vlen);
  d = __riscv_vmerge (d, scaled, small, vlen);
  VUINT q = __riscv_vsrl (F_AS_U (d), MAN_LEN_F, vlen);
  q = __riscv_vand (q, 0xff, vlen);
  VINT bias = __riscv_vmerge (VMVI_VX (0x7f, vlen),
                              VMVI_VX (64 + 0x7f, vlen), small, vlen);
  return __riscv_vsub (U_AS_I (q), bias, vlen);
}

static inline VFLOAT
RVV_LIBMVEC_ATANHF_VPOW2I (VINT q, size_t vlen)
{
  VINT biased = __riscv_vadd (q, EXP_BIAS_F, vlen);
  return I_AS_F (__riscv_vsll (biased, MAN_LEN_F, vlen));
}

static inline void
RVV_LIBMVEC_ATANHF_LOGK2F
  (VFLOAT d_hi, VFLOAT d_lo, VFLOAT *r_hi, VFLOAT *r_lo,
   VBOOL active, size_t vlen)
{
  /* This is the pinned SLEEF logk2f graph.  Only the first normalization
     multiply is masked: clean has already sanitized discarded lanes to d=1,
     and that multiply is the only neutral-path operation that is inexact for
     such a lane.  The remaining d=1 path is exact and supplies harmless
     values to the common vector operations.  */
  VFLOAT scaled_hi = VFMUL_VF_M (active, d_hi, 1.0f / 0.75f, vlen);
  VINT e = RVV_LIBMVEC_ATANHF_VILOGBK (scaled_hi, vlen);
  VFLOAT scale = RVV_LIBMVEC_ATANHF_VPOW2I
    (__riscv_vneg (e, vlen), vlen);

  VFLOAT m_hi, m_lo;
  RVV_SLEEF_ORIG_DFSCALE_VF2_VF
    (d_hi, d_lo, scale, m_hi, m_lo, vlen);

  VFLOAT n_hi, n_lo, q_hi, q_lo;
  RVV_SLEEF_ORIG_DFADD2_VF2_VF
    (m_hi, m_lo, VFMV_VF (-1.0f, vlen), n_hi, n_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF
    (m_hi, m_lo, VFMV_VF (1.0f, vlen), q_hi, q_lo, vlen);

  VFLOAT x_hi, x_lo;
  RVV_SLEEF_ORIG_DFDIV_VF2_VF2
    (n_hi, n_lo, q_hi, q_lo, x_hi, x_lo, vlen);
  VFLOAT x2_hi, x2_lo;
  RVV_SLEEF_ORIG_DFSQU_VF2 (x_hi, x_lo, x2_hi, x2_lo, vlen);

  VFLOAT t = VFMV_VF (RVV_ATANHF_LOG_POLY0_F, vlen);
  t = VFMADD_VV (t, x2_hi, VFMV_VF (RVV_ATANHF_LOG_POLY1_F, vlen), vlen);
  t = VFMADD_VV (t, x2_hi, VFMV_VF (RVV_ATANHF_LOG_POLY2_F, vlen), vlen);
  t = VFMADD_VV (t, x2_hi, VFMV_VF (RVV_ATANHF_LOG_POLY3_F, vlen), vlen);

  VFLOAT e_f = __riscv_vfcvt_f (e, vlen);
  VFLOAT s_hi, s_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (VFMV_VF (RVV_LOG1PF_LN2_HI_F, vlen),
     VFMV_VF (RVV_LOG1PF_LN2_LO_F, vlen), e_f,
     s_hi, s_lo, vlen);

  VFLOAT two_x_hi, two_x_lo;
  RVV_SLEEF_ORIG_DFSCALE_VF2_VF
    (x_hi, x_lo, VFMV_VF (2.0f, vlen), two_x_hi, two_x_lo, vlen);
  RVV_SLEEF_ORIG_DFADD_VF2_VF2
    (s_hi, s_lo, two_x_hi, two_x_lo, s_hi, s_lo, vlen);

  VFLOAT x2x_hi, x2x_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (x2_hi, x2_lo, x_hi, x_lo, x2x_hi, x2x_lo, vlen);
  VFLOAT correction_hi, correction_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (x2x_hi, x2x_lo, t, correction_hi, correction_lo, vlen);
  RVV_SLEEF_ORIG_DFADD_VF2_VF2
    (s_hi, s_lo, correction_hi, correction_lo, *r_hi, *r_lo, vlen);
}

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT vx, vy, vy_special;
  VBOOL special_args;

  SET_ROUNDTONEAREST;
  vx = x;

  VFLOAT y = __riscv_vfabs (vx, vlen);
  VUINT vclass = __riscv_vfclass (vx, vlen);
  VBOOL x_nan;
  IDENTIFY_F (vclass, class_NaN, x_nan, vlen);

  VFLOAT y_cmp = __riscv_vfmerge (y, fp_posZero, x_nan, vlen);
  VBOOL y_gt_1 = __riscv_vmfgt (y_cmp, 1.0f, vlen);
  VBOOL y_eq_1 = __riscv_vmfeq (y_cmp, 1.0f, vlen);

  special_args = __riscv_vmor (y_gt_1, y_eq_1, vlen);
  special_args = __riscv_vmor (special_args, x_nan, vlen);

  vy_special = vx;
  if (__riscv_vcpop (special_args, vlen) > 0)
    {
      /* Quiet invalid-domain lanes explicitly instead of relying on core
         divide/log accidents.  */
      VFLOAT invalid_operand = VFMV_VF (fp_sNaN, vlen);
      VFLOAT invalid_result = __riscv_vfadd (y_gt_1, invalid_operand,
                                             invalid_operand, vlen);
      vy_special = __riscv_vmerge (vy_special, invalid_result, y_gt_1, vlen);

      /* Preserve the pole at |x| == 1 and raise FE_DIVBYZERO only in those
         lanes.  The signed zero carries the sign of the original input, so
         the provider also constructs the correct signed infinity.  */
      VFLOAT pole_zero
        = __riscv_vfsgnj (VFMV_VF (fp_posZero, vlen), vx, vlen);
      VFLOAT pole_den
        = VMERGE_VVM_F (VFMV_VF (fp_posOne, vlen), pole_zero,
                         y_eq_1, vlen);
      VFLOAT pole_result
        = __riscv_vfdiv (VFMV_VF (fp_posOne, vlen), pole_den, vlen);
      vy_special = __riscv_vmerge (vy_special, pole_result, y_eq_1, vlen);

      /* Propagate NaNs via x + x.  */
      VFLOAT vy_nan = __riscv_vfadd (x_nan, vx, vx, vlen);
      vy_special = __riscv_vmerge (vy_special, vy_nan, x_nan, vlen);

      /* Feed benign zeros to the finite core for repaired lanes.  */
      vx = __riscv_vfmerge (vx, fp_posZero, special_args, vlen);
      y = __riscv_vfabs (vx, vlen);
    }

  /* Do not run the sanitized finite core when every lane is already
     resolved.  This is important for mixed-lane fenv isolation and keeps the
     pole/NaN provider operations above as the only exceptional operations.  */
  if (__riscv_vcpop (special_args, vlen) == vlen)
    {
      RESTORE_FRM;
      return vy_special;
    }

  /* Compute (1+y)/(1-y) using double-float arithmetic */
  VFLOAT one_plus_y_hi, one_plus_y_lo;
  RVV_SLEEF_ORIG_DFADD2_VF_VF
    (VFMV_VF (1.0f, vlen), y, one_plus_y_hi, one_plus_y_lo, vlen);

  VFLOAT one_minus_y_hi, one_minus_y_lo;
  RVV_SLEEF_ORIG_DFADD2_VF_VF
    (VFMV_VF (1.0f, vlen), __riscv_vfneg (y, vlen),
     one_minus_y_hi, one_minus_y_lo, vlen);

  /* d = (1+y) / (1-y) as double-float */
  VFLOAT d_hi, d_lo;
  RVV_SLEEF_ORIG_DFDIV_VF2_VF2
    (one_plus_y_hi, one_plus_y_lo, one_minus_y_hi, one_minus_y_lo,
     d_hi, d_lo, vlen);

  VFLOAT final_hi, final_lo;
  RVV_LIBMVEC_ATANHF_LOGK2F
    (d_hi, d_lo, &final_hi, &final_lo,
     __riscv_vmnot (special_args, vlen), vlen);

  /* atanh(x) = 0.5 * log((1+x)/(1-x)) */
  vy = __riscv_vfmul (__riscv_vfadd (final_hi, final_lo, vlen), 0.5f, vlen);

  /* Restore sign */
  vy = __riscv_vfsgnj (vy, x, vlen);
  vy = __riscv_vmerge (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
