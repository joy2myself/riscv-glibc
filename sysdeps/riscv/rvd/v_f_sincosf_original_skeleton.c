/* Single-precision source-original SLEEF sincosf skeleton.

   This is a target-owned translation of SLEEF XSINCOSF_U1, which is exported
   as Sleef_sincosf_u10 by the pinned source rename table.  The direct source
   reduction (2/pi), rempif large-argument path, double-float reconstruction,
   sine/cosine polynomials, quadrant/sign reconstruction, and late
   special-value behavior are kept in source order.  Only the RVV
   representation, LMUL/VL plumbing, target-owned double-float helpers, and
   glibc importer ABI are adapted.  No SLEEF source file is included,
   compiled, or linked.

   SLEEF is compiled with -ffp-contract=off.  The wrapper's build rule keeps
   that policy for this comparison implementation.

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

   You should have received a copy of the GNU Lesser General Public License
   along with the GNU C Library; if not, see
   <https://www.gnu.org/licenses/>.  */

#include <vecmath_config.h>
#include <v_math.h>

#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

static void
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT d,
                                           VFLOAT *out_sin,
                                           VFLOAT *out_cos, size_t vlen)
{
  VINT q;
  VFLOAT u, v;
  VFLOAT s_hi, s_lo, t_hi, t_lo;
  VFLOAT rx, ry;

  /* Source XSINCOSF_U1 starts with vrint(d * 2/pi), then converts that
     already-rounded value to the integer quadrant.  */
  u = __riscv_vfcvt_f
    (RVV_VFCVT_X_F_RNE
       (VFMUL_VF (d, RVV_SLEEF_TRIG_TWO_OVER_PI, vlen), vlen), vlen);
  q = RVV_VFCVT_X_F_RNE (u, vlen);

  v = VFMADD_VV
    (u, VFMV_VF (RVV_SLEEF_TRIG_NEG_PIBY2_A2_F, vlen), d, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF_VF
    (v, VFMUL_VF (u, RVV_SLEEF_TRIG_NEG_PIBY2_B2_F, vlen),
     s_hi, s_lo, vlen);
  RVV_SLEEF_ORIG_DFADD_VF2_VF2_VF
    (s_hi, s_lo,
     VFMUL_VF (u, RVV_SLEEF_TRIG_NEG_PIBY2_C2_F, vlen),
     s_hi, s_lo, vlen);

  VBOOL g = __riscv_vmflt
    (__riscv_vfabs (d, vlen),
     RVV_SLEEF_TRIG_RANGE_MAX2, vlen);

  if (__riscv_vcpop (g, vlen) != vlen)
    {
      VFLOAT rem_hi, rem_lo;
      VINT rem_q;
      RVV_LIBMVEC_IMPL_NAME (sinf_original_reduce_large, LMUL)
        (d, &rem_hi, &rem_lo, &rem_q, vlen);

      /* XSINCOSF_U1 poisons only the high component for Inf/NaN lanes
         after rempif, preserving its source-native late handling.  */
      VBOOL input_inf = __riscv_vmfeq
        (__riscv_vfabs (d, vlen), fp_posInf_f, vlen);
      VBOOL input_nan = __riscv_vmfne (d, d, vlen);
      VBOOL input_nonfinite = __riscv_vmor (input_inf, input_nan, vlen);
      VUINT poison = __riscv_vmerge
        (VMVU_VX (0, vlen), VMVU_VX (0xffffffffu, vlen),
         input_nonfinite, vlen);
      rem_hi = U_AS_F
        (__riscv_vor (F_AS_U (rem_hi), poison, vlen));

      VBOOL large = __riscv_vmnot (g, vlen);
      q = __riscv_vmerge (q, rem_q, large, vlen);
      s_hi = __riscv_vmerge (s_hi, rem_hi, large, vlen);
      s_lo = __riscv_vmerge (s_lo, rem_lo, large, vlen);
    }

  /* Source keeps the unreduced pair in t and replaces only s.x with
     dfsqu_vf_vf2(s); s.y is intentionally left as the source pair low part.  */
  t_hi = s_hi;
  t_lo = s_lo;
  /* Replace only the high component with the source pair square.  */
  {
    VFLOAT xhxl = VFMUL_VV (s_hi, s_lo, vlen);
    VFLOAT twice_xhxl = VFADD_VV (xhxl, xhxl, vlen);
    s_hi = VFMADD_VV (s_hi, s_hi, twice_xhxl, vlen);
  }

  /* sin(r) = r + r^3 * P(r^2), in source Horner/order.  */
  u = VFMV_VF (RVV_SINCOSF_SIN_POLY0_F, vlen);
  u = VFMADD_VV
    (u, s_hi, VFMV_VF (RVV_SINCOSF_SIN_POLY1_F, vlen), vlen);
  u = VFMADD_VV
    (u, s_hi, VFMV_VF (RVV_SINCOSF_SIN_POLY2_F, vlen), vlen);
  VFLOAT sxt = VFMUL_VV (s_hi, t_hi, vlen);
  sxt = VFMUL_VV (u, sxt, vlen);

  VFLOAT sin_hi, sin_lo;
  RVV_SLEEF_ORIG_DFADD_VF2_VF2_VF
    (t_hi, t_lo, sxt, sin_hi, sin_lo, vlen);
  rx = VFADD_VV (sin_hi, sin_lo, vlen);

  VBOOL negzero = __riscv_vmseq
    (F_AS_U (d), 0x80000000u, vlen);
  rx = VMERGE_VVM_F (rx, VFMV_VF (-0.0f, vlen), negzero, vlen);

  /* cos(r) = 1 + r^2 * Q(r^2), retaining the source pair operations.  */
  u = VFMV_VF (RVV_SINCOSF_COS_POLY0_F, vlen);
  u = VFMADD_VV
    (u, s_hi, VFMV_VF (RVV_SINCOSF_COS_POLY1_F, vlen), vlen);
  u = VFMADD_VV
    (u, s_hi, VFMV_VF (RVV_SINCOSF_COS_POLY2_F, vlen), vlen);
  u = VFMADD_VV
    (u, s_hi, VFMV_VF (RVV_SINCOSF_COS_POLY3_F, vlen), vlen);
  u = VFMADD_VV
    (u, s_hi, VFMV_VF (RVV_SINCOSF_COS_POLY4_F, vlen), vlen);

  VFLOAT cos_term_hi, cos_term_lo;
  RVV_SLEEF_ORIG_DFMUL_VF_VF
    (s_hi, u, cos_term_hi, cos_term_lo, vlen);
  VFLOAT cos_hi, cos_lo;
  RVV_SLEEF_ORIG_DFADD_VF_VF2
    (VFMV_VF (1.0f, vlen), cos_term_hi, cos_term_lo,
     cos_hi, cos_lo, vlen);
  ry = VFADD_VV (cos_hi, cos_lo, vlen);

  /* Quadrant reconstruction and source bitwise sign changes.  */
  VBOOL even = __riscv_vmseq
    (__riscv_vand (q, 1, vlen), 0, vlen);
  VFLOAT sin_result = VMERGE_VVM_F (ry, rx, even, vlen);
  VFLOAT cos_result = VMERGE_VVM_F (rx, ry, even, vlen);

  VBOOL sin_flip = __riscv_vmseq
    (__riscv_vand (q, 2, vlen), 2, vlen);
  VUINT sin_sign = __riscv_vmerge
    (VMVU_VX (0, vlen), VMVU_VX (0x80000000u, vlen), sin_flip, vlen);
  sin_result = U_AS_F
    (__riscv_vxor (F_AS_U (sin_result), sin_sign, vlen));

  VINT q_plus_one = __riscv_vadd (q, 1, vlen);
  VBOOL cos_flip = __riscv_vmseq
    (__riscv_vand (q_plus_one, 2, vlen), 2, vlen);
  VUINT cos_sign = __riscv_vmerge
    (VMVU_VX (0, vlen), VMVU_VX (0x80000000u, vlen), cos_flip, vlen);
  cos_result = U_AS_F
    (__riscv_vxor (F_AS_U (cos_result), cos_sign, vlen));

  *out_sin = sin_result;
  *out_cos = cos_result;
}

#undef LMUL
#undef SIMDLEN
