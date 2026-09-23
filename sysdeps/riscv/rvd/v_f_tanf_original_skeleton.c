/* Single-precision source-original SLEEF tanf skeleton.

   This is a target-owned translation of SLEEF xtanf_u1, which is exported as
   Sleef_tanf_u10 by the pinned source rename table.  The direct source
   reduction (2/pi), rempif large-argument path, double-float reconstruction,
   tangent polynomial, sign reconstruction, and late special-value
   behavior are kept in source order.  Only the RVV representation, LMUL/VL
   plumbing, target-owned double-float helpers, and glibc importer ABI are
   adapted.  No SLEEF source file is included, compiled, or linked.

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

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT d, size_t vlen)
{
  VINT q;
  VFLOAT u, v;
  VFLOAT s_hi, s_lo, t_hi, t_lo;
  VFLOAT x_hi, x_lo;

  /* xtanf_u1: source vrint(2*d/pi), then integer quotient conversion.  */
  u = VFMUL_VF (d, RVV_SLEEF_TRIG_TWO_OVER_PI, vlen);
  u = __riscv_vfcvt_f
    (RVV_VFCVT_X_F_RNE (u, vlen), vlen);
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

      /* xtanf_u1 performs late source-native Inf/NaN propagation by ORing
         the all-one mask into both double-float components.  */
      VBOOL input_inf = __riscv_vmfeq
        (__riscv_vfabs (d, vlen), fp_posInf_f, vlen);
      VBOOL input_nan = __riscv_vmfne (d, d, vlen);
      VBOOL input_nonfinite = __riscv_vmor (input_inf, input_nan, vlen);
      VUINT poison = __riscv_vmerge
        (VMVU_VX (0, vlen), VMVU_VX (~0u, vlen), input_nonfinite, vlen);
      rem_hi = U_AS_F
        (__riscv_vor (F_AS_U (rem_hi), poison, vlen));
      rem_lo = U_AS_F
        (__riscv_vor (F_AS_U (rem_lo), poison, vlen));

      VBOOL large = __riscv_vmnot (g, vlen);
      q = __riscv_vmerge (q, rem_q, large, vlen);
      s_hi = __riscv_vmerge (s_hi, rem_hi, large, vlen);
      s_lo = __riscv_vmerge (s_lo, rem_lo, large, vlen);
    }

  /* q odd: source flips both components of the reduced angle before saving
     t and evaluating the tangent polynomial.  */
  VBOOL odd = __riscv_vmseq
    (__riscv_vand (q, 1, vlen), 1, vlen);
  VUINT sign_mask = __riscv_vmerge
    (VMVU_VX (0, vlen), VMVU_VX (0x80000000u, vlen), odd, vlen);
  s_hi = U_AS_F
    (__riscv_vxor (F_AS_U (s_hi), sign_mask, vlen));
  s_lo = U_AS_F
    (__riscv_vxor (F_AS_U (s_lo), sign_mask, vlen));

  t_hi = s_hi;
  t_lo = s_lo;
  RVV_SLEEF_ORIG_DFSQU_VF2 (s_hi, s_lo, s_hi, s_lo, vlen);
  RVV_SLEEF_ORIG_DFNORMALIZE_VF2 (s_hi, s_lo, s_hi, s_lo, vlen);

  /* Source tangent polynomial, retaining its Horner order and literals.  */
  u = VFMV_VF (RVV_TANF_POLY0, vlen);
  u = VFMADD_VV
    (u, s_hi, VFMV_VF (RVV_TANF_POLY1, vlen), vlen);
  u = VFMADD_VV
    (u, s_hi, VFMV_VF (RVV_TANF_POLY2, vlen), vlen);
  u = VFMADD_VV
    (u, s_hi, VFMV_VF (RVV_TANF_POLY3, vlen), vlen);
  u = VFMADD_VV
    (u, s_hi, VFMV_VF (RVV_TANF_POLY4, vlen), vlen);

  RVV_SLEEF_ORIG_DFADD_VF_VF
    (VFMV_VF (RVV_TANF_DF_C0, vlen),
     VFMUL_VV (u, s_hi, vlen), x_hi, x_lo, vlen);

  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (s_hi, s_lo, x_hi, x_lo, x_hi, x_lo, vlen);
  RVV_SLEEF_ORIG_DFADD_VF_VF2
    (VFMV_VF (RVV_TANF_DF_C1, vlen),
     x_hi, x_lo, x_hi, x_lo, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (s_hi, s_lo, x_hi, x_lo, x_hi, x_lo, vlen);
  RVV_SLEEF_ORIG_DFADD_VF_VF2
    (VFMV_VF (1.0f, vlen), x_hi, x_lo, x_hi, x_lo, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (t_hi, t_lo, x_hi, x_lo, x_hi, x_lo, vlen);

  VFLOAT rec_hi, rec_lo;
  RVV_SLEEF_ORIG_DFREC_VF2_VF2
    (x_hi, x_lo, rec_hi, rec_lo, vlen);
  x_hi = __riscv_vmerge (x_hi, rec_hi, odd, vlen);
  x_lo = __riscv_vmerge (x_lo, rec_lo, odd, vlen);

  VFLOAT result = VFADD_VV (x_hi, x_lo, vlen);
  VBOOL negative_zero = __riscv_vmseq
    (F_AS_U (d), 0x80000000u, vlen);
  result = VMERGE_VVM_F (result, d, negative_zero, vlen);

  return result;
}

#undef LMUL
#undef SIMDLEN
