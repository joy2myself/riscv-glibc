/* Single-precision source-original SLEEF sinf skeleton.

   This is a target-owned translation of SLEEF xsinf_u1, which is exported as
   Sleef_sinf_u10 by the pinned source rename table.  The direct source
   reduction, rempif large-argument path, scalar sine polynomial, sign
   reconstruction, and late special-value behavior remain source-shaped.
   Only the RVV representation, LMUL/VL plumbing, target-owned data, and
   glibc importer ABI are adapted.  No SLEEF source file is included,
   compiled, or linked.

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
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT d, size_t vlen)
{
  VFLOAT input = d;

  /* Source xsinf_u1 starts with vrint(d * M_1_PI), followed by the source
     integer quotient conversion and an integer-to-float cast.  */
  VFLOAT rounded_product = RVV_LIBMVEC_IMPL_NAME
    (rvv_vrint_f32, LMUL)
    (VFMUL_VF (d, RVV_ONE_OVER_PI_F, vlen), vlen);
  VINT q = RVV_VFCVT_X_F_RNE (rounded_product, vlen);
  VFLOAT u = __riscv_vfcvt_f (q, vlen);

  VFLOAT v = VFMADD_VV
    (u, VFMV_VF (RVV_SLEEF_TRIG_NEG_PI_A2_F, vlen), d, vlen);
  VFLOAT s_hi, s_lo;
  RVV_SLEEF_ORIG_DFADD2_VF_VF
    (v, VFMUL_VF (u, RVV_SLEEF_TRIG_NEG_PI_B2_F, vlen), s_hi, s_lo, vlen);
  RVV_SLEEF_ORIG_DFADD_VF2_VF2_VF
    (s_hi, s_lo, VFMUL_VF (u, RVV_SLEEF_TRIG_NEG_PI_C2_F, vlen),
     s_hi, s_lo, vlen);
  VBOOL g = __riscv_vmflt
    (__riscv_vfabs (input, vlen), RVV_SLEEF_TRIG_RANGE_MAX2, vlen);

  if (__riscv_vcpop (g, vlen) != vlen)
    {
      VFLOAT rem_hi, rem_lo;
      VINT rem_q;
      RVV_LIBMVEC_IMPL_NAME (sinf_original_reduce_large, LMUL)
        (input, &rem_hi, &rem_lo, &rem_q, vlen);

      VINT q2 = __riscv_vand (rem_q, 3, vlen);
      VBOOL rem_positive = __riscv_vmfgt
        (rem_hi, fp_posZero_f, vlen);
      VINT q2_add = VMERGE_VXM_I
        (VMVI_VX (1, vlen), 2, rem_positive, vlen);
      q2 = __riscv_vadd (__riscv_vadd (q2, q2, vlen), q2_add, vlen);
      q2 = __riscv_vsra (q2, 2, vlen);

      VBOOL rem_q_odd = __riscv_vmseq
        (__riscv_vand (rem_q, 1, vlen), 1, vlen);
      VFLOAT corr_hi = __riscv_vfsgnjx
        (VFMV_VF ((-RVV_PIBY2_F), vlen), rem_hi, vlen);
      VFLOAT corr_lo = __riscv_vfsgnjx
        (VFMV_VF ((-RVV_PIBY2_LO_F), vlen), rem_hi, vlen);
      VFLOAT alt_hi, alt_lo;
      RVV_SLEEF_ORIG_DFADD2_VF2_VF2
        (rem_hi, rem_lo, corr_hi, corr_lo, alt_hi, alt_lo, vlen);
      VFLOAT selected_hi = VMERGE_VVM_F
        (rem_hi, alt_hi, rem_q_odd, vlen);
      VFLOAT selected_lo = VMERGE_VVM_F
        (rem_lo, alt_lo, rem_q_odd, vlen);
      RVV_SLEEF_ORIG_DFNORMALIZE_VF2
        (selected_hi, selected_lo, selected_hi, selected_lo, vlen);

      /* Source ORs all bits into the high part for Inf/NaN inputs after the
         rempif result is formed.  */
      VBOOL input_inf = __riscv_vmfeq
        (__riscv_vfabs (input, vlen), fp_posInf_f, vlen);
      VBOOL input_nan = __riscv_vmfne (input, input, vlen);
      VBOOL input_nonfinite = __riscv_vmor
        (input_inf, input_nan, vlen);
      VUINT poison = __riscv_vmerge
        (VMVU_VX (0, vlen), VMVU_VX (~0u, vlen), input_nonfinite, vlen);
      selected_hi = U_AS_F
        (__riscv_vor (F_AS_U (selected_hi), poison, vlen));

      VBOOL large = __riscv_vmnot (g, vlen);
      q = __riscv_vmerge (q, q2, large, vlen);
      s_hi = VMERGE_VVM_F (s_hi, selected_hi, large, vlen);
      s_lo = VMERGE_VVM_F (s_lo, selected_lo, large, vlen);
    }

  VFLOAT t_hi = s_hi;
  VFLOAT t_lo = s_lo;
  VFLOAT s2_hi, s2_lo;
  RVV_SLEEF_ORIG_DFSQU_VF2
    (s_hi, s_lo, s2_hi, s2_lo, vlen);

  /* This is the source's scalar polynomial on the high part of s^2.  */
  u = VFMV_VF (RVV_SLEEF_TRIG_SIN_POLY0_F, vlen);
  u = VFMADD_VV
    (u, s2_hi, VFMV_VF (RVV_SLEEF_TRIG_SIN_POLY1_F, vlen), vlen);
  u = VFMADD_VV
    (u, s2_hi, VFMV_VF (RVV_SLEEF_TRIG_SIN_POLY2_F, vlen), vlen);

  VFLOAT inner_hi, inner_lo;
  RVV_SLEEF_ORIG_DFADD_VF_VF
    (VFMV_VF (RVV_SLEEF_TRIG_SIN_POLY3_F, vlen),
     VFMUL_VV (u, s2_hi, vlen), inner_hi, inner_lo, vlen);
  VFLOAT prod_hi, prod_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (inner_hi, inner_lo, s2_hi, s2_lo, prod_hi, prod_lo, vlen);
  VFLOAT x_hi, x_lo;
  RVV_SLEEF_ORIG_DFADD_VF_VF2
    (VFMV_VF (1.0f, vlen), prod_hi, prod_lo, x_hi, x_lo, vlen);

  VFLOAT result;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2_SCALAR
    (t_hi, t_lo, x_hi, x_lo, result, vlen);

  VBOOL q_odd = __riscv_vmseq
    (__riscv_vand (q, 1, vlen), 1, vlen);
  VUINT sign = __riscv_vmerge
    (VMVU_VX (0, vlen), VMVU_VX (0x80000000u, vlen), q_odd, vlen);
  result = U_AS_F (__riscv_vxor (F_AS_U (result), sign, vlen));

  VBOOL negative_zero = __riscv_vmseq
    (F_AS_U (input), 0x80000000u, vlen);
  return VMERGE_VVM_F (result, input, negative_zero, vlen);
}

#ifndef RVV_SLEEF_ORIG_KEEP_LMUL
# undef LMUL
# undef SIMDLEN
#endif
