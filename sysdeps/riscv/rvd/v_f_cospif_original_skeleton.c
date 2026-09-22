/* Single-precision source-original SLEEF cospif skeleton.

   This is a target-owned translation of SLEEF xcospif_u05/cospifk from the
   pinned source tree.  The source signed quotient correction, polynomial
   selection, double-float reconstruction order, quadrant sign operation, and
   late range/infinity handling are kept separate from the clean cospif
   implementation.  Only the RVV representation, LMUL/VL plumbing, and the
   glibc importer ABI are adapted.  No SLEEF source file is included,
   compiled, or linked.

   SLEEF's source build uses explicit fused multiply-add operations and
   -ffp-contract=off; the wrapper's build rule preserves that contract for
   this comparison implementation.

   Copyright (C) 2026 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public License as
   published by the Free Software Foundation; either version 2.1, or (at
   your option) any later version.

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

/* Reusable target-owned translation of SLEEF xcospif_u05.  */
static VFLOAT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_cospif_kernel_m, LMUL)
  (VFLOAT d, size_t vlen)
{
  /* SLEEF cospifk: u = 4*d and q is the even truncation with the source's
     signed correction.  The correction matters for negative arguments; the
     clean implementation deliberately uses a different abs-domain form.  */
  VFLOAT u = VFMUL_VF (d, 4.0f, vlen);
  VINT q = __riscv_vfcvt_rtz_x (u, vlen);
  VUINT q_bits = I_AS_U (q);
  VUINT q_sign = __riscv_vsrl (q_bits, BIT_WIDTH_F - 1, vlen);
  VUINT q_correction = __riscv_vxor
    (q_sign, VMVU_VX (1, vlen), vlen);
  q = __riscv_vand
    (__riscv_vadd (q, U_AS_I (q_correction), vlen),
     VMVI_VX (~1, vlen), vlen);

  VBOOL use_cos = __riscv_vmseq
    (__riscv_vand (q, 2, vlen), 0, vlen);

  VFLOAT q_float = __riscv_vfcvt_f (q, vlen);
  VFLOAT s = VFSUB_VV (u, q_float, vlen);
  VFLOAT t = s;
  s = VFMUL_VV (s, s, vlen);

  VFLOAT s2_hi, s2_lo;
  RVV_SLEEF_ORIG_DFMUL_VF_VF
    (t, t, s2_hi, s2_lo, vlen);

  /* Source vsel_vf_vo_f_f(o, cos_value, sin_value).  */
  VFLOAT p = VFMERGE_VFM
    (VFMV_VF (RVV_COSPIF_POLY_COS0_F, vlen), RVV_COSPIF_POLY_SIN0_F,
     use_cos, vlen);
  p = VFMADD_VV
    (p, s,
     VFMERGE_VFM (VFMV_VF (RVV_COSPIF_POLY_COS1_F, vlen),
                  RVV_COSPIF_POLY_SIN1_F, use_cos, vlen), vlen);
  p = VFMADD_VV
    (p, s,
     VFMERGE_VFM (VFMV_VF (RVV_COSPIF_POLY_COS2_F, vlen),
                  RVV_COSPIF_POLY_SIN2_F, use_cos, vlen), vlen);

  VFLOAT c_hi = VFMERGE_VFM
    (VFMV_VF (RVV_COSPIF_DF_FIRST_COS_HI_F, vlen),
     RVV_COSPIF_DF_FIRST_SIN_HI_F, use_cos, vlen);
  VFLOAT c_lo = VFMERGE_VFM
    (VFMV_VF (RVV_COSPIF_DF_FIRST_COS_LO_F, vlen),
     RVV_COSPIF_DF_FIRST_SIN_LO_F, use_cos, vlen);
  VFLOAT x_hi, x_lo;
  RVV_SLEEF_ORIG_DFADD2_VF_VF2
    (VFMUL_VV (p, s, vlen), c_hi, c_lo, x_hi, x_lo, vlen);

  VFLOAT tmp_hi, tmp_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (s2_hi, s2_lo, x_hi, x_lo, tmp_hi, tmp_lo, vlen);
  c_hi = VFMERGE_VFM
    (VFMV_VF (RVV_COSPIF_DF_SECOND_COS_HI_F, vlen),
     RVV_COSPIF_DF_SECOND_SIN_HI_F, use_cos, vlen);
  c_lo = VFMERGE_VFM
    (VFMV_VF (RVV_COSPIF_DF_SECOND_COS_LO_F, vlen),
     RVV_COSPIF_DF_SECOND_SIN_LO_F, use_cos, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF2
    (tmp_hi, tmp_lo, c_hi, c_lo, x_hi, x_lo, vlen);

  /* x *= (s2 or t,0), followed by the source cos-only +1 reconstruction.  */
  VFLOAT mul_hi = VMERGE_VVM_F
    (t, s2_hi, use_cos, vlen);
  VFLOAT mul_lo = VMERGE_VVM_F
    (VFMV_VF (0.0f, vlen), s2_lo, use_cos, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (x_hi, x_lo, mul_hi, mul_lo, tmp_hi, tmp_lo, vlen);
  x_hi = tmp_hi;
  x_lo = tmp_lo;

  VFLOAT add_hi, add_lo;
  RVV_SLEEF_ORIG_DFADD2_VF2_VF
    (x_hi, x_lo, VFMV_VF (1.0f, vlen), add_hi, add_lo, vlen);
  x_hi = VMERGE_VVM_F (x_hi, add_hi, use_cos, vlen);
  x_lo = VMERGE_VVM_F (x_lo, add_lo, use_cos, vlen);

  /* Source flips both double-float components when ((q + 2) & 4) == 4.  */
  VUINT q_plus2 = I_AS_U (__riscv_vadd (q, 2, vlen));
  VBOOL flip = __riscv_vmseq
    (__riscv_vand (q_plus2, 4, vlen), 4, vlen);
  VUINT sign_bit = VMVU_VX (0x80000000u, vlen);
  VUINT flip_bits = __riscv_vmerge
    (VMVU_VX (0, vlen), sign_bit, flip, vlen);
  x_hi = U_AS_F (__riscv_vxor (F_AS_U (x_hi), flip_bits, vlen));
  x_lo = U_AS_F (__riscv_vxor (F_AS_U (x_lo), flip_bits, vlen));
  VFLOAT result = VFADD_VV (x_hi, x_lo, vlen);

  /* xcospif_u05's late source-native range and infinity operations.  */
  VBOOL out_of_range = __riscv_vmfgt
    (__riscv_vfabs (d, vlen), RVV_COSPIF_TRIG_RANGE_MAX4_F,
                    vlen);
  result = VFMERGE_VFM (result, 1.0f, out_of_range, vlen);

  VBOOL input_inf = __riscv_vmfeq
    (__riscv_vfabs (d, vlen), fp_posInf_f, vlen);
  VUINT inf_bits = __riscv_vmerge
    (VMVU_VX (0, vlen), VMVU_VX (~0u, vlen), input_inf, vlen);
  result = U_AS_F (__riscv_vor (F_AS_U (result), inf_bits, vlen));

  return result;
}

#ifndef RVV_SLEEF_ORIG_COSPIF_HELPERS_ONLY
static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT d, size_t vlen)
{
  return RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_cospif_kernel_m, LMUL) (d, vlen);
}
#endif

#ifndef RVV_SLEEF_ORIG_COSPIF_HELPERS_ONLY
#undef LMUL
#undef SIMDLEN
#endif
