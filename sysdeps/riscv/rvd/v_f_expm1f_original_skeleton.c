/* Single-precision source-original SLEEF expm1f skeleton.

   This is a target-owned translation of SLEEF expk2f and xexpm1f.  The
   source double-float reduction, polynomial, reconstruction, vldexp2
   sequence, and late threshold selections are kept separate from the clean
   implementation.  Only the RVV representation, LMUL/VL plumbing, and
   glibc importer ABI are adapted; no SLEEF source file is included or linked.

   SLEEF's source file requests -ffp-contract=off.  The corresponding build
   flag is attached to this wrapper in the RVV mathvec Makefile.

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

/* SLEEF vpow2i_vf_vi2.  The source exponent range is bounded by expk2f's
   reduction, so the source bit construction is sufficient here.  */
static inline VFLOAT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_expm1f_vpow2i_m, LMUL)
  (VINT q, size_t vlen)
{
  VINT biased = __riscv_vadd (q, EXP_BIAS_F, vlen);
  return I_AS_F (__riscv_vsll (biased, MAN_LEN_F, vlen));
}

/* SLEEF vldexp2_vf_vf_vi2: split the exponent and multiply by both source
   power-of-two factors in source order.  */
static inline VFLOAT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_expm1f_vldexp2_m, LMUL)
  (VFLOAT d, VINT e, size_t vlen)
{
  VINT e_half = __riscv_vsra (e, 1, vlen);
  VINT e_rest = __riscv_vsub (e, e_half, vlen);
  VFLOAT p_half = RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_expm1f_vpow2i_m, LMUL) (e_half, vlen);
  VFLOAT p_rest = RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_expm1f_vpow2i_m, LMUL) (e_rest, vlen);
  return VFMUL_VV (VFMUL_VV (d, p_half, vlen), p_rest, vlen);
}

/* SLEEF expk2f, retaining its source reduction and reconstruction order.  */
static inline void
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_expm1f_expk2f_m, LMUL)
  (VFLOAT d_hi, VFLOAT d_lo, VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT d = VFADD_VV (d_hi, d_lo, vlen);
  VFLOAT u = VFMUL_VF (d, RVV_SLEEF_EXPM1F_R_LN2, vlen);
  VINT q = RVV_VFCVT_X_F_RNE (u, vlen);
  VFLOAT q_f = __riscv_vfcvt_f (q, vlen);

  VFLOAT s_hi, s_lo;
  VFLOAT q_l2u = VFMUL_VF (q_f, -RVV_SLEEF_EXPM1F_L2U, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF (d_hi, d_lo, q_l2u,
                                s_hi, s_lo, vlen);
  VFLOAT q_l2l = VFMUL_VF (q_f, -RVV_SLEEF_EXPM1F_L2L, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF (s_hi, s_lo, q_l2l,
                                s_hi, s_lo, vlen);

  u = VFMV_VF (RVV_SLEEF_EXPM1F_POLY0, vlen);
  u = VFMADD_VV (u, s_hi, VFMV_VF (RVV_SLEEF_EXPM1F_POLY1, vlen), vlen);
  u = VFMADD_VV (u, s_hi, VFMV_VF (RVV_SLEEF_EXPM1F_POLY2, vlen), vlen);
  u = VFMADD_VV (u, s_hi, VFMV_VF (RVV_SLEEF_EXPM1F_POLY3, vlen), vlen);

  VFLOAT t_hi, t_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF (s_hi, s_lo, u, t_hi, t_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF_S
    (t_hi, t_lo, RVV_SLEEF_EXPM1F_TAYLOR_C, t_hi, t_lo, vlen);

  VFLOAT tmp_hi, tmp_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (s_hi, s_lo, t_hi, t_lo, tmp_hi, tmp_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF_S
    (tmp_hi, tmp_lo, 0.5, t_hi, t_lo, vlen);

  VFLOAT s2_hi, s2_lo;
  RVV_SLEEF_ORIG_DFSQU_VF2 (s_hi, s_lo, s2_hi, s2_lo, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (s2_hi, s2_lo, t_hi, t_lo, tmp_hi, tmp_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF2
    (s_hi, s_lo, tmp_hi, tmp_lo, t_hi, t_lo, vlen);

  /* This is source dfadd_vf2_vf_vf2, not the compensated dfadd2 variant.  */
  RVV_SLEEF_ORIG_DFADD_VF_VF2
    (VFMV_VF (1.0f, vlen), t_hi, t_lo, t_hi, t_lo, vlen);

  t_hi = RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_expm1f_vldexp2_m, LMUL) (t_hi, q, vlen);
  t_lo = RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_expm1f_vldexp2_m, LMUL) (t_lo, q, vlen);

  VBOOL underflow = __riscv_vmflt (d_hi, RVV_SLEEF_EXPM1F_UNDERFLOW, vlen);
  t_hi = VFMERGE_VFM (t_hi, fp_posZero_f, underflow, vlen);
  t_lo = VFMERGE_VFM (t_lo, fp_posZero_f, underflow, vlen);

  *r_hi = t_hi;
  *r_lo = t_lo;
}

#ifndef RVV_SLEEF_ORIG_EXPM1F_HELPERS_ONLY
static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT a, size_t vlen)
{
  VFLOAT exp_hi, exp_lo;
  RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_expm1f_expk2f_m, LMUL)
    (a, VFMV_VF (0.0f, vlen), &exp_hi, &exp_lo, vlen);

  VFLOAT d_hi, d_lo;
  RVV_SLEEF_ORIG_DFADD2_VF2_VF_S
    (exp_hi, exp_lo, -1.0f, d_hi, d_lo, vlen);
  VFLOAT result = VFADD_VV (d_hi, d_lo, vlen);

  VBOOL overflow = __riscv_vmfgt (a, RVV_SLEEF_EXPM1F_MAX, vlen);
  result = VFMERGE_VFM (result, fp_posInf_f, overflow, vlen);

  VBOOL underflow = __riscv_vmflt (a, RVV_SLEEF_EXPM1F_MIN, vlen);
  result = VFMERGE_VFM (result, fp_negOne_f, underflow, vlen);

  VBOOL negzero = __riscv_vmseq
    (F_AS_U (a), VMVU_VX (0x80000000u, vlen), vlen);
  result = VFMERGE_VFM (result, fp_negZero_f, negzero, vlen);

  return result;
}
#endif

#ifndef RVV_SLEEF_ORIG_EXPM1F_HELPERS_ONLY
#undef LMUL
#undef SIMDLEN
#endif
