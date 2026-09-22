/* Single-precision source-original SLEEF cosf skeleton.

   This is a target-owned translation of SLEEF xcosf_u1, which is exported as
   Sleef_cosf_u10 by the pinned source rename table.  The direct source
   reduction (1/pi), rempif large-argument path, double-float reconstruction,
   dedicated sine polynomial, sign reconstruction, and late special-value
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

/* SLEEF vilogb2k for the normal/large finite range used by rempif.  */
static inline VINT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_cosf_vilogb2k, LMUL) (VFLOAT a,
                                                             size_t vlen)
{
  VUINT bits = F_AS_U (a);
  VUINT exponent = __riscv_vsrl (bits, MAN_LEN_F, vlen);
  exponent = __riscv_vand (exponent, 0xff, vlen);
  return __riscv_vsub (U_AS_I (exponent), EXP_BIAS_F, vlen);
}

/* SLEEF vldexp3_vf_vf_vi2: exponent-field adjustment for normalized values.  */
static inline VFLOAT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_cosf_vldexp3, LMUL) (VFLOAT a, VINT q,
                                                            size_t vlen)
{
  return I_AS_F (__riscv_vadd
                 (F_AS_I (a), __riscv_vsll (q, MAN_LEN_F, vlen), vlen));
}

/* SLEEF vorsign_vf_vf_vf: preserve the magnitude bits and OR in the sign
   bit of the second operand.  vfsgnj is not equivalent for a source-shaped
   port because it overwrites all sign semantics of the first result.  */
static inline VFLOAT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_cosf_vorsign, LMUL) (VFLOAT x, VFLOAT y,
                                                            size_t vlen)
{
  VUINT sign = __riscv_vand
    (F_AS_U (y), 0x80000000u, vlen);
  return U_AS_F (__riscv_vor (F_AS_U (x), sign, vlen));
}

/* SLEEF vrint_vf_vf: explicit RNE conversion, independent of the enclosing
   thread's current rounding mode.  */
static inline VFLOAT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_cosf_vrint, LMUL) (VFLOAT x,
                                                          size_t vlen)
{
  VINT rounded = RVV_VFCVT_X_F_RNE (x, vlen);
  return __riscv_vfcvt_f (rounded, vlen);
}

/* Source rempisubf in the RVV helper uses this non-FULL_FP_ROUNDING branch.  */
static inline VFLOAT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_cosf_rempisubf, LMUL)
  (VFLOAT x, VINT *q_out, size_t vlen)
{
  /* SLEEF vmulsign(2^23, x): preserve the source sign-XOR primitive.  */
  VFLOAT c = __riscv_vfsgnjx
    (VFMV_VF (RVV_SLEEF_TRIG_REMPI_TWO23, vlen), x, vlen);
  VFLOAT x4 = VFMUL_VF (x, RVV_SLEEF_TRIG_REMPI_FOUR, vlen);
  VFLOAT abs_x4 = __riscv_vfabs (x4, vlen);
  VFLOAT abs_x = __riscv_vfabs (x, vlen);
  /* The source predicate is strict `> 2^23`; equality must remain on the
     round-with-a-sign-bit path rather than taking the direct path.  */
  VBOOL x4_large = __riscv_vmfgt
    (abs_x4, RVV_SLEEF_TRIG_REMPI_TWO23, vlen);
  VBOOL x_large = __riscv_vmfgt
    (abs_x, RVV_SLEEF_TRIG_REMPI_TWO23, vlen);

  VFLOAT rint4x = VFMADD_VV
    (VFMV_VF (RVV_SLEEF_TRIG_REMPI_FOUR, vlen), x, c, vlen);
  rint4x = VFSUB_VV (rint4x, c, vlen);
  rint4x = RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_cosf_vorsign, LMUL)
    (rint4x, x, vlen);
  rint4x = VMERGE_VVM_F (rint4x, x4, x4_large, vlen);

  VFLOAT rintx = VFADD_VV (x, c, vlen);
  rintx = VFSUB_VV (rintx, c, vlen);
  rintx = RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_cosf_vorsign, LMUL)
    (rintx, x, vlen);
  rintx = VMERGE_VVM_F (rintx, x, x_large, vlen);

  VFLOAT rem = VFMADD_VV
    (VFMV_VF (RVV_SLEEF_TRIG_REMPI_NEG_QUARTER, vlen), rint4x, x,
              vlen);
  *q_out = __riscv_vfcvt_rtz_x
    (VFMADD_VV (VFMV_VF (RVV_SLEEF_TRIG_REMPI_NEG_FOUR, vlen), rintx,
                rint4x, vlen), vlen);
  return rem;
}

/* Target-owned translation of SLEEF rempif.  Its output is the normalized
   double-float remainder and the quotient used by xcosf_u1.  */
static inline void
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_cosf_rempif, LMUL)
  (VFLOAT a, VFLOAT *out_hi, VFLOAT *out_lo, VINT *out_q, size_t vlen)
{
  const struct v_sleef_trig_data *trig_data =
    ptr_barrier (&__v_sleef_trig_data);

  VINT ex = RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_cosf_vilogb2k, LMUL) (a, vlen);
  ex = __riscv_vsub (ex, 25, vlen);

  VBOOL scale_mask = __riscv_vmsgt (ex, 90 - 25, vlen);
  VINT scale = VMERGE_VXM_I (VMVI_VX (0, vlen), -64, scale_mask, vlen);
  a = RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_cosf_vldexp3, LMUL) (a, scale, vlen);

  VBOOL ex_negative = __riscv_vmslt (ex, 0, vlen);
  ex = VMERGE_VXM_I (ex, 0, ex_negative, vlen);
  /* SLEEF shifts by two before vgather, whose helper multiplies by
     sizeof(float).  This translation calls vluxei32 directly, so the
     equivalent byte offset is ex << 4.  */
  ex = __riscv_vsll (ex, 4, vlen);
  VUINT ex_u = I_AS_U (ex);

  VFLOAT tab0 = RVV_VLUXEI32_F (trig_data->rempitabsp, ex_u, vlen);
  VFLOAT tab1 = RVV_VLUXEI32_F (trig_data->rempitabsp + 1, ex_u, vlen);
  VFLOAT tab2 = RVV_VLUXEI32_F (trig_data->rempitabsp + 2, ex_u, vlen);
  VFLOAT tab3 = RVV_VLUXEI32_F (trig_data->rempitabsp + 3, ex_u, vlen);

  VFLOAT x_hi, x_lo;
  RVV_SLEEF_ORIG_DFMUL_VF_VF (a, tab0, x_hi, x_lo, vlen);
  VINT q;
  x_hi = RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_cosf_rempisubf, LMUL)
    (x_hi, &q, vlen);
  RVV_SLEEF_ORIG_DFNORMALIZE_VF2 (x_hi, x_lo, x_hi, x_lo, vlen);

  VFLOAT y_hi, y_lo;
  RVV_SLEEF_ORIG_DFMUL_VF_VF (a, tab1, y_hi, y_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF2
    (x_hi, x_lo, y_hi, y_lo, x_hi, x_lo, vlen);
  VINT q_part;
  x_hi = RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_cosf_rempisubf, LMUL)
    (x_hi, &q_part, vlen);
  q = __riscv_vadd (q, q_part, vlen);
  RVV_SLEEF_ORIG_DFNORMALIZE_VF2 (x_hi, x_lo, x_hi, x_lo, vlen);

  y_hi = tab2;
  y_lo = tab3;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (y_hi, y_lo, a, y_hi, y_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF2
    (x_hi, x_lo, y_hi, y_lo, x_hi, x_lo, vlen);
  RVV_SLEEF_ORIG_DFNORMALIZE_VF2 (x_hi, x_lo, x_hi, x_lo, vlen);

  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (x_hi, x_lo,
     VFMV_VF (RVV_SLEEF_TRIG_REMPI_TWO_PI_HI, vlen),
     VFMV_VF (RVV_SLEEF_TRIG_REMPI_TWO_PI_LO, vlen),
     x_hi, x_lo, vlen);

  VBOOL a_small = __riscv_vmflt
    (__riscv_vfabs (a, vlen), RVV_SLEEF_TRIG_REMPI_SMALL, vlen);
  x_hi = VMERGE_VVM_F (x_hi, a, a_small, vlen);
  x_lo = VFMERGE_VFM (x_lo, fp_posZero_f, a_small, vlen);

  *out_hi = x_hi;
  *out_lo = x_lo;
  *out_q = q;
}

#if !RVV_MATH_IS_COSF_HELPERS_ONLY (V_MATH_FUNC)
static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT d, size_t vlen)
{
  /* xcosf_u1: round (d/pi - 1/2), then form dq = 2*n + 1.  */
  VFLOAT reduced = VFMADD_VV
    (d, VFMV_VF (RVV_SLEEF_TRIG_ONE_OVER_PI, vlen),
     VFMV_VF (RVV_SLEEF_COSF_REDUCE_NEG_HALF, vlen), vlen);
  VFLOAT u = RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_cosf_vrint, LMUL)
    (reduced, vlen);
  VFLOAT dq = VFMADD_VV
    (u, VFMV_VF (RVV_SLEEF_COSF_REDUCE_TWO, vlen),
     VFMV_VF (fp_posOne_f, vlen), vlen);
  VINT q = RVV_VFCVT_X_F_RNE (dq, vlen);

  VFLOAT s_hi, s_lo;
  RVV_SLEEF_ORIG_DFADD2_VF_VF
    (d, VFMUL_VF (dq, RVV_SLEEF_COSF_NEG_PI_A2_HALF, vlen),
     s_hi, s_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF
    (s_hi, s_lo,
     VFMUL_VF (dq, RVV_SLEEF_COSF_NEG_PI_B2_HALF, vlen),
     s_hi, s_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF
    (s_hi, s_lo,
     VFMUL_VF (dq, RVV_SLEEF_COSF_NEG_PI_C2_HALF, vlen),
     s_hi, s_lo, vlen);

  VBOOL g = __riscv_vmflt
    (__riscv_vfabs (d, vlen), RVV_SLEEF_TRIG_RANGE_MAX2, vlen);

  if (__riscv_vcpop (g, vlen) != vlen)
    {
      VFLOAT rem_hi, rem_lo;
      VINT rem_q;
      RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_cosf_rempif, LMUL)
        (d, &rem_hi, &rem_lo, &rem_q, vlen);

      VINT q2 = __riscv_vand (rem_q, 3, vlen);
      VBOOL rem_positive = __riscv_vmfgt (rem_hi, fp_posZero_f, vlen);
      VINT q2_add = VMERGE_VXM_I
        (VMVI_VX (7, vlen), 8, rem_positive, vlen);
      q2 = __riscv_vadd
        (__riscv_vadd (q2, q2, vlen), q2_add, vlen);
      q2 = __riscv_vsra (q2, 1, vlen);

      VBOOL even = __riscv_vmseq
        (__riscv_vand (rem_q, 1, vlen), 0, vlen);
      VFLOAT y = VFMERGE_VFM
        (VFMV_VF (fp_negOne_f, vlen), fp_posZero_f, rem_positive, vlen);
      VFLOAT corr_hi = __riscv_vfsgnjx
        (VFMV_VF (RVV_SLEEF_COSF_REMPI_HALF_PI_HI, vlen), y, vlen);
      VFLOAT corr_lo = __riscv_vfsgnjx
        (VFMV_VF (RVV_SLEEF_COSF_REMPI_HALF_PI_LO, vlen), y, vlen);

      VFLOAT alt_hi, alt_lo;
      RVV_SLEEF_ORIG_DFADD2_VF2_VF2
        (rem_hi, rem_lo, corr_hi, corr_lo, alt_hi, alt_lo, vlen);
      VFLOAT t_hi = VMERGE_VVM_F (rem_hi, alt_hi, even, vlen);
      VFLOAT t_lo = VMERGE_VVM_F (rem_lo, alt_lo, even, vlen);
      RVV_SLEEF_ORIG_DFNORMALIZE_VF2
        (t_hi, t_lo, t_hi, t_lo, vlen);

      /* This is SLEEF's late non-finite bit propagation, not clean EH
         quarantine.  It deliberately stays after rempif.  */
      VBOOL input_inf = __riscv_vmfeq
        (__riscv_vfabs (d, vlen), fp_posInf_f, vlen);
      VBOOL input_nan = __riscv_vmfne (d, d, vlen);
      VBOOL input_nonfinite = __riscv_vmor
        (input_inf, input_nan, vlen);
      VUINT poison = __riscv_vmerge
        (VMVU_VX (0, vlen), VMVU_VX (~0u, vlen),
         input_nonfinite, vlen);
      t_hi = U_AS_F
        (__riscv_vor (F_AS_U (t_hi), poison, vlen));

      /* RVV vmerge selects its second operand when the mask bit is one.
         `g` denotes the small direct-reduction path, whereas rempif is
         required for the complementary large-input lanes.  Keep the
         source's late selection, but select the rempif result with the
         inverted mask.  */
      VBOOL large = __riscv_vmnot (g, vlen);
      q = __riscv_vmerge (q, q2, large, vlen);
      s_hi = __riscv_vmerge (s_hi, t_hi, large, vlen);
      s_lo = __riscv_vmerge (s_lo, t_lo, large, vlen);
    }

  /* xcosf_u1 polynomial and final quadrant sign reconstruction.  */
  VFLOAT t_hi = s_hi;
  VFLOAT t_lo = s_lo;
  RVV_SLEEF_ORIG_DFSQU_VF2
    (s_hi, s_lo, s_hi, s_lo, vlen);

  VFLOAT poly = VFMV_VF (RVV_SLEEF_COSF_POLY_0, vlen);
  poly = VFMADD_VV
    (poly, s_hi,
     VFMV_VF (RVV_SLEEF_COSF_POLY_1, vlen), vlen);
  poly = VFMADD_VV
    (poly, s_hi,
     VFMV_VF (RVV_SLEEF_COSF_POLY_2, vlen), vlen);

  VFLOAT inner_hi, inner_lo;
  RVV_SLEEF_ORIG_DFADD_VF_VF
    (VFMV_VF (RVV_SLEEF_COSF_POLY_3, vlen),
     VFMUL_VV (poly, s_hi, vlen), inner_hi, inner_lo, vlen);

  VFLOAT product_hi, product_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (inner_hi, inner_lo, s_hi, s_lo,
     product_hi, product_lo, vlen);

  VFLOAT x_hi, x_lo;
  RVV_SLEEF_ORIG_DFADD_VF_VF2
    (VFMV_VF (fp_posOne_f, vlen), product_hi, product_lo,
     x_hi, x_lo, vlen);

  VFLOAT result;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2_SCALAR
    (t_hi, t_lo, x_hi, x_lo, result, vlen);

  VBOOL flip = __riscv_vmseq
    (__riscv_vand (q, 2, vlen), 0, vlen);
  VUINT sign_mask = __riscv_vmerge
    (VMVU_VX (0, vlen), VMVU_VX (0x80000000u, vlen),
     flip, vlen);
  result = U_AS_F
    (__riscv_vxor (F_AS_U (result), sign_mask, vlen));

  return result;
}
#endif

#if !RVV_MATH_IS_COSF_HELPERS_ONLY (V_MATH_FUNC)
#undef LMUL
#undef SIMDLEN
#endif
