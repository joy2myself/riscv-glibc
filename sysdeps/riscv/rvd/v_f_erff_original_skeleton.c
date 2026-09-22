/* Single-precision source-original SLEEF erff skeleton.

   This is a target-owned translation of SLEEF xerff_u1 from the pinned
   source tree.  The source Estrin POLY6, double-float polynomial helpers,
   expkf path, late special handling, and source operation order are kept
   separate from the clean erff implementation.  Only the RVV
   representation, LMUL/VL plumbing, target-owned helpers, and glibc
   importer ABI are adapted.  No SLEEF source file is included, compiled, or
   linked.

   SLEEF's source build requires -ffp-contract=off; the Makefile comparison
   rule preserves that contract for this translation unit.  */

#include <vecmath_config.h>
#include <v_math.h>

#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

/* Source vpow2i_vf_vi2.  The source erff expkf range stays within the
   exponent-field range used by this direct bit construction.  */
static inline VFLOAT
RVV_SLEEF_ORIG_ERFF_HELPER (vpow2i) (VINT e, size_t vlen)
{
  VINT biased = __riscv_vadd (e, EXP_BIAS_F, vlen);
  return I_AS_F (__riscv_vsll (biased, MAN_LEN_F, vlen));
}

/* Source vldexp_vf_vf_vi2, including its four-multiply split scaling.  This
   is deliberately not FAST_LDEXP: the latter is a clean-only adaptation.  */
static inline VFLOAT
RVV_SLEEF_ORIG_ERFF_HELPER (vldexp) (VFLOAT x, VINT q, size_t vlen)
{
  VINT m = __riscv_vsra (q, 31, vlen);
  VINT tmp = __riscv_vadd (m, q, vlen);
  tmp = __riscv_vsra (tmp, 6, vlen);
  tmp = __riscv_vsub (tmp, m, vlen);
  m = __riscv_vsll (tmp, 4, vlen);
  q = __riscv_vsub (q, __riscv_vsll (m, 2, vlen), vlen);

  m = __riscv_vadd (m, EXP_BIAS_F, vlen);
  VBOOL m_positive = __riscv_vmsgt (m, 0, vlen);
  m = __riscv_vmerge (VMVI_VX (0, vlen), m, m_positive, vlen);
  VBOOL m_too_large = __riscv_vmsgt (m, 0xff, vlen);
  m = __riscv_vmerge (m, 0xff, m_too_large, vlen);

  VFLOAT scale = I_AS_F (__riscv_vsll (m, MAN_LEN_F, vlen));
  x = VFMUL_VV (x, scale, vlen);
  x = VFMUL_VV (x, scale, vlen);
  x = VFMUL_VV (x, scale, vlen);
  x = VFMUL_VV (x, scale, vlen);

  scale = I_AS_F (__riscv_vsll
                  (__riscv_vadd (q, EXP_BIAS_F, vlen), MAN_LEN_F, vlen));
  return VFMUL_VV (x, scale, vlen);
}

/* Source dfmla_vf2_vf_vf2_vf2(x,y,z) = dfadd(z,dfmul(y,x)), using the
   non-compensated source dfadd rather than the clean DFADD2 variant.  */
static inline void
RVV_SLEEF_ORIG_ERFF_HELPER (poly2df) (VFLOAT x, VFLOAT c1,
                                      VFLOAT c0_hi, VFLOAT c0_lo,
                                      VFLOAT *r_hi, VFLOAT *r_lo,
                                      size_t vlen)
{
  VFLOAT t_hi, t_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (c1, VFMV_VF (0.0f, vlen), x, t_hi, t_lo, vlen);
  RVV_SLEEF_ORIG_DFADD_VF2_VF2
    (c0_hi, c0_lo, t_hi, t_lo, *r_hi, *r_lo, vlen);
}

static inline void
RVV_SLEEF_ORIG_ERFF_HELPER (poly2df_b) (VFLOAT x, VFLOAT c1_hi,
                                        VFLOAT c1_lo, VFLOAT c0_hi,
                                        VFLOAT c0_lo, VFLOAT *r_hi,
                                        VFLOAT *r_lo, size_t vlen)
{
  VFLOAT t_hi, t_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (c1_hi, c1_lo, x, t_hi, t_lo, vlen);
  RVV_SLEEF_ORIG_DFADD_VF2_VF2
    (c0_hi, c0_lo, t_hi, t_lo, *r_hi, *r_lo, vlen);
}

static inline void
RVV_SLEEF_ORIG_ERFF_HELPER (poly4df) (VFLOAT x, VFLOAT c3,
                                      VFLOAT c2_hi, VFLOAT c2_lo,
                                      VFLOAT c1_hi, VFLOAT c1_lo,
                                      VFLOAT c0_hi, VFLOAT c0_lo,
                                      VFLOAT *r_hi, VFLOAT *r_lo,
                                      size_t vlen)
{
  VFLOAT x2 = VFMUL_VV (x, x, vlen);
  VFLOAT p1_hi, p1_lo, p0_hi, p0_lo;
  RVV_SLEEF_ORIG_ERFF_HELPER (poly2df)
    (x, c3, c2_hi, c2_lo, &p1_hi, &p1_lo, vlen);
  RVV_SLEEF_ORIG_ERFF_HELPER (poly2df_b)
    (x, c1_hi, c1_lo, c0_hi, c0_lo, &p0_hi, &p0_lo, vlen);

  VFLOAT t_hi, t_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (p1_hi, p1_lo, x2, t_hi, t_lo, vlen);
  RVV_SLEEF_ORIG_DFADD_VF2_VF2
    (p0_hi, p0_lo, t_hi, t_lo, *r_hi, *r_lo, vlen);
}

/* Target-owned translation of POLY6/POLY4 from common/estrin.h.  The first
   operand of VFMADD_VV is the accumulated polynomial, matching SLEEF's MLA
   association.  */
static inline VFLOAT
RVV_SLEEF_ORIG_ERFF_HELPER (poly6) (VFLOAT x, VFLOAT x2, VFLOAT x4,
                                    VFLOAT c5, VFLOAT c4, VFLOAT c3,
                                    VFLOAT c2, VFLOAT c1, VFLOAT c0,
                                    size_t vlen)
{
  VFLOAT poly2 = VFMADD_VV (x, c5, c4, vlen);
  VFLOAT poly4_left = VFMADD_VV (x, c3, c2, vlen);
  VFLOAT poly4_right = VFMADD_VV (x, c1, c0, vlen);
  VFLOAT poly4 = VFMADD_VV (x2, poly4_left, poly4_right, vlen);
  return VFMADD_VV (x4, poly2, poly4, vlen);
}

/* Source expkf(vfloat2), retaining source DFADD/DFADD2 choices, normalize,
   high-plus-low reconstruction, vldexp scaling, and late underflow clear.  */
static inline VFLOAT
RVV_SLEEF_ORIG_ERFF_HELPER (expkf) (VFLOAT d_hi, VFLOAT d_lo,
                                    size_t vlen)
{
  VFLOAT u = VFMUL_VF (VFADD_VV (d_hi, d_lo, vlen),
                       RVV_ERFF_R_LN2_F, vlen);
  VINT q = RVV_SLEEF_ORIG_ERFF_VFCVT_X_RNE (u, vlen);
  VFLOAT q_f = __riscv_vfcvt_f (q, vlen);

  VFLOAT s_hi, s_lo;
  VFLOAT q_l2u = VFMUL_VF (q_f, -RVV_ERFF_L2U_F, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF
    (d_hi, d_lo, q_l2u, s_hi, s_lo, vlen);
  VFLOAT q_l2l = VFMUL_VF (q_f, -RVV_ERFF_L2L_F, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF
    (s_hi, s_lo, q_l2l, s_hi, s_lo, vlen);
  RVV_SLEEF_ORIG_DFNORMALIZE_VF2 (s_hi, s_lo, s_hi, s_lo, vlen);

  u = VFMV_VF (RVV_ERFF_EXP_POLY_0_F, vlen);
  u = VFMADD_VV (u, s_hi,
                 VFMV_VF (RVV_ERFF_EXP_POLY_1_F, vlen), vlen);
  u = VFMADD_VV (u, s_hi,
                 VFMV_VF (RVV_ERFF_EXP_POLY_2_F, vlen), vlen);
  u = VFMADD_VV (u, s_hi,
                 VFMV_VF (RVV_ERFF_EXP_POLY_3_F, vlen), vlen);
  u = VFMADD_VV (u, s_hi,
                 VFMV_VF (RVV_ERFF_EXP_POLY_4_F, vlen), vlen);

  VFLOAT s2_hi, s2_lo;
  RVV_SLEEF_ORIG_DFSQU_VF2 (s_hi, s_lo, s2_hi, s2_lo, vlen);
  VFLOAT t_hi, t_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (s2_hi, s2_lo, u, t_hi, t_lo, vlen);
  RVV_SLEEF_ORIG_DFADD_VF2_VF2
    (s_hi, s_lo, t_hi, t_lo, t_hi, t_lo, vlen);
  RVV_SLEEF_ORIG_DFADD_VF_VF2
    (VFMV_VF (1.0f, vlen), t_hi, t_lo, t_hi, t_lo, vlen);

  u = VFADD_VV (t_hi, t_lo, vlen);
  u = RVV_SLEEF_ORIG_ERFF_HELPER (vldexp) (u, q, vlen);

  VBOOL underflow = __riscv_vmflt (d_hi, RVV_ERFF_EXP_UNDERFLOW_F, vlen);
  return VFMERGE_VFM (u, fp_posZero_f, underflow, vlen);
}

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT a, size_t vlen)
{
  VFLOAT x = __riscv_vfabs (a, vlen);
  VFLOAT x2 = VFMUL_VV (x, x, vlen);
  VFLOAT x4 = VFMUL_VV (x2, x2, vlen);
  VBOOL o25 = __riscv_vmfle (x, RVV_ERFF_RANGE_SPLIT_F, vlen);

  VFLOAT t2_hi, t2_lo;
  if (__riscv_vcpop (o25, vlen) == vlen)
    {
      VFLOAT t = RVV_SLEEF_ORIG_ERFF_HELPER (poly6)
        (x, x2, x4,
         VFMV_VF (RVV_ERFF_POLY25_0_F, vlen),
         VFMV_VF (RVV_ERFF_POLY25_1_F, vlen),
         VFMV_VF (RVV_ERFF_POLY25_2_F, vlen),
         VFMV_VF (RVV_ERFF_POLY25_3_F, vlen),
         VFMV_VF (RVV_ERFF_POLY25_4_F, vlen),
         VFMV_VF (RVV_ERFF_POLY25_5_F, vlen), vlen);

      RVV_SLEEF_ORIG_ERFF_HELPER (poly4df)
        (x, t,
         VFMV_VF (RVV_ERFF_DF25_0_0_F, vlen),
         VFMV_VF (RVV_ERFF_DF25_0_1_F, vlen),
         VFMV_VF (RVV_ERFF_DF25_1_0_F, vlen),
         VFMV_VF (RVV_ERFF_DF25_1_1_F, vlen),
         VFMV_VF (RVV_ERFF_DF25_2_0_F, vlen),
         VFMV_VF (RVV_ERFF_DF25_2_1_F, vlen),
         &t2_hi, &t2_lo, vlen);
      RVV_SLEEF_ORIG_DFMUL_VF2_VF
        (t2_hi, t2_lo, x, t2_hi, t2_lo, vlen);
      RVV_SLEEF_ORIG_DFADD_VF_VF2
        (VFMV_VF (1.0f, vlen), t2_hi, t2_lo, t2_hi, t2_lo, vlen);
      RVV_SLEEF_ORIG_DFSQU_VF2 (t2_hi, t2_lo, t2_hi, t2_lo, vlen);
      RVV_SLEEF_ORIG_DFSQU_VF2 (t2_hi, t2_lo, t2_hi, t2_lo, vlen);
      RVV_SLEEF_ORIG_DFSQU_VF2 (t2_hi, t2_lo, t2_hi, t2_lo, vlen);
      RVV_SLEEF_ORIG_DFSQU_VF2 (t2_hi, t2_lo, t2_hi, t2_lo, vlen);
      RVV_SLEEF_ORIG_DFREC_VF2_VF2
        (t2_hi, t2_lo, t2_hi, t2_lo, vlen);
    }
  else
    {
      VFLOAT p5 = VFMERGE_VFM (VFMV_VF (RVV_ERFF_POLY_GT_0_F, vlen),
                               RVV_ERFF_POLY25_0_F, o25, vlen);
      VFLOAT p4 = VFMERGE_VFM (VFMV_VF (RVV_ERFF_POLY_GT_1_F, vlen),
                               RVV_ERFF_POLY25_1_F, o25, vlen);
      VFLOAT p3 = VFMERGE_VFM (VFMV_VF (RVV_ERFF_POLY_GT_2_F, vlen),
                               RVV_ERFF_POLY25_2_F, o25, vlen);
      VFLOAT p2 = VFMERGE_VFM (VFMV_VF (RVV_ERFF_POLY_GT_3_F, vlen),
                               RVV_ERFF_POLY25_3_F, o25, vlen);
      VFLOAT p1 = VFMERGE_VFM (VFMV_VF (RVV_ERFF_POLY_GT_4_F, vlen),
                               RVV_ERFF_POLY25_4_F, o25, vlen);
      VFLOAT p0 = VFMERGE_VFM (VFMV_VF (RVV_ERFF_POLY_GT_5_F, vlen),
                               RVV_ERFF_POLY25_5_F, o25, vlen);
      VFLOAT t = RVV_SLEEF_ORIG_ERFF_HELPER (poly6)
        (x, x2, x4, p5, p4, p3, p2, p1, p0, vlen);

      VFLOAT c2_hi = VFMERGE_VFM (VFMV_VF (RVV_ERFF_DF_GT_0_0_F, vlen),
                                  RVV_ERFF_DF25_0_0_F, o25, vlen);
      VFLOAT c2_lo = VFMERGE_VFM (VFMV_VF (RVV_ERFF_DF_GT_0_1_F, vlen),
                                  RVV_ERFF_DF25_0_1_F, o25, vlen);
      VFLOAT c1_hi = VFMERGE_VFM (VFMV_VF (RVV_ERFF_DF_GT_1_0_F, vlen),
                                  RVV_ERFF_DF25_1_0_F, o25, vlen);
      VFLOAT c1_lo = VFMERGE_VFM (VFMV_VF (RVV_ERFF_DF_GT_1_1_F, vlen),
                                  RVV_ERFF_DF25_1_1_F, o25, vlen);
      VFLOAT c0_hi = VFMERGE_VFM (VFMV_VF (RVV_ERFF_DF_GT_2_0_F, vlen),
                                  RVV_ERFF_DF25_2_0_F, o25, vlen);
      VFLOAT c0_lo = VFMERGE_VFM (VFMV_VF (RVV_ERFF_DF_GT_2_1_F, vlen),
                                  RVV_ERFF_DF25_2_1_F, o25, vlen);
      RVV_SLEEF_ORIG_ERFF_HELPER (poly4df)
        (x, t, c2_hi, c2_lo, c1_hi, c1_lo, c0_hi, c0_lo,
         &t2_hi, &t2_lo, vlen);
      RVV_SLEEF_ORIG_DFMUL_VF2_VF
        (t2_hi, t2_lo, x, t2_hi, t2_lo, vlen);

      VFLOAT s2_hi, s2_lo;
      RVV_SLEEF_ORIG_DFADD_VF_VF2
        (VFMV_VF (1.0f, vlen), t2_hi, t2_lo, s2_hi, s2_lo, vlen);
      RVV_SLEEF_ORIG_DFSQU_VF2 (s2_hi, s2_lo, s2_hi, s2_lo, vlen);
      RVV_SLEEF_ORIG_DFSQU_VF2 (s2_hi, s2_lo, s2_hi, s2_lo, vlen);
      RVV_SLEEF_ORIG_DFSQU_VF2 (s2_hi, s2_lo, s2_hi, s2_lo, vlen);
      RVV_SLEEF_ORIG_DFSQU_VF2 (s2_hi, s2_lo, s2_hi, s2_lo, vlen);
      RVV_SLEEF_ORIG_DFREC_VF2_VF2
        (s2_hi, s2_lo, s2_hi, s2_lo, vlen);

      VFLOAT exp_hi = RVV_SLEEF_ORIG_ERFF_HELPER
        (expkf) (t2_hi, t2_lo, vlen);
      t2_hi = VMERGE_VVM_F (exp_hi, s2_hi, o25, vlen);
      t2_lo = VMERGE_VVM_F (VFMV_VF (0.0f, vlen), s2_lo, o25, vlen);
    }

  RVV_SLEEF_ORIG_DFADD2_VF2_VF_S
    (t2_hi, t2_lo, -1.0f, t2_hi, t2_lo, vlen);

  VBOOL very_small = __riscv_vmflt (x, RVV_ERFF_TINY_THRESHOLD_F, vlen);
  VFLOAT small_hi, small_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (VFMV_VF (RVV_ERFF_TWO_OVER_SQRT_PI_F, vlen),
     VFMV_VF (RVV_ERFF_TWO_OVER_SQRT_PI_LO_F, vlen),
     x, small_hi, small_lo, vlen);
  t2_hi = VMERGE_VVM_F (t2_hi, small_hi, very_small, vlen);
  t2_lo = VMERGE_VVM_F (t2_lo, small_lo, very_small, vlen);

  VFLOAT z = VFNEG_V (VFADD_VV (t2_hi, t2_lo, vlen), vlen);
  VBOOL ge6 = __riscv_vmfge (x, RVV_ERFF_SATURATION_THRESHOLD_F, vlen);
  z = VFMERGE_VFM (z, fp_posOne_f, ge6, vlen);

  VUINT vclass = __riscv_vfclass (a, vlen);
  VBOOL is_inf;
  IDENTIFY (vclass, class_Inf, is_inf, vlen);
  z = VFMERGE_VFM (z, fp_posOne_f, is_inf, vlen);

  VBOOL eq0 = __riscv_vmfeq (a, 0.0f, vlen);
  z = VFMERGE_VFM (z, fp_posZero_f, eq0, vlen);
  z = __riscv_vfsgnjx (z, a, vlen);

  return z;
}

#undef LMUL
#undef SIMDLEN
