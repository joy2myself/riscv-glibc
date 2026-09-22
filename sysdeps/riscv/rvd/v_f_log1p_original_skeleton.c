/* Single-precision source-original SLEEF log1pf skeleton.

   This is a target-owned translation of SLEEF xlog1pf and its xlogf_u1
   large-input fallback.  The source reduction, polynomial, double-float
   operation order, and final fixup order are kept in this comparison entry.
   Clean-only quarantine, large-input repair, and tiny-input repair are not
   imported here.  */

#include <float.h>
#include <vecmath_config.h>
#include <v_math.h>

#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

/* The clean mixed-lane path uses the masked DF helpers from v_math.h.  */

static inline VINT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_log1pf_vilogb2k_m, LMUL)
  (VFLOAT d, size_t vlen)
{
  VUINT q = __riscv_vsrl (F_AS_U (d), MAN_LEN_F, vlen);
  q = __riscv_vand (q, 0xff, vlen);
  return __riscv_vsub (U_AS_I (q), EXP_BIAS_F, vlen);
}

static inline VFLOAT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_log1pf_vldexp3_m, LMUL)
  (VFLOAT d, VINT q, size_t vlen)
{
  /* SLEEF vldexp3 is an exponent-field add for a normalized float.  */
  return I_AS_F (__riscv_vadd (
      F_AS_I (d), __riscv_vsll (q, MAN_LEN_F, vlen), vlen));
}

static VFLOAT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_log1pf_xlogf_u1_m, LMUL)
  (VFLOAT d, size_t vlen)
{
  VBOOL subnormal = __riscv_vmflt (d, RVV_LOG1PF_FLT_MIN_F, vlen);
  VFLOAT d_scaled = VFMUL_VF (d, RVV_LOG1PF_SCALE_F, vlen);
  d = VMERGE_VVM_F (d, d_scaled, subnormal, vlen);

  VFLOAT d_adj = VFMUL_VF (d, RVV_LOG1PF_INV_075_F, vlen);
  VINT e = RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_log1pf_vilogb2k_m, LMUL)
    (d_adj, vlen);
  VFLOAT m = RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_log1pf_vldexp3_m, LMUL)
    (d, __riscv_vneg (e, vlen), vlen);
  e = VMERGE_VVM_I (e, __riscv_vsub (e, 64, vlen), subnormal, vlen);

  VFLOAT s_hi, s_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF (
    VFMV_VF (RVV_LOG1PF_LN2_HI_F, vlen),
    VFMV_VF (RVV_LOG1PF_LN2_LO_F, vlen),
    __riscv_vfcvt_f (e, vlen), s_hi, s_lo, vlen);

  VFLOAT num_hi, num_lo;
  VFLOAT den_hi, den_lo;
  RVV_SLEEF_ORIG_DFADD2_VF_VF (VFMV_VF (-1.0f, vlen), m,
                               num_hi, num_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF_VF (VFMV_VF (1.0f, vlen), m,
                               den_hi, den_lo, vlen);

  VFLOAT x_hi, x_lo;
  RVV_SLEEF_ORIG_DFDIV_VF2_VF2 (num_hi, num_lo, den_hi, den_lo,
                                x_hi, x_lo, vlen);
  VFLOAT x2 = VFMUL_VV (x_hi, x_hi, vlen);

  VFLOAT poly = VFMV_VF (RVV_SLEEF_LOG_U1_POLY0_F, vlen);
  poly = VFMADD_VV (poly, x2,
                    VFMV_VF (RVV_SLEEF_LOG_U1_POLY1_F, vlen), vlen);
  poly = VFMADD_VV (poly, x2,
                    VFMV_VF (RVV_SLEEF_LOG_U1_POLY2_F, vlen), vlen);

  VFLOAT poly_tail = VFMUL_VV (
    VFMUL_VV (x2, x_hi, vlen), poly, vlen);
  VFLOAT two_x_hi, two_x_lo;
  RVV_SLEEF_ORIG_DFSCALE_VF2_VF (x_hi, x_lo,
                                 VFMV_VF (2.0f, vlen),
                                 two_x_hi, two_x_lo, vlen);
  RVV_SLEEF_ORIG_DFADD_VF2_VF2 (s_hi, s_lo, two_x_hi, two_x_lo,
                                s_hi, s_lo, vlen);
  RVV_SLEEF_ORIG_DFADD_VF2_VF2_VF (s_hi, s_lo, poly_tail,
                                   s_hi, s_lo, vlen);
  VFLOAT r = VFADD_VV (s_hi, s_lo, vlen);

  /* SLEEF xlogf_u1 fixup order: +Inf, negative/NaN, then zero.  */
  /* SLEEF applies these fixups to its source-mutated d, after the
     subnormal scaling above.  Keep that observation point intact.  */
  VUINT bits = F_AS_U (d);
  VUINT abs_bits = __riscv_vand (bits, 0x7fffffff, vlen);
  VBOOL pos_inf = __riscv_vmseq (bits, 0x7f800000, vlen);
  VBOOL zero = __riscv_vmseq (abs_bits, 0, vlen);
  VBOOL negative = __riscv_vmflt (d, fp_posZero_f, vlen);
  VUINT vclass = __riscv_vfclass (d, vlen);
  VBOOL nan;
  IDENTIFY_F (vclass, class_NaN, nan, vlen);
  VBOOL invalid = __riscv_vmor (negative, nan, vlen);
  r = VFMERGE_VFM (r, fp_posInf_f, pos_inf, vlen);
  r = VFMERGE_VFM (r, fp_qNaN_f, invalid, vlen);
  r = VFMERGE_VFM (r, fp_negInf_f, zero, vlen);
  return r;
}

#if RVV_MATH_IS_LOG1PF_HELPERS_ONLY (V_MATH_FUNC)
/* Masked source-shaped xlogf_u1 fallback used by the clean mixed-lane path.
   The active mask is deliberately kept on every floating operation.  Integer
   exponent extraction and bit fixups are harmless on inactive lanes, but the
   inactive exponent lanes are normalized before the integer-to-float
   conversion so that no implementation-defined masked-off value can reach a
   floating operation.  */
static VFLOAT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_log1pf_xlogf_u1_masked_m, LMUL)
  (VFLOAT d, VBOOL active, size_t vlen)
{
  VBOOL subnormal = __riscv_vmflt (d, RVV_LOG1PF_FLT_MIN_F, vlen);
  subnormal = __riscv_vmand (subnormal, active, vlen);
  VFLOAT d_scaled = __riscv_vfmul (active, d,
                                    RVV_LOG1PF_SCALE_F, vlen);
  d = __riscv_vmerge (d, d_scaled, subnormal, vlen);

  VFLOAT d_adj = __riscv_vfmul
    (active, d, RVV_LOG1PF_INV_075_F, vlen);
  VINT e = RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_log1pf_vilogb2k_m, LMUL)
    (d_adj, vlen);
  VFLOAT m = RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_log1pf_vldexp3_m, LMUL)
    (d, __riscv_vneg (e, vlen), vlen);
  VINT e_sub = __riscv_vsub (e, 64, vlen);
  e = VMERGE_VVM_I (e, e_sub, subnormal, vlen);
  e = VMERGE_VVM_I (VMVI_VX (0, vlen), e, active, vlen);

  VFLOAT s_hi, s_lo;
  RVV_SLEEF_LOG1PF_DFMUL_VF2_VF_M
    (active,
     VFMV_VF (RVV_LOG1PF_LN2_HI_F, vlen),
     VFMV_VF (RVV_LOG1PF_LN2_LO_F, vlen),
     __riscv_vfcvt_f (e, vlen), s_hi, s_lo, vlen);

  VFLOAT num_hi, num_lo;
  VFLOAT den_hi, den_lo;
  RVV_SLEEF_LOG1PF_DFADD2_VF_VF_M
    (active, VFMV_VF (-1.0f, vlen), m, num_hi, num_lo, vlen);
  RVV_SLEEF_LOG1PF_DFADD2_VF_VF_M
    (active, VFMV_VF (1.0f, vlen), m, den_hi, den_lo, vlen);

  VFLOAT x_hi, x_lo;
  RVV_SLEEF_LOG1PF_DFDIV_VF2_VF2_M
    (active, num_hi, num_lo, den_hi, den_lo, x_hi, x_lo, vlen);
  VFLOAT x2 = RVV_SLEEF_LOG1PF_MUL_VV_M (active, x_hi, x_hi, vlen);

  VFLOAT poly = VFMV_VF (RVV_SLEEF_LOG_U1_POLY0_F, vlen);
  poly = RVV_SLEEF_LOG1PF_MADD_VV_M
    (active, poly, x2, VFMV_VF (RVV_SLEEF_LOG_U1_POLY1_F, vlen), vlen);
  poly = RVV_SLEEF_LOG1PF_MADD_VV_M
    (active, poly, x2, VFMV_VF (RVV_SLEEF_LOG_U1_POLY2_F, vlen), vlen);

  VFLOAT poly_tail = RVV_SLEEF_LOG1PF_MUL_VV_M
    (active,
     RVV_SLEEF_LOG1PF_MUL_VV_M (active, x2, x_hi, vlen), poly, vlen);
  VFLOAT two_x_hi, two_x_lo;
  RVV_SLEEF_LOG1PF_DFSCALE_VF2_VF_M
    (active, x_hi, x_lo, VFMV_VF (2.0f, vlen), two_x_hi, two_x_lo, vlen);
  RVV_SLEEF_LOG1PF_DFADD_VF2_VF2_M
    (active, s_hi, s_lo, two_x_hi, two_x_lo, s_hi, s_lo, vlen);
  RVV_SLEEF_LOG1PF_DFADD_VF2_VF2_VF_M
    (active, s_hi, s_lo, poly_tail, s_hi, s_lo, vlen);
  VFLOAT r = RVV_SLEEF_LOG1PF_ADD_VV_M (active, s_hi, s_lo, vlen);

  VUINT bits = F_AS_U (d);
  VUINT abs_bits = __riscv_vand (bits, 0x7fffffff, vlen);
  VBOOL pos_inf = __riscv_vmseq (bits, 0x7f800000, vlen);
  pos_inf = __riscv_vmand (pos_inf, active, vlen);
  VBOOL zero = __riscv_vmseq (abs_bits, 0, vlen);
  zero = __riscv_vmand (zero, active, vlen);
  VBOOL negative = __riscv_vmflt (d, fp_posZero_f, vlen);
  negative = __riscv_vmand (negative, active, vlen);
  VUINT vclass = __riscv_vfclass (d, vlen);
  VBOOL nan;
  IDENTIFY_F (vclass, class_NaN, nan, vlen);
  nan = __riscv_vmand (nan, active, vlen);
  VBOOL invalid = __riscv_vmor (negative, nan, vlen);
  r = VFMERGE_VFM (r, fp_posInf_f, pos_inf, vlen);
  r = VFMERGE_VFM (r, fp_qNaN_f, invalid, vlen);
  r = VFMERGE_VFM (r, fp_negInf_f, zero, vlen);
  return r;
}

#endif /* RVV_MATH_IS_LOG1PF_HELPERS_ONLY */

static VFLOAT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_log1pf_kernel_m, LMUL)
  (VFLOAT d, size_t vlen)
{
  /* SLEEF xlog1pf starts with the rounded dp1 operation.  */
  VFLOAT dp1 = VFADD_VF (d, 1.0f, vlen);
  VBOOL subnormal = __riscv_vmflt (dp1, RVV_LOG1PF_FLT_MIN_F, vlen);
  dp1 = VMERGE_VVM_F (
    dp1, VFMUL_VF (dp1, RVV_LOG1PF_SCALE_F, vlen),
    subnormal, vlen);

  VFLOAT dp1_adj = VFMUL_VF (dp1, RVV_LOG1PF_INV_075_F, vlen);
  VINT e = RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_log1pf_vilogb2k_m, LMUL)
    (dp1_adj, vlen);
  VFLOAT t = RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_log1pf_vldexp3_m, LMUL)
    (VFMV_VF (1.0f, vlen), __riscv_vneg (e, vlen), vlen);
  VFLOAT m = VFMADD_VV (d, t, VFSUB_VF (t, 1.0f, vlen), vlen);
  e = VMERGE_VVM_I (e, __riscv_vsub (e, 64, vlen), subnormal, vlen);

  VFLOAT s_hi, s_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF (
    VFMV_VF (RVV_LOG1PF_LN2_HI_F, vlen),
    VFMV_VF (RVV_LOG1PF_LN2_LO_F, vlen),
    __riscv_vfcvt_f (e, vlen), s_hi, s_lo, vlen);

  VFLOAT den_hi, den_lo;
  RVV_SLEEF_ORIG_DFADD_VF_VF (VFMV_VF (2.0f, vlen), m,
                              den_hi, den_lo, vlen);
  VFLOAT x_hi, x_lo;
  RVV_SLEEF_ORIG_DFDIV_VF2_VF2 (
    m, VFMV_VF (0.0f, vlen), den_hi, den_lo, x_hi, x_lo, vlen);
  VFLOAT x2 = VFMUL_VV (x_hi, x_hi, vlen);

  VFLOAT poly = VFMV_VF (RVV_SLEEF_LOG_U1_POLY0_F, vlen);
  poly = VFMADD_VV (poly, x2,
                    VFMV_VF (RVV_SLEEF_LOG_U1_POLY1_F, vlen), vlen);
  poly = VFMADD_VV (poly, x2,
                    VFMV_VF (RVV_SLEEF_LOG_U1_POLY2_F, vlen), vlen);
  VFLOAT poly_tail = VFMUL_VV (
    VFMUL_VV (x2, x_hi, vlen), poly, vlen);

  VFLOAT two_x_hi, two_x_lo;
  RVV_SLEEF_ORIG_DFSCALE_VF2_VF (x_hi, x_lo,
                                 VFMV_VF (2.0f, vlen),
                                 two_x_hi, two_x_lo, vlen);
  RVV_SLEEF_ORIG_DFADD_VF2_VF2 (s_hi, s_lo, two_x_hi, two_x_lo,
                                s_hi, s_lo, vlen);
  RVV_SLEEF_ORIG_DFADD_VF2_VF2_VF (s_hi, s_lo, poly_tail,
                                   s_hi, s_lo, vlen);
  VFLOAT r = VFADD_VV (s_hi, s_lo, vlen);

  /* Keep SLEEF's vector xlogf_u1 fallback as a target-owned helper.  */
  VBOOL core = __riscv_vmfle (d, RVV_LOG1PF_BOUND_F, vlen);
  /* Match SLEEF's vtestallones guard: do not evaluate xlogf_u1 when every
     lane is in the main xlog1pf core.  This also preserves the source
     operation/fenv boundary for the all-core case.  */
  if (__riscv_vcpop (core, vlen) != vlen)
    {
      VFLOAT fallback = RVV_SLEEF_ORIG_MANGLE
        (rvv_sleef_orig_log1pf_xlogf_u1_m, LMUL) (d, vlen);
      r = VMERGE_VVM_F (r, fallback, __riscv_vmnot (core, vlen), vlen);
    }

  VBOOL less_than_neg_one = __riscv_vmflt (d, fp_negOne_f, vlen);
  /* Source xlog1pf uses vor_vm_vo32_vm here, which sets every result bit in
     the domain-error lanes instead of selecting the framework qNaN value.  */
  VUINT all_ones = VMVU_VX (~(uint32_t) 0, vlen);
  VUINT domain_bits = __riscv_vmerge
    (VMVU_VX (0, vlen), all_ones, less_than_neg_one, vlen);
  r = U_AS_F (__riscv_vor (F_AS_U (r), domain_bits, vlen));
  VBOOL equal_neg_one = __riscv_vmfeq (d, fp_negOne_f, vlen);
  r = VFMERGE_VFM (r, fp_negInf_f, equal_neg_one, vlen);
  VUINT bits = F_AS_U (d);
  VBOOL neg_zero = __riscv_vmseq (bits, 0x80000000u, vlen);
  r = VFMERGE_VFM (r, fp_negZero_f, neg_zero, vlen);

  return r;
}

#if RVV_MATH_IS_LOG1PF_HELPERS_ONLY (V_MATH_FUNC)
/* Source-shaped finite kernel with an explicit active mask.  This is used
   only for clean vectors containing a mixture of finite and special lanes;
   it is not exposed as the _original entry point.  */
static VFLOAT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_log1pf_kernel_masked_m, LMUL)
  (VFLOAT d, VBOOL active, size_t vlen)
{
  VFLOAT dp1 = RVV_SLEEF_LOG1PF_ADD_VF_M (active, d, 1.0f, vlen);
  dp1 = VFMERGE_VFM (dp1, 1.0f, __riscv_vmnot (active, vlen), vlen);

  VBOOL subnormal = __riscv_vmflt (dp1, RVV_LOG1PF_FLT_MIN_F, vlen);
  subnormal = __riscv_vmand (subnormal, active, vlen);
  VFLOAT dp1_scaled = RVV_SLEEF_LOG1PF_MUL_VF_M
    (active, dp1, RVV_LOG1PF_SCALE_F, vlen);
  dp1 = __riscv_vmerge (dp1, dp1_scaled, subnormal, vlen);

  VFLOAT dp1_adj = RVV_SLEEF_LOG1PF_MUL_VF_M
    (active, dp1, RVV_LOG1PF_INV_075_F, vlen);
  VINT e = RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_log1pf_vilogb2k_m, LMUL)
    (dp1_adj, vlen);
  VFLOAT t = RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_log1pf_vldexp3_m, LMUL)
    (VFMV_VF (1.0f, vlen), __riscv_vneg (e, vlen), vlen);
  VFLOAT t_minus_one = RVV_SLEEF_LOG1PF_SUB_VF_M
    (active, t, 1.0f, vlen);
  VFLOAT m = RVV_SLEEF_LOG1PF_MADD_VV_M
    (active, d, t, t_minus_one, vlen);
  e = VMERGE_VVM_I (e, __riscv_vsub (e, 64, vlen), subnormal, vlen);
  e = VMERGE_VVM_I (VMVI_VX (0, vlen), e, active, vlen);

  VFLOAT s_hi, s_lo;
  RVV_SLEEF_LOG1PF_DFMUL_VF2_VF_M
    (active,
     VFMV_VF (RVV_LOG1PF_LN2_HI_F, vlen),
     VFMV_VF (RVV_LOG1PF_LN2_LO_F, vlen),
     __riscv_vfcvt_f (e, vlen), s_hi, s_lo, vlen);

  VFLOAT den_hi, den_lo;
  RVV_SLEEF_LOG1PF_DFADD_VF_VF_M
    (active, VFMV_VF (2.0f, vlen), m, den_hi, den_lo, vlen);
  VFLOAT x_hi, x_lo;
  RVV_SLEEF_LOG1PF_DFDIV_VF2_VF2_M
    (active, m, VFMV_VF (0.0f, vlen), den_hi, den_lo, x_hi, x_lo, vlen);
  VFLOAT x2 = RVV_SLEEF_LOG1PF_MUL_VV_M (active, x_hi, x_hi, vlen);

  VFLOAT poly = VFMV_VF (RVV_SLEEF_LOG_U1_POLY0_F, vlen);
  poly = RVV_SLEEF_LOG1PF_MADD_VV_M
    (active, poly, x2, VFMV_VF (RVV_SLEEF_LOG_U1_POLY1_F, vlen), vlen);
  poly = RVV_SLEEF_LOG1PF_MADD_VV_M
    (active, poly, x2, VFMV_VF (RVV_SLEEF_LOG_U1_POLY2_F, vlen), vlen);
  VFLOAT poly_tail = RVV_SLEEF_LOG1PF_MUL_VV_M
    (active,
     RVV_SLEEF_LOG1PF_MUL_VV_M (active, x2, x_hi, vlen), poly, vlen);

  VFLOAT two_x_hi, two_x_lo;
  RVV_SLEEF_LOG1PF_DFSCALE_VF2_VF_M
    (active, x_hi, x_lo, VFMV_VF (2.0f, vlen), two_x_hi, two_x_lo, vlen);
  RVV_SLEEF_LOG1PF_DFADD_VF2_VF2_M
    (active, s_hi, s_lo, two_x_hi, two_x_lo, s_hi, s_lo, vlen);
  RVV_SLEEF_LOG1PF_DFADD_VF2_VF2_VF_M
    (active, s_hi, s_lo, poly_tail, s_hi, s_lo, vlen);
  VFLOAT r = RVV_SLEEF_LOG1PF_ADD_VV_M (active, s_hi, s_lo, vlen);

  VBOOL core = __riscv_vmfle (d, RVV_LOG1PF_BOUND_F, vlen);
  core = __riscv_vmand (core, active, vlen);
  VBOOL fallback_mask
    = __riscv_vmand (active, __riscv_vmnot (core, vlen), vlen);
  if (__riscv_vcpop (fallback_mask, vlen) != 0)
    {
      VFLOAT fallback = RVV_SLEEF_ORIG_MANGLE
        (rvv_sleef_orig_log1pf_xlogf_u1_masked_m, LMUL)
        (d, fallback_mask, vlen);
      r = VMERGE_VVM_F (r, fallback, fallback_mask, vlen);
    }

  VBOOL less_than_neg_one = __riscv_vmflt (d, fp_negOne_f, vlen);
  less_than_neg_one = __riscv_vmand (less_than_neg_one, active, vlen);
  VUINT all_ones = VMVU_VX (~(uint32_t) 0, vlen);
  VUINT domain_bits = __riscv_vmerge
    (VMVU_VX (0, vlen), all_ones, less_than_neg_one, vlen);
  r = U_AS_F (__riscv_vor (F_AS_U (r), domain_bits, vlen));
  VBOOL equal_neg_one = __riscv_vmfeq (d, fp_negOne_f, vlen);
  equal_neg_one = __riscv_vmand (equal_neg_one, active, vlen);
  r = VFMERGE_VFM (r, fp_negInf_f, equal_neg_one, vlen);
  VUINT bits = F_AS_U (d);
  VBOOL neg_zero = __riscv_vmseq (bits, 0x80000000u, vlen);
  neg_zero = __riscv_vmand (neg_zero, active, vlen);
  r = VFMERGE_VFM (r, fp_negZero_f, neg_zero, vlen);

  return r;
}

#endif /* RVV_MATH_IS_LOG1PF_HELPERS_ONLY */

#if !RVV_MATH_IS_LOG1PF_HELPERS_ONLY (V_MATH_FUNC)
static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT d, size_t vlen)
{
  return RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_log1pf_kernel_m, LMUL) (d, vlen);
}
#endif

#if !RVV_MATH_IS_LOG1PF_HELPERS_ONLY (V_MATH_FUNC)
#undef LMUL
#undef SIMDLEN
#endif
