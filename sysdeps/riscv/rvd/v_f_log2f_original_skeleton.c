/* Single-precision source-original SLEEF log2f skeleton.

   This is a target-owned translation of SLEEF xlog2f.  The source finite
   reduction, constants, double-float operations, and late fixup order are
   kept separate from the clean implementation's EH and range adaptations.
   SLEEF headers and source files are not included or linked.  */

#include <vecmath_config.h>
#include <v_math.h>

#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

/* Reusable target-owned translation of SLEEF xlog2f.  */
static VFLOAT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_log2f_kernel_m, LMUL)
  (VFLOAT d, size_t vlen)
{
  /* Keep SLEEF's mutated d operand: its late fixups inspect this value.  */
  VBOOL scaled_input = __riscv_vmflt
    (d, RVV_LOG2F_FLT_MIN_F, vlen);
#if RVV_MATH_IS_LOG2P1F_CLEAN (V_MATH_FUNC)
  /* Clean callers quarantine the scale operation itself so a large finite
     lane cannot overflow in the inactive source subnormal path.  The
     source-original translation leaves this operation unmasked.  */
  VFLOAT d_scaled = VFMUL_VF_M
    (scaled_input, d, RVV_LOG2F_SCALE_F, vlen);
#else
  VFLOAT d_scaled = VFMUL_VF (d, RVV_LOG2F_SCALE_F, vlen);
#endif
  d = VMERGE_VVM_F (d, d_scaled, scaled_input, vlen);

  VFLOAT d_adj = VFMUL_VF (d, RVV_LOG2F_INV_075_F, vlen);
  VINT e = RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_log2f_vilogb2k_m, LMUL) (d_adj, vlen);
  VFLOAT m = RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_log2f_vldexp3_m, LMUL)
    (d, __riscv_vneg (e, vlen), vlen);
  e = VMERGE_VVM_I
    (e, __riscv_vsub (e, RVV_LOG2F_SUBNORMAL_EXP_ADJUST, vlen), scaled_input,
     vlen);

  /* x = (m - 1) / (m + 1), preserving SLEEF's scalar-first dfadd2 order.  */
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

  /* SLEEF source order: t=C0; t=vmla(t,x2,C1); t=vmla(t,x2,C2).
     With vmla(x,y,z) = z + x*y this is C2 + x2*(C1 + x2*C0).  */
  VFLOAT t = VFMV_VF (RVV_LOG2F_POLY0_F, vlen);
  t = VFMADD_VV (t, x2,
                 VFMV_VF (RVV_LOG2F_POLY1_F, vlen), vlen);
  t = VFMADD_VV (t, x2,
                 VFMV_VF (RVV_LOG2F_POLY2_F, vlen), vlen);

  VFLOAT ef = __riscv_vfcvt_f (e, vlen);
  VFLOAT scale_hi, scale_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
     (x_hi, x_lo,
     VFMV_VF (RVV_LOG2F_X_HI_F, vlen),
     VFMV_VF (RVV_LOG2F_X_LO_F, vlen),
     scale_hi, scale_lo, vlen);

  /* SLEEF's non-AVX512 branch does dfadd2(float(e), pair).  */
  VFLOAT s_hi, s_lo;
  RVV_SLEEF_ORIG_DFADD2_VF_VF2 (ef, scale_hi, scale_lo,
                                s_hi, s_lo, vlen);

  VFLOAT poly_tail = VFMUL_VV
    (VFMUL_VV (x2, x_hi, vlen), t, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF (s_hi, s_lo, poly_tail,
                                s_hi, s_lo, vlen);
  VFLOAT result = VFADD_VV (s_hi, s_lo, vlen);

  /* SLEEF xlog2f fixup order: +Inf, negative/NaN, then zero.  */
  VUINT vclass = __riscv_vfclass (d, vlen);
  VBOOL is_inf, is_nan, is_zero;
  IDENTIFY_F (vclass, class_posInf, is_inf, vlen);
  IDENTIFY_F (vclass, class_NaN, is_nan, vlen);
  VBOOL is_negative = __riscv_vmflt (d, fp_posZero_f, vlen);
  VBOOL invalid = __riscv_vmor (is_negative, is_nan, vlen);
  is_zero = __riscv_vmfeq (d, fp_posZero_f, vlen);

  result = VFMERGE_VFM (result, fp_posInf_f, is_inf, vlen);
  result = VFMERGE_VFM (result, fp_qNaN_f, invalid, vlen);
  result = VFMERGE_VFM (result, fp_negInf_f, is_zero, vlen);

  return result;
}

#if RVV_MATH_IS_LOG2F_ORIGINAL (V_MATH_FUNC)
static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT d, size_t vlen)
{
  return RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_log2f_kernel_m, LMUL) (d, vlen);
}
#endif

#if RVV_MATH_IS_LOG2F_ORIGINAL (V_MATH_FUNC)
#undef LMUL
#undef SIMDLEN
#endif
