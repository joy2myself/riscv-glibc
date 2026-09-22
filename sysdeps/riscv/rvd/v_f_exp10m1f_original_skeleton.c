/* Single-precision SLEEF-component composite for exp10m1f.

   SLEEF has no native whole-routine exp10m1f entry.  The comparison contract
   is a source-shaped composition: form x*ln(10) with the source split and
   use the target-owned translations of SLEEF xexpm1f and xexp10f.  The
   small-input branch retains the transformed argument as a double-float pair
   for xexpm1f; the complementary branch computes xexp10f and subtracts one.
   This remains explicitly distinct from a native whole-routine SLEEF entry
   and from the clean implementation's EH quarantine and endpoint repairs.

   Copyright (C) 2026 Free Software Foundation, Inc.  */

#include <vecmath_config.h>
#include <v_math.h>

#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

#define RVV_SLEEF_ORIG_EXPM1F_HELPERS_ONLY
#include "v_f_expm1f_original_skeleton.c"
#undef RVV_SLEEF_ORIG_EXPM1F_HELPERS_ONLY

#define RVV_SLEEF_ORIG_EXP10F_HELPERS_ONLY
#include "v_f_exp10_original_skeleton.c"
#undef RVV_SLEEF_ORIG_EXP10F_HELPERS_ONLY

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  /* Source xexp10f carries this split ln(10) into its double-float
     construction.  Keep the transformed-argument arithmetic in source
     order; the reused expk2f helper supplies its own explicit RNE quotient
     conversion while ordinary arithmetic observes caller FRM.  */
  VFLOAT abs_x = __riscv_vfabs (x, vlen);
  VBOOL small = __riscv_vmflt (abs_x, 0.5f, vlen);

  /* Keep x*ln(10) as a double-float pair until expk2f consumes it.  The
     earlier composite collapsed this product to one binary32 value before
     entering xexpm1f; that is source-shaped at the API boundary but loses
     enough argument precision for the error to be amplified at large
     positive x.  The expm1 component is evaluated only on the small-
     magnitude side of the composite; the other side uses the target-owned
     xexp10f component before subtracting one, avoiding cancellation and
     transformed-argument overflow.  */
  VFLOAT x_expm1 = VFMERGE_VFM (x, 0.0f, __riscv_vmnot (small, vlen), vlen);
  VFLOAT a_hi, a_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (x_expm1, VFMV_VF (0.0f, vlen),
     VFMV_VF
       (RVV_SLEEF_ORIG_F64_TO_F32 (RVV_EXP10F_SOURCE_LN10_HI), vlen),
     VFMV_VF
       (RVV_SLEEF_ORIG_F64_TO_F32 (RVV_EXP10F_SOURCE_LN10_LO), vlen),
     a_hi, a_lo, vlen);
  VFLOAT a = VFADD_VV (a_hi, a_lo, vlen);

  VFLOAT exp_hi, exp_lo;
  RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_expm1f_expk2f_m, LMUL)
    (a_hi, a_lo, &exp_hi, &exp_lo, vlen);

  VFLOAT d_hi, d_lo;
  RVV_SLEEF_ORIG_DFADD2_VF2_VF_S
    (exp_hi, exp_lo, -1.0f, d_hi, d_lo, vlen);
  VFLOAT expm1_result = VFADD_VV (d_hi, d_lo, vlen);

  VFLOAT exp10_value = RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_exp10f_kernel_m, LMUL) (x, vlen);
  VFLOAT exp10_result = VFSUB_VF (exp10_value, 1.0f, vlen);
  /* Helper calls use v0 for their own predicates.  Re-form the selection
     predicate after those calls instead of carrying a live mask across the
     source-component graph.  This preserves the source selection and avoids
     an RVV mask live-range hazard under non-RNE caller FRM.  */
  VBOOL select_expm1 = __riscv_vmflt
    (__riscv_vfabs (x, vlen), 0.5f, vlen);
  VFLOAT result = VMERGE_VVM_F
    (exp10_result, expm1_result, select_expm1, vlen);

  /* xexpm1f's source late-selection order, now in transformed-argument
     space.  */
  VBOOL overflow = __riscv_vmfgt
    (a, RVV_SLEEF_EXPM1F_MAX, vlen);
  result = VFMERGE_VFM (result, fp_posInf_f, overflow, vlen);

  VBOOL underflow = __riscv_vmflt
    (a, RVV_SLEEF_EXPM1F_MIN, vlen);
  result = VFMERGE_VFM (result, fp_negOne_f, underflow, vlen);

  /* The transformed argument is meaningful only for the small-input
     expm1 branch.  In the complementary branch x_expm1 is deliberately
     sanitized to +0; under a downward caller FRM the split product can
     nevertheless carry a -0 low word, which must not select the expm1
     endpoint for the exp10 branch.  */
  VBOOL negzero = __riscv_vmseq
    (F_AS_U (a), VMVU_VX (0x80000000u, vlen), vlen);
  negzero = __riscv_vmand (negzero, small, vlen);
  result = VFMERGE_VFM (result, fp_negZero_f, negzero, vlen);

  /* The source-shaped split x*ln(10) is a composite boundary rather than
     SLEEF's native xexpm1f argument.  For raw infinities, the two split
     products can become Inf + -Inf and manufacture NaN.  Preserve the
     mathematical endpoint of the composite using the original input.  The
     same raw-input check is needed for -0 because the low split product may
     erase its sign when added to the high product.  Finite lanes retain the
     source component order above.  */
  VBOOL x_inf = __riscv_vmfeq
    (__riscv_vfabs (x, vlen), fp_posInf_f, vlen);
  VBOOL x_pos_inf = __riscv_vmfeq (x, fp_posInf_f, vlen);
  VFLOAT inf_endpoint = VFMERGE_VFM
    (VFMV_VF (fp_negOne_f, vlen), fp_posInf_f, x_pos_inf, vlen);
  result = VMERGE_VVM_F (result, inf_endpoint, x_inf, vlen);

  VBOOL raw_negzero = __riscv_vmseq
    (F_AS_U (x), VMVU_VX (0x80000000u, vlen), vlen);
  result = VFMERGE_VFM (result, fp_negZero_f, raw_negzero, vlen);

  return result;
}

#undef LMUL
#undef SIMDLEN
