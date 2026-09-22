/* Single-precision source-original SLEEF lgammaf skeleton.

   This is a target-owned translation of SLEEF gammafk/xlgammaf_u1.  The
   gamma reduction, source double-float helpers, reflection denominator and
   late special-result selection are shared with the independently translated
   tgammaf source kernel, but the public lgammaf comparison entry is separate
   and has no dependency on a SLEEF runtime symbol or source include.

   The helper-only import is deliberate: it reuses only target-owned gammafk,
   logk2f and their source constants, while omitting the tgammaf public body.
   Clean-side EH quarantine and finite-core repairs are not imported into this
   source-original comparator.

   Copyright (C) 2026 Free Software Foundation, Inc.  */

#include <vecmath_config.h>
#include <v_math.h>

#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT a, size_t vlen)
{
  VFLOAT clc_hi, clc_lo, clln_hi, clln_lo, clld_hi, clld_lo;
  RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_tgammaf_gammafk_m, LMUL)
    (a, &clc_hi, &clc_lo, &clln_hi, &clln_lo, &clld_hi, &clld_lo,
     vlen);

  /* gammafk returns clc, clln and clld.  xlgammaf_u1 first forms the
     double-float quotient clln/clld (the value returned as gammafk's second
     pair), then applies dfabs to that quotient before logk2f.  Applying
     logk2f directly to clld loses the numerator and is not the SLEEF
     algorithm.  */
  VFLOAT quotient_hi, quotient_lo;
  RVV_SLEEF_ORIG_DFDIV_VF2_VF2
    (clln_hi, clln_lo, clld_hi, clld_lo, quotient_hi, quotient_lo, vlen);

  /* SLEEF dfabs flips only the low component's sign when the high component
     is negative; the low component is a residual and must not be independently
     fabs'ed.  */
  VFLOAT abs_quotient_hi = __riscv_vfabs (quotient_hi, vlen);
  VUINT quotient_hi_sign = __riscv_vand
    (F_AS_U (quotient_hi), 0x80000000u, vlen);
  VFLOAT abs_quotient_lo = U_AS_F
    (__riscv_vxor (F_AS_U (quotient_lo), quotient_hi_sign, vlen));
  VFLOAT log_den_hi, log_den_lo;
  RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_tgammaf_logk2f_m, LMUL)
    (abs_quotient_hi, abs_quotient_lo, &log_den_hi, &log_den_lo, vlen);

  VFLOAT y_hi, y_lo;
  RVV_SLEEF_ORIG_DFADD2_VF2_VF2
    (clc_hi, clc_lo, log_den_hi, log_den_lo, y_hi, y_lo, vlen);
  VFLOAT result = VFADD_VV (y_hi, y_lo, vlen);

  /* SLEEF xlgammaf_u1's source final selection is:
       +/-Inf input, non-positive integer input, or a numeric finite input
       whose gamma reduction produced NaN -> +Inf.
     The RTZ integer test intentionally remains source-shaped and is not
     moved behind a clean-side conversion quarantine.  */
  VUINT aclass = __riscv_vfclass (a, vlen);
  VBOOL a_inf, a_nan;
  IDENTIFY_F (aclass, class_Inf, a_inf, vlen);
  IDENTIFY_F (aclass, class_NaN, a_nan, vlen);
  VBOOL a_number = __riscv_vmnot
    (__riscv_vmor (a_inf, a_nan, vlen), vlen);

  VINT a_trunc_i = RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_tgammaf_truncate_m, LMUL) (a, vlen);
  VFLOAT a_trunc = __riscv_vfcvt_f (a_trunc_i, vlen);
  VBOOL a_is_int = __riscv_vmfeq (a_trunc, a, vlen);
  VBOOL a_nonpositive = __riscv_vmfle (a, 0.0f, vlen);
  VBOOL negative_integer = __riscv_vmand (a_nonpositive, a_is_int, vlen);
  VBOOL reduction_nan = __riscv_vmfne (result, result, vlen);
  VBOOL finite_nan = __riscv_vmand
    (a_number, reduction_nan, vlen);

  VBOOL inf_result = __riscv_vmor
    (a_inf, __riscv_vmor (negative_integer, finite_nan, vlen),
     vlen);
  result = VFMERGE_VFM (result, fp_posInf_f, inf_result, vlen);

  return result;
}

#undef LMUL
#undef SIMDLEN
