/* Single-precision source-original SLEEF hypotf skeleton.

   This is a target-owned translation of SLEEF xhypotf_u05.  The source
   finite-kernel dataflow and source fixup order are kept separate from the
   clean implementation's EH quarantine and precision adaptations.  */

#include <vecmath_config.h>
#include <v_math.h>

#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT vx, VFLOAT vy,
                                           size_t vlen)
{
  /* SLEEF xhypotf_u05 starts with absolute values and retains the
     pre-scaling min/max values for its later fixups.  */
  VFLOAT x = __riscv_vfabs (vx, vlen);
  VFLOAT y = __riscv_vfabs (vy, vlen);
  VFLOAT vmin = __riscv_vfmin (x, y, vlen);
  VFLOAT n = vmin;
  VFLOAT vmax = __riscv_vfmax (x, y, vlen);
  VFLOAT d = vmax;

  /* Source xhypotf_u05 scales both ratio operands for subnormal max.  */
  VBOOL scale_mask
      = __riscv_vmflt (vmax, RVV_HYPOTF_SOURCE_FLT_MIN_F, vlen);
  VFLOAT n_scaled
      = VFMUL_VF (n, RVV_HYPOTF_SOURCE_SCALE_F, vlen);
  VFLOAT d_scaled
      = VFMUL_VF (d, RVV_HYPOTF_SOURCE_SCALE_F, vlen);
  n = VMERGE_VVM_F (n, n_scaled, scale_mask, vlen);
  d = VMERGE_VVM_F (d, d_scaled, scale_mask, vlen);

  /* t = n/d in the source double-float representation.  */
  VFLOAT zero = VFMV_VF (fp_posZero_f, vlen);
  VFLOAT t_hi, t_lo;
  RVV_SLEEF_ORIG_DFDIV_VF2_VF2
    (n, zero, d, zero, t_hi, t_lo, vlen);

  /* t = max * sqrt(t*t + 1), preserving the source operation order.  */
  VFLOAT t2_hi, t2_lo;
  RVV_SLEEF_ORIG_DFSQU_VF2 (t_hi, t_lo, t2_hi, t2_lo, vlen);

  VFLOAT sum_hi, sum_lo;
  RVV_SLEEF_ORIG_DFADD2_VF2_VF_S
    (t2_hi, t2_lo, fp_posOne_f, sum_hi, sum_lo, vlen);

  VFLOAT sqrt_hi, sqrt_lo;
  RVV_SLEEF_ORIG_DFSQRT_VF2_VF2
    (sum_hi, sum_lo, sqrt_hi, sqrt_lo, vlen);

  VFLOAT result_hi, result_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (sqrt_hi, sqrt_lo, vmax, result_hi, result_lo, vlen);
  VFLOAT result = VFADD_VV (result_hi, result_lo, vlen);

  /* Keep SLEEF's post-core selection order:
       computed NaN -> +Inf, min == 0 -> max,
       any input NaN -> NaN, any input Inf -> +Inf.  */
  VBOOL result_nan = __riscv_vmfne (result, result, vlen);
  result = VFMERGE_VFM (result, fp_posInf_f, result_nan, vlen);

  VBOOL min_zero = __riscv_vmfeq (vmin, fp_posZero_f, vlen);
  result = VMERGE_VVM_F (result, vmax, min_zero, vlen);

  VBOOL x_nan = __riscv_vmfne (x, x, vlen);
  VBOOL y_nan = __riscv_vmfne (y, y, vlen);
  VBOOL any_nan = __riscv_vmor (x_nan, y_nan, vlen);
  result = VFMERGE_VFM (result, fp_qNaN_f, any_nan, vlen);

  VBOOL x_inf = __riscv_vmfeq (x, fp_posInf_f, vlen);
  VBOOL y_inf = __riscv_vmfeq (y, fp_posInf_f, vlen);
  VBOOL any_inf = __riscv_vmor (x_inf, y_inf, vlen);
  result = VFMERGE_VFM (result, fp_posInf_f, any_inf, vlen);

  return result;
}

#undef LMUL
#undef SIMDLEN
