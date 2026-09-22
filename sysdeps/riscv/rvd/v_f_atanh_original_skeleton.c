/* Single-precision source-original SLEEF atanhf skeleton.

   The finite path follows SLEEF xatanhf/logk2f.  Special values are handled
   in the source entry's post-core bitwise style; clean's early quarantine and
   repair policy is intentionally not included.  */

#include <vecmath_config.h>
#include <v_math.h>

#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

static inline VINT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_atanh_vilogbk_m, LMUL)
  (VFLOAT d, size_t vlen)
{
  VBOOL small = __riscv_vmflt (d, 5.421010862427522e-20f, vlen);
  VFLOAT scaled = VFMUL_VF (d, 1.8446744073709552e19f, vlen);
  d = __riscv_vmerge (d, scaled, small, vlen);
  VUINT q = __riscv_vsrl (F_AS_U (d), MAN_LEN_F, vlen);
  q = __riscv_vand (q, 0xff, vlen);
  VINT bias = __riscv_vmerge (VMVI_VX (0x7f, vlen),
                              VMVI_VX (64 + 0x7f, vlen), small, vlen);
  return __riscv_vsub (U_AS_I (q), bias, vlen);
}

static inline VFLOAT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_atanh_vpow2i_m, LMUL)
  (VINT q, size_t vlen)
{
  VINT biased = __riscv_vadd (q, EXP_BIAS_F, vlen);
  return I_AS_F (__riscv_vsll (biased, MAN_LEN_F, vlen));
}

static inline void
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_atanh_logk2f_m, LMUL)
  (VFLOAT d_hi, VFLOAT d_lo, VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT scaled_hi = VFMUL_VF (d_hi, 1.0f / 0.75f, vlen);
  VINT e = RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_atanh_vilogbk_m, LMUL)
    (scaled_hi, vlen);
  VFLOAT scale = RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_atanh_vpow2i_m, LMUL)
    (__riscv_vneg (e, vlen), vlen);

  VFLOAT m_hi, m_lo;
  RVV_SLEEF_ORIG_DFSCALE_VF2_VF (d_hi, d_lo, scale, m_hi, m_lo, vlen);

  VFLOAT n_hi, n_lo, q_hi, q_lo;
  RVV_SLEEF_ORIG_DFADD2_VF2_VF (m_hi, m_lo,
                                VFMV_VF (-1.0f, vlen), n_hi, n_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF (m_hi, m_lo,
                                VFMV_VF (1.0f, vlen), q_hi, q_lo, vlen);

  VFLOAT x_hi, x_lo;
  RVV_SLEEF_ORIG_DFDIV_VF2_VF2 (n_hi, n_lo, q_hi, q_lo,
                                x_hi, x_lo, vlen);
  VFLOAT x2_hi, x2_lo;
  RVV_SLEEF_ORIG_DFSQU_VF2 (x_hi, x_lo, x2_hi, x2_lo, vlen);

  VFLOAT t = VFMV_VF (RVV_ATANHF_LOG_POLY0_F, vlen);
  t = VFMADD_VV (t, x2_hi,
                 VFMV_VF (RVV_ATANHF_LOG_POLY1_F, vlen), vlen);
  t = VFMADD_VV (t, x2_hi,
                 VFMV_VF (RVV_ATANHF_LOG_POLY2_F, vlen), vlen);
  t = VFMADD_VV (t, x2_hi,
                 VFMV_VF (RVV_ATANHF_LOG_POLY3_F, vlen), vlen);

  VFLOAT s_hi, s_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF (
    VFMV_VF (RVV_LOG1PF_LN2_HI_F, vlen),
    VFMV_VF (RVV_LOG1PF_LN2_LO_F, vlen),
    __riscv_vfcvt_f (e, vlen), s_hi, s_lo, vlen);

  VFLOAT two_x_hi, two_x_lo;
  RVV_SLEEF_ORIG_DFSCALE_VF2_VF (x_hi, x_lo, VFMV_VF (2.0f, vlen),
                                 two_x_hi, two_x_lo, vlen);
  RVV_SLEEF_ORIG_DFADD_VF2_VF2 (s_hi, s_lo, two_x_hi, two_x_lo,
                                s_hi, s_lo, vlen);

  VFLOAT x2x_hi, x2x_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2 (x2_hi, x2_lo, x_hi, x_lo,
                                x2x_hi, x2x_lo, vlen);
  VFLOAT correction_hi, correction_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF (x2x_hi, x2x_lo, t,
                               correction_hi, correction_lo, vlen);
  RVV_SLEEF_ORIG_DFADD_VF2_VF2 (s_hi, s_lo, correction_hi, correction_lo,
                                *r_hi, *r_lo, vlen);
}

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT y = __riscv_vfabs (x, vlen);
  VFLOAT neg_y = VFNEG_V (y, vlen);

  VFLOAT num_hi, num_lo, den_hi, den_lo;
  RVV_SLEEF_ORIG_DFADD2_VF_VF (VFMV_VF (1.0f, vlen), y,
                               num_hi, num_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF_VF (VFMV_VF (1.0f, vlen), neg_y,
                               den_hi, den_lo, vlen);

  VFLOAT log_hi, log_lo;
  VFLOAT ratio_hi, ratio_lo;
  RVV_SLEEF_ORIG_DFDIV_VF2_VF2 (num_hi, num_lo, den_hi, den_lo,
                                ratio_hi, ratio_lo, vlen);
  RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_atanh_logk2f_m, LMUL)
    (ratio_hi, ratio_lo, &log_hi, &log_lo, vlen);

  VFLOAT finite = VFMUL_VF (VFADD_VV (log_hi, log_lo, vlen), 0.5f, vlen);
  VUINT ax = __riscv_vand (F_AS_U (y), 0x7fffffff, vlen);
  VBOOL y_gt_one = __riscv_vmsgtu (ax, 0x3f800000, vlen);
  VBOOL y_eq_one = __riscv_vmseq (ax, 0x3f800000, vlen);
  VFLOAT result = __riscv_vmerge (finite, VFMV_VF (fp_posInf, vlen),
                                  y_eq_one, vlen);

  VUINT all_ones = VMVU_VX (~(uint32_t) 0, vlen);
  VUINT domain_bits = __riscv_vmerge (VMVU_VX (0, vlen), all_ones,
                                      y_gt_one, vlen);
  result = U_AS_F (__riscv_vor (F_AS_U (result), domain_bits, vlen));

  VUINT ax_result = __riscv_vand (F_AS_U (result), 0x7fffffff, vlen);
  VBOOL result_nan = __riscv_vmsgtu (ax_result, 0x7f800000, vlen);
  VUINT inf_or_nan_bits = __riscv_vmerge (VMVU_VX (0, vlen), all_ones,
                                          __riscv_vmor (
                                            __riscv_vmseq (ax, 0x7f800000,
                                                           vlen),
                                            result_nan, vlen), vlen);
  result = U_AS_F (__riscv_vor (F_AS_U (result), inf_or_nan_bits, vlen));
  result = __riscv_vfsgnjx (result, x, vlen);

  VUINT ax_input = __riscv_vand (F_AS_U (x), 0x7fffffff, vlen);
  VBOOL x_nan = __riscv_vmsgtu (ax_input, 0x7f800000, vlen);
  VUINT input_nan_bits = __riscv_vmerge (VMVU_VX (0, vlen), all_ones,
                                         x_nan, vlen);
  return U_AS_F (__riscv_vor (F_AS_U (result), input_nan_bits, vlen));
}

#undef LMUL
#undef SIMDLEN
