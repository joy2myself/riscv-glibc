/* Single-precision source-original SLEEF asinhf skeleton.

   This is a target-owned translation of SLEEF xasinhf/logk2f from the pinned
   source checkout.  The finite double-float operations, source constants,
   reduction, and late special-value selections are kept separate from the
   clean implementation's early EH quarantine and large-input repair.  Only
   RVV types, LMUL/VL plumbing, and the glibc importer ABI are adapted; SLEEF
   source files are evidence only and are not included, compiled, or linked.

   The source build uses -ffp-contract=off.  The corresponding flag is
   attached to this translation unit in the RVV mathvec Makefile.  */

#include <stdint.h>
#include <vecmath_config.h>
#include <v_math.h>


#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

/* SLEEF vilogbk_vi2_vf, including the source subnormal normalization.  The
   logk2f input is positive in the finite xasinhf path, but the helper is kept
   source-shaped so the late source fixups see the same intermediate form.  */
static inline VINT
RVV_SLEEF_ORIG_ASINHF_MANGLE (rvv_sleef_orig_asinhf_vilogbk_m, LMUL)
  (VFLOAT d, size_t vlen)
{
  VBOOL small = __riscv_vmflt
    (d, 5.421010862427522e-20f, vlen);
  VFLOAT scaled = VFMUL_VF
    (d, 1.8446744073709552e19f, vlen);
  d = VMERGE_VVM_F (d, scaled, small, vlen);

  VUINT q = __riscv_vsrl (F_AS_U (d), MAN_LEN_F, vlen);
  q = __riscv_vand (q, 0xff, vlen);
  VINT bias = VMERGE_VVM_I
    (VMVI_VX (EXP_BIAS_F, vlen),
     VMVI_VX (EXP_BIAS_F + 64, vlen), small, vlen);
  return __riscv_vsub (U_AS_I (q), bias, vlen);
}

/* SLEEF vpow2i_vf_vi2: construct an exact binary32 power of two by bits.  */
static inline VFLOAT
RVV_SLEEF_ORIG_ASINHF_MANGLE (rvv_sleef_orig_asinhf_vpow2i_m, LMUL)
  (VINT q, size_t vlen)
{
  VINT biased = __riscv_vadd (q, EXP_BIAS_F, vlen);
  return I_AS_F (__riscv_vsll (biased, MAN_LEN_F, vlen));
}

/* SLEEF dfsqrt_vf2_vf2, using the pinned non-reciprocal-sqrt source branch:
     t = sqrt(d.x + d.y)
     q = dfmul(dfadd2(d, dfmul(t,t)), dfrec(t))
     result = q * 0.5.  */
static inline void
RVV_SLEEF_ORIG_ASINHF_MANGLE (rvv_sleef_orig_asinhf_dfsqrt_m, LMUL)
  (VFLOAT d_hi, VFLOAT d_lo, VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT t = __riscv_vfsqrt (VFADD_VV (d_hi, d_lo, vlen), vlen);

  VFLOAT tsq_hi, tsq_lo;
  RVV_SLEEF_ORIG_DFMUL_VF_VF (t, t, tsq_hi, tsq_lo, vlen);

  VFLOAT sum_hi, sum_lo;
  RVV_SLEEF_ORIG_DFADD2_VF2_VF2
    (d_hi, d_lo, tsq_hi, tsq_lo, sum_hi, sum_lo, vlen);

  VFLOAT rec_hi, rec_lo;
  VFLOAT one = VFMV_VF (1.0f, vlen);
  rec_hi = VFDIV_VV (one, t, vlen);
  VFLOAT correction = VFNMSUB_VV (t, rec_hi, one, vlen);
  rec_lo = VFMUL_VV (rec_hi, correction, vlen);

  VFLOAT q_hi, q_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (sum_hi, sum_lo, rec_hi, rec_lo, q_hi, q_lo, vlen);

  *r_hi = VFMUL_VF (q_hi, 0.5f, vlen);
  *r_lo = VFMUL_VF (q_lo, 0.5f, vlen);
}

/* SLEEF logk2f, preserving its source reduction, polynomial, and
   double-float accumulation order.  */
static inline void
RVV_SLEEF_ORIG_ASINHF_MANGLE (rvv_sleef_orig_asinhf_logk2f_m, LMUL)
  (VFLOAT d_hi, VFLOAT d_lo, VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT d_adj = VFMUL_VF (d_hi, 1.0f / 0.75f, vlen);
  VINT e = RVV_SLEEF_ORIG_ASINHF_MANGLE
    (rvv_sleef_orig_asinhf_vilogbk_m, LMUL) (d_adj, vlen);

  VFLOAT scale = RVV_SLEEF_ORIG_ASINHF_MANGLE
    (rvv_sleef_orig_asinhf_vpow2i_m, LMUL)
    (__riscv_vneg (e, vlen), vlen);
  VFLOAT m_hi, m_lo;
  RVV_SLEEF_ORIG_DFSCALE_VF2_VF
    (d_hi, d_lo, scale, m_hi, m_lo, vlen);

  VFLOAT num_hi, num_lo;
  RVV_SLEEF_ORIG_DFADD2_VF2_VF_S
    (m_hi, m_lo, -1.0f, num_hi, num_lo, vlen);
  VFLOAT den_hi, den_lo;
  RVV_SLEEF_ORIG_DFADD2_VF2_VF_S
    (m_hi, m_lo, 1.0f, den_hi, den_lo, vlen);

  VFLOAT x_hi, x_lo;
  RVV_SLEEF_ORIG_DFDIV_VF2_VF2
    (num_hi, num_lo, den_hi, den_lo, x_hi, x_lo, vlen);

  VFLOAT x2_hi, x2_lo;
  RVV_SLEEF_ORIG_DFSQU_VF2 (x_hi, x_lo, x2_hi, x2_lo, vlen);

  VFLOAT t = VFMV_VF (RVV_SLEEF_ORIG_ASINHF_LOG_P0, vlen);
  t = VFMADD_VV
    (t, x2_hi, VFMV_VF (RVV_SLEEF_ORIG_ASINHF_LOG_P1, vlen), vlen);
  t = VFMADD_VV
    (t, x2_hi, VFMV_VF (RVV_SLEEF_ORIG_ASINHF_LOG_P2, vlen), vlen);
  t = VFMADD_VV
    (t, x2_hi, VFMV_VF (RVV_SLEEF_ORIG_ASINHF_LOG_P3, vlen), vlen);

  VFLOAT s_hi, s_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (VFMV_VF (RVV_SLEEF_ORIG_ASINHF_L2U, vlen),
     VFMV_VF (RVV_SLEEF_ORIG_ASINHF_L2L, vlen),
     __riscv_vfcvt_f (e, vlen), s_hi, s_lo, vlen);

  VFLOAT two_x_hi, two_x_lo;
  RVV_SLEEF_ORIG_DFSCALE_VF2_VF
    (x_hi, x_lo, VFMV_VF (2.0f, vlen), two_x_hi, two_x_lo, vlen);
  RVV_SLEEF_ORIG_DFADD_VF2_VF2
    (s_hi, s_lo, two_x_hi, two_x_lo, s_hi, s_lo, vlen);

  VFLOAT x2x_hi, x2x_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (x2_hi, x2_lo, x_hi, x_lo, x2x_hi, x2x_lo, vlen);
  VFLOAT tail_hi, tail_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (x2x_hi, x2x_lo, t, tail_hi, tail_lo, vlen);
  RVV_SLEEF_ORIG_DFADD_VF2_VF2
    (s_hi, s_lo, tail_hi, tail_lo, *r_hi, *r_lo, vlen);
}

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  /* xasinhf source body: form d from either 1/x or |x|, take the pair
     square root, then add x before logk2f.  Special values are intentionally
     left to the source's late selection sequence.  */
  VFLOAT y = __riscv_vfabs (x, vlen);
  VBOOL use_recip = __riscv_vmfgt (y, 1.0f, vlen);

  VFLOAT rec_hi, rec_lo;
  VFLOAT one = VFMV_VF (1.0f, vlen);
  rec_hi = VFDIV_VV (one, x, vlen);
  VFLOAT correction = VFNMSUB_VV (x, rec_hi, one, vlen);
  rec_lo = VFMUL_VV (rec_hi, correction, vlen);
  VFLOAT d_hi = VMERGE_VVM_F (y, rec_hi, use_recip, vlen);
  VFLOAT d_lo = VMERGE_VVM_F
    (VFMV_VF (0.0f, vlen), rec_lo, use_recip, vlen);

  /* The source calls dfsqrt(dfadd2(dfsqu(d), 1)), rather than taking the
     square root of d itself.  Keep the pair square/add order explicit here;
     dropping it changes the finite function, not merely its rounding.  */
  VFLOAT d_sq_hi, d_sq_lo;
  RVV_SLEEF_ORIG_DFSQU_VF2
    (d_hi, d_lo, d_sq_hi, d_sq_lo, vlen);

  VFLOAT d_one_hi, d_one_lo;
  RVV_SLEEF_ORIG_DFADD2_VF2_VF_S
    (d_sq_hi, d_sq_lo, 1.0f, d_one_hi, d_one_lo, vlen);

  VFLOAT sqrt_hi, sqrt_lo;
  RVV_SLEEF_ORIG_ASINHF_MANGLE (rvv_sleef_orig_asinhf_dfsqrt_m, LMUL)
    (d_one_hi, d_one_lo, &sqrt_hi, &sqrt_lo, vlen);

  VFLOAT scaled_hi, scaled_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (sqrt_hi, sqrt_lo, y, scaled_hi, scaled_lo, vlen);
  d_hi = VMERGE_VVM_F (sqrt_hi, scaled_hi, use_recip, vlen);
  d_lo = VMERGE_VVM_F (sqrt_lo, scaled_lo, use_recip, vlen);

  VFLOAT arg_hi, arg_lo;
  RVV_SLEEF_ORIG_DFADD2_VF2_VF
    (d_hi, d_lo, x, arg_hi, arg_lo, vlen);
  VFLOAT norm_hi, norm_lo;
  RVV_SLEEF_ORIG_DFNORMALIZE_VF2
    (arg_hi, arg_lo, norm_hi, norm_lo, vlen);

  VFLOAT log_hi, log_lo;
  RVV_SLEEF_ORIG_ASINHF_MANGLE (rvv_sleef_orig_asinhf_logk2f_m, LMUL)
    (norm_hi, norm_lo, &log_hi, &log_lo, vlen);
  VFLOAT result = VFADD_VV (log_hi, log_lo, vlen);

  /* xasinhf late source selects: large/result-NaN -> signed Inf, then input
     NaN bitwise propagation, then the negative-zero repair.  */
  VBOOL y_nan = __riscv_vmfne (result, result, vlen);
  VBOOL x_large = __riscv_vmfgt
    (y, RVV_SLEEF_ORIG_ASINHF_SQRT_FLT_MAX, vlen);
  VBOOL inf_result = __riscv_vmor (x_large, y_nan, vlen);
  /* SLEEF vmulsign(+Inf, x): preserve the source sign-XOR primitive.  */
  VFLOAT signed_inf = __riscv_vfsgnjx
    (VFMV_VF (fp_posInf, vlen), x, vlen);
  result = __riscv_vmerge (result, signed_inf, inf_result, vlen);

  VBOOL x_nan = __riscv_vmfne (x, x, vlen);
  VUINT nan_bits = __riscv_vmerge
    (VMVU_VX (0, vlen), VMVU_VX (~(uint32_t) 0, vlen), x_nan, vlen);
  result = U_AS_F (__riscv_vor (F_AS_U (result), nan_bits, vlen));

  VBOOL negzero = __riscv_vmseq (F_AS_U (x), 0x80000000u, vlen);
  result = __riscv_vfmerge (result, fp_negZero, negzero, vlen);
  return result;
}

#undef LMUL
#undef SIMDLEN
