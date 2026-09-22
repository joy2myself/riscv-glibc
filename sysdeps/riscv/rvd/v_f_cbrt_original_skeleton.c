/* Single-precision source-original SLEEF cbrtf skeleton.

   The finite path follows SLEEF xcbrtf_u1, including its explicit RTZ
   quotient/remainder conversion and double-float refinement order.  No
   clean-side EH quarantine, integer quotient repair, or FRM envelope is
   included here.  */

#include <vecmath_config.h>
#include <v_math.h>

#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

static inline VINT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_cbrt_vilogbk_m, LMUL)
  (VFLOAT d, size_t vlen)
{
  VBOOL small = __riscv_vmflt
    (d, VFMV_VF (RVV_SLEEF_CBRTF_VILOGBK_THRESH, vlen), vlen);
  VFLOAT scaled = VFMUL_VF (d, RVV_SLEEF_CBRTF_VILOGBK_SCALE, vlen);
  d = __riscv_vmerge (d, scaled, small, vlen);
  VUINT q = __riscv_vsrl (F_AS_U (d), MAN_LEN_F, vlen);
  q = __riscv_vand (q, RVV_SLEEF_CBRTF_VILOGBK_MASK, vlen);
  VINT bias = __riscv_vmerge
    (VMVI_VX (RVV_SLEEF_CBRTF_VILOGBK_SMALL_BIAS, vlen),
     VMVI_VX (RVV_SLEEF_CBRTF_VILOGBK_LARGE_BIAS, vlen), small, vlen);
  return __riscv_vsub (U_AS_I (q), bias, vlen);
}

static inline VFLOAT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_cbrt_vpow2i_m, LMUL)
  (VINT q, size_t vlen)
{
  VINT biased = __riscv_vadd (q, EXP_BIAS_F, vlen);
  return I_AS_F (__riscv_vsll (biased, MAN_LEN_F, vlen));
}

static inline VFLOAT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_cbrt_vldexp2_m, LMUL)
  (VFLOAT d, VINT e, size_t vlen)
{
  VINT e_half = __riscv_vsra (e, 1, vlen);
  VINT e_rest = __riscv_vsub (e, e_half, vlen);
  VFLOAT s1 = RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_cbrt_vpow2i_m, LMUL) (e_half, vlen);
  VFLOAT s2 = RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_cbrt_vpow2i_m, LMUL) (e_rest, vlen);
  return VFMUL_VV (VFMUL_VV (d, s1, vlen), s2, vlen);
}

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT d, size_t vlen)
{
  VFLOAT input = d;
  VINT e = __riscv_vadd (
    RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_cbrt_vilogbk_m, LMUL)
      (__riscv_vfabs (d, vlen), vlen),
    RVV_SLEEF_CBRTF_EXPONENT_OFFSET, vlen);
  d = RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_cbrt_vldexp2_m, LMUL)
    (d, __riscv_vneg (e, vlen), vlen);

  VFLOAT t = VFADD_VV (__riscv_vfcvt_f (e, vlen),
                       VFMV_VF (RVV_SLEEF_CBRTF_QUOTIENT_OFFSET_F, vlen),
                       vlen);
  VINT qu = __riscv_vfcvt_rtz_x
    (VFMUL_VF (t, RVV_SLEEF_CBRTF_ONE_THIRD, vlen), vlen);
  VFLOAT qu_f = __riscv_vfcvt_f (qu, vlen);
  VFLOAT rem = VFSUB_VV
    (t, VFMUL_VF (qu_f, RVV_SLEEF_CBRTF_QUOTIENT_DIVISOR_F, vlen), vlen);
  VINT re = __riscv_vfcvt_rtz_x (rem, vlen);

  VFLOAT q2_hi = VFMV_VF (fp_posOne_f, vlen);
  VFLOAT q2_lo = VFMV_VF (fp_posZero_f, vlen);
  VBOOL re_one = __riscv_vmseq (re, RVV_SLEEF_CBRTF_REMAINDER_ONE, vlen);
  VBOOL re_two = __riscv_vmseq (re, RVV_SLEEF_CBRTF_REMAINDER_TWO, vlen);
  q2_hi = __riscv_vmerge (q2_hi,
                          VFMV_VF (RVV_SLEEF_CBRTF_Q2_1_HI, vlen),
                          re_one, vlen);
  q2_lo = __riscv_vmerge (q2_lo,
                          VFMV_VF (RVV_SLEEF_CBRTF_Q2_1_LO, vlen),
                          re_one, vlen);
  q2_hi = __riscv_vmerge (q2_hi,
                          VFMV_VF (RVV_SLEEF_CBRTF_Q2_2_HI, vlen),
                          re_two, vlen);
  q2_lo = __riscv_vmerge (q2_lo,
                          VFMV_VF (RVV_SLEEF_CBRTF_Q2_2_LO, vlen),
                          re_two, vlen);

  q2_hi = __riscv_vfsgnjx (q2_hi, d, vlen);
  q2_lo = __riscv_vfsgnjx (q2_lo, d, vlen);
  d = __riscv_vfabs (d, vlen);

  VFLOAT x = VFMV_VF (RVV_SLEEF_CBRTF_POLY0_F, vlen);
  x = VFMADD_VV (x, d, VFMV_VF (RVV_SLEEF_CBRTF_POLY1_F, vlen), vlen);
  x = VFMADD_VV (x, d, VFMV_VF (RVV_SLEEF_CBRTF_POLY2_F, vlen), vlen);
  x = VFMADD_VV (x, d, VFMV_VF (RVV_SLEEF_CBRTF_POLY3_F, vlen), vlen);
  x = VFMADD_VV (x, d, VFMV_VF (RVV_SLEEF_CBRTF_POLY4_F, vlen), vlen);
  x = VFMADD_VV (x, d, VFMV_VF (RVV_SLEEF_CBRTF_POLY5_F, vlen), vlen);

  VFLOAT y = VFMUL_VV (x, x, vlen);
  y = VFMUL_VV (y, y, vlen);
  VFLOAT residual = VFNMSUB_VV (d, y, x, vlen);
  VFLOAT correction = VFMUL_VF
    (residual, RVV_SLEEF_CBRTF_NEG_ONE_THIRD, vlen);
  x = VFSUB_VV (x, correction, vlen);
  VFLOAT z = x;

  VFLOAT u_hi, u_lo;
  RVV_SLEEF_ORIG_DFMUL_VF_VF (x, x, u_hi, u_lo, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2 (u_hi, u_lo, u_hi, u_lo,
                                u_hi, u_lo, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF (u_hi, u_lo, d, u_hi, u_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF (u_hi, u_lo, VFNEG_V (x, vlen),
                                u_hi, u_lo, vlen);
  y = VFADD_VV (u_hi, u_lo, vlen);

  y = VFMUL_VF (y, RVV_SLEEF_CBRTF_NEG_TWO_THIRDS, vlen);
  y = VFMUL_VV (y, z, vlen);

  VFLOAT v_hi, v_lo;
  RVV_SLEEF_ORIG_DFMUL_VF_VF (z, z, v_hi, v_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF (v_hi, v_lo, y, v_hi, v_lo, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF (v_hi, v_lo, d, v_hi, v_lo, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2 (v_hi, v_lo, q2_hi, q2_lo,
                                v_hi, v_lo, vlen);

  VINT exp_adj = __riscv_vsub (qu, RVV_SLEEF_CBRTF_EXP_ADJUST, vlen);
  z = RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_cbrt_vldexp2_m, LMUL)
    (VFADD_VV (v_hi, v_lo, vlen), exp_adj, vlen);

  VUINT abs_input = __riscv_vand (F_AS_U (input), 0x7fffffff, vlen);
  VBOOL input_inf = __riscv_vmseq (abs_input, 0x7f800000, vlen);
  VBOOL input_zero = __riscv_vmseq (abs_input, 0, vlen);
  /* SLEEF vmulsign(+Inf, q2_hi): preserve the source sign-XOR primitive.  */
  VFLOAT inf_signed = __riscv_vfsgnjx
    (VFMV_VF (fp_posInf_f, vlen), q2_hi, vlen);
  VFLOAT zero_signed = __riscv_vfsgnj
    (VFMV_VF (fp_posZero_f, vlen), q2_hi, vlen);
  z = __riscv_vmerge (z, inf_signed, input_inf, vlen);
  z = __riscv_vmerge (z, zero_signed, input_zero, vlen);
  return z;
}

#undef LMUL
#undef SIMDLEN
