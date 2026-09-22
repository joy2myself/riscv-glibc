/* Single-precision source-original SLEEF expf skeleton.

   The finite kernel is a target-owned translation of SLEEF xexpf.  It keeps
   the source split-ln(2) reduction, ordinary-float polynomial, and vldexp2
   reconstruction.  Clean-only EH quarantine and endpoint repairs are not
   silently imported into this comparison implementation.  */

#include <vecmath_config.h>
#include <v_math.h>

#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (expf_original_vpow2i, LMUL) (VINT e, size_t vlen)
{
  VINT biased = __riscv_vadd (e, EXP_BIAS_F, vlen);
  return I_AS_F (__riscv_vsll (biased, MAN_LEN_F, vlen));
}

static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (expf_original_vldexp2, LMUL)
  (VFLOAT d, VINT e, size_t vlen)
{
  /* Source vldexp2: split e into two powers and multiply twice.  */
  VINT e_half = __riscv_vsra (e, 1, vlen);
  VINT e_rest = __riscv_vsub (e, e_half, vlen);
  VFLOAT p_half = RVV_LIBMVEC_IMPL_NAME
    (expf_original_vpow2i, LMUL) (e_half, vlen);
  VFLOAT p_rest = RVV_LIBMVEC_IMPL_NAME
    (expf_original_vpow2i, LMUL) (e_rest, vlen);
  return VFMUL_VV (VFMUL_VV (d, p_half, vlen), p_rest, vlen);
}

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT d, size_t vlen)
{
  VFLOAT scaled = VFMUL_VF (d, RVV_EXPF_SOURCE_R_LN2_F, vlen);
  VINT q = RVV_VFCVT_X_F_RNE (scaled, vlen);
  VFLOAT q_f = __riscv_vfcvt_f (q, vlen);

  VFLOAT s = VFMADD_VV (q_f,
                        VFMV_VF (-RVV_EXPF_SOURCE_L2U_F, vlen), d,
                        vlen);
  s = VFMADD_VV (q_f,
                 VFMV_VF (-RVV_EXPF_SOURCE_L2L_F, vlen), s, vlen);

  VFLOAT u = VFMV_VF
    (RVV_SLEEF_ORIG_F64_TO_F32 (RVV_EXPF_SOURCE_POLY0_D), vlen);
  u = VFMADD_VV (u, s, VFMV_VF
    (RVV_SLEEF_ORIG_F64_TO_F32 (RVV_EXPF_SOURCE_POLY1_D), vlen),
                 vlen);
  u = VFMADD_VV (u, s, VFMV_VF
    (RVV_SLEEF_ORIG_F64_TO_F32 (RVV_EXPF_SOURCE_POLY2_D), vlen),
                 vlen);
  u = VFMADD_VV (u, s, VFMV_VF
    (RVV_SLEEF_ORIG_F64_TO_F32 (RVV_EXPF_SOURCE_POLY3_D), vlen),
                 vlen);
  u = VFMADD_VV (u, s, VFMV_VF
    (RVV_SLEEF_ORIG_F64_TO_F32 (RVV_EXPF_SOURCE_POLY4_D), vlen),
                 vlen);
  u = VFMADD_VV (u, s, VFMV_VF
    (RVV_SLEEF_ORIG_F64_TO_F32 (RVV_EXPF_SOURCE_POLY5_D), vlen), vlen);

  VFLOAT s2 = VFMUL_VV (s, s, vlen);
  u = VFMADD_VV (s2, u, s, vlen);
  u = VFADD_VV (VFMV_VF (1.0f, vlen), u, vlen);
  u = RVV_LIBMVEC_IMPL_NAME (expf_original_vldexp2, LMUL) (u, q, vlen);

  VBOOL underflow = __riscv_vmflt (d, RVV_EXPF_SOURCE_MIN_INPUT_F,
                                   vlen);
  VBOOL overflow = __riscv_vmfgt (d, RVV_EXPF_SOURCE_MAX_INPUT_F,
                                  vlen);
  u = VFMERGE_VFM (u, fp_posZero_f, underflow, vlen);
  u = VFMERGE_VFM (u, fp_posInf_f, overflow, vlen);

  return u;
}

#undef LMUL
#undef SIMDLEN
