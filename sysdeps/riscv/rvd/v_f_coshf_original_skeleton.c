/* Single-precision source-original SLEEF coshf skeleton.

   This is a target-owned translation of SLEEF expk2f and xcoshf.  The
   source double-float reduction, polynomial, reconstruction, scaling, and
   late special-result selections remain separate from clean's EH path.
   Only the RVV representation, LMUL/VL plumbing, and glibc importer ABI are
   adapted; no SLEEF source file is included or linked.

   SLEEF requests -ffp-contract=off for this source.  The corresponding
   build flag is attached to the wrapper in the RVV mathvec Makefile.  */

#include <vecmath_config.h>
#include <v_math.h>

#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

/* SLEEF vpow2i_vf_vi2.  The expk2f quotient range keeps the source bit
   construction valid for this comparison implementation.  */
static inline VFLOAT
RVV_SLEEF_ORIG_COSHF_VPOW2I (VINT e, size_t vlen)
{
  VINT biased = __riscv_vadd (e, EXP_BIAS_F, vlen);
  return I_AS_F (__riscv_vsll (biased, MAN_LEN_F, vlen));
}

/* SLEEF vldexp2_vf_vf_vi2: split the exponent and preserve its two
   source-ordered power-of-two multiplies.  */
static inline VFLOAT
RVV_SLEEF_ORIG_COSHF_VLDEXP2 (VFLOAT d, VINT e, size_t vlen)
{
  VINT e_half = __riscv_vsra (e, 1, vlen);
  VINT e_rest = __riscv_vsub (e, e_half, vlen);
  VFLOAT p_half = RVV_SLEEF_ORIG_COSHF_VPOW2I (e_half, vlen);
  VFLOAT p_rest = RVV_SLEEF_ORIG_COSHF_VPOW2I (e_rest, vlen);
  return VFMUL_VV (VFMUL_VV (d, p_half, vlen), p_rest, vlen);
}

/* Source expk2f, including its double-float reconstruction order.  */
static inline void
RVV_SLEEF_ORIG_COSHF_EXPK2F (VFLOAT d_hi, VFLOAT d_lo,
                             VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT u = VFMUL_VF (VFADD_VV (d_hi, d_lo, vlen),
                       RVV_SLEEF_EXPM1F_R_LN2, vlen);
  VINT q = RVV_VFCVT_X_F_RNE (u, vlen);
  VFLOAT q_f = __riscv_vfcvt_f (q, vlen);

  VFLOAT s_hi, s_lo;
  VFLOAT q_l2u = VFMUL_VF (q_f, -RVV_SLEEF_EXPM1F_L2U, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF (d_hi, d_lo, q_l2u,
                                s_hi, s_lo, vlen);
  VFLOAT q_l2l = VFMUL_VF (q_f, -RVV_SLEEF_EXPM1F_L2L, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF (s_hi, s_lo, q_l2l,
                                s_hi, s_lo, vlen);

  u = VFMV_VF (RVV_SLEEF_EXPM1F_POLY0, vlen);
  u = VFMADD_VV (u, s_hi, VFMV_VF (RVV_SLEEF_EXPM1F_POLY1, vlen), vlen);
  u = VFMADD_VV (u, s_hi, VFMV_VF (RVV_SLEEF_EXPM1F_POLY2, vlen), vlen);
  u = VFMADD_VV (u, s_hi, VFMV_VF (RVV_SLEEF_EXPM1F_POLY3, vlen), vlen);

  VFLOAT t_hi, t_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF (s_hi, s_lo, u, t_hi, t_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF_S
    (t_hi, t_lo, RVV_SLEEF_EXPM1F_TAYLOR_C,
     t_hi, t_lo, vlen);

  VFLOAT tmp_hi, tmp_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (s_hi, s_lo, t_hi, t_lo, tmp_hi, tmp_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF_S
    (tmp_hi, tmp_lo, RVV_SLEEF_EXPM1F_HALF, t_hi, t_lo, vlen);

  VFLOAT s2_hi, s2_lo;
  RVV_SLEEF_ORIG_DFSQU_VF2 (s_hi, s_lo, s2_hi, s2_lo, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (s2_hi, s2_lo, t_hi, t_lo, tmp_hi, tmp_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF2
    (s_hi, s_lo, tmp_hi, tmp_lo, t_hi, t_lo, vlen);

  /* Source dfadd_vf2_vf_vf2 is the non-compensated scalar-plus-pair add.  */
  RVV_SLEEF_ORIG_DFADD_VF_VF2
    (VFMV_VF (1.0f, vlen), t_hi, t_lo, t_hi, t_lo, vlen);

  t_hi = RVV_SLEEF_ORIG_COSHF_VLDEXP2 (t_hi, q, vlen);
  t_lo = RVV_SLEEF_ORIG_COSHF_VLDEXP2 (t_lo, q, vlen);

  VBOOL underflow = __riscv_vmflt (d_hi, RVV_SLEEF_EXPM1F_UNDERFLOW, vlen);
  t_hi = VFMERGE_VFM (t_hi, fp_posZero_f, underflow, vlen);
  t_lo = VFMERGE_VFM (t_lo, fp_posZero_f, underflow, vlen);

  *r_hi = t_hi;
  *r_lo = t_lo;
}

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT y = __riscv_vfabs (x, vlen);
  VFLOAT d_hi, d_lo;
  RVV_SLEEF_ORIG_COSHF_EXPK2F
    (y, VFMV_VF (0.0f, vlen), &d_hi, &d_lo, vlen);

  VFLOAT e_hi, e_lo;
  RVV_SLEEF_ORIG_DFREC_VF2_VF2
    (d_hi, d_lo, e_hi, e_lo, vlen);

  /* xcoshf uses source dfadd_vf2_vf2_vf2, not compensated dfadd2.  */
  VFLOAT sum_hi, sum_lo;
  RVV_SLEEF_ORIG_DFADD_VF2_VF2
    (d_hi, d_lo, e_hi, e_lo, sum_hi, sum_lo, vlen);
  y = VFMUL_VF
    (VFADD_VV (sum_hi, sum_lo, vlen), RVV_SLEEF_EXPM1F_HALF, vlen);

  VBOOL result_special = __riscv_vmor
    (__riscv_vmfgt (__riscv_vfabs (x, vlen),
                    RVV_SLEEF_COSHF_SOURCE_THRESHOLD, vlen),
     __riscv_vmfne (y, y, vlen), vlen);
  y = VFMERGE_VFM (y, fp_posInf_f, result_special, vlen);

  /* Source's final OR preserves its NaN payload/quieting behavior.  */
  VBOOL x_nan = __riscv_vmfne (x, x, vlen);
  VINT nan_bits = VMERGE_VXM_I
    (VMVI_VX (0, vlen), -1, x_nan, vlen);
  y = I_AS_F
    (__riscv_vor (F_AS_I (y), nan_bits, vlen));

  return y;
}

#undef LMUL
#undef SIMDLEN
