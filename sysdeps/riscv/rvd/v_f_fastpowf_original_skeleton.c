/* Single-precision source-original SLEEF fastpowf skeleton.

   Target-owned translation of SLEEF xfastpowf_u3500.  The source finite
   contract is kept in order: logk3f(abs(x)), multiply by y, expk3f, then the
   source integer/sign and zero identities.  Clean-only EH quarantine and
   special-result preparation are intentionally not imported into this
   comparison implementation.

   SLEEF is compiled with -ffp-contract=off; the source-original object uses
   the same build policy.  */

#include <vecmath_config.h>
#include <v_math.h>

#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif


#define RVV_SLEEF_ORIG_FASTPOWF_MANGLE(base) \
  RVV_SLEEF_ORIG_MANGLE (base, LMUL)
#define RVV_SLEEF_ORIG_FASTPOWF_VFCVT_X_RNE(value, vlen) \
  RVV_VFCVT_X_F_RNE ((value), (vlen))

#ifndef RVV_SLEEF_ORIG_FASTPOWF_CONSTANTS
# define RVV_SLEEF_ORIG_FASTPOWF_CONSTANTS
# define RVV_SLEEF_ORIG_FASTPOWF_R_LN2 \
  1.442695040888963407359924681001892137426645954152985934135449406931f
# define RVV_SLEEF_ORIG_FASTPOWF_L2U 0.693145751953125f
# define RVV_SLEEF_ORIG_FASTPOWF_L2L 1.428606765330187045e-06f
# define RVV_SLEEF_ORIG_FASTPOWF_LN2 \
  0.693147180559945286226764f
# define RVV_SLEEF_ORIG_FASTPOWF_INV_075 (1.0f / 0.75f)
# define RVV_SLEEF_ORIG_FASTPOWF_SCALE_SMALL 0x1.0p64f
# define RVV_SLEEF_ORIG_FASTPOWF_FLT_MIN 0x1.0p-126f
# define RVV_SLEEF_ORIG_FASTPOWF_INT_BOUND 16777216.0f
#endif

static inline VINT
RVV_SLEEF_ORIG_FASTPOWF_MANGLE (rvv_sleef_orig_fastpowf_vilogb2k)
  (VFLOAT d, size_t vlen)
{
  VUINT bits = F_AS_U (d);
  VUINT exponent = __riscv_vsrl (bits, MAN_LEN_F, vlen);
  exponent = __riscv_vand (exponent, 0xff, vlen);
  return __riscv_vsub (U_AS_I (exponent), EXP_BIAS_F, vlen);
}

static inline VFLOAT
RVV_SLEEF_ORIG_FASTPOWF_MANGLE (rvv_sleef_orig_fastpowf_vldexp3)
  (VFLOAT d, VINT q, size_t vlen)
{
  return I_AS_F (__riscv_vadd
                 (F_AS_I (d), __riscv_vsll (q, MAN_LEN_F, vlen), vlen));
}

static inline VFLOAT
RVV_SLEEF_ORIG_FASTPOWF_MANGLE (rvv_sleef_orig_fastpowf_vpow2i)
  (VINT q, size_t vlen)
{
  VINT biased = __riscv_vadd (q, EXP_BIAS_F, vlen);
  return I_AS_F (__riscv_vsll (biased, MAN_LEN_F, vlen));
}

static inline VFLOAT
RVV_SLEEF_ORIG_FASTPOWF_MANGLE (rvv_sleef_orig_fastpowf_vldexp2)
  (VFLOAT d, VINT e, size_t vlen)
{
  /* Source vldexp2: split the exponent and multiply both powers.  */
  VINT e_half = __riscv_vsra (e, 1, vlen);
  VINT e_rest = __riscv_vsub (e, e_half, vlen);
  VFLOAT p_half = RVV_SLEEF_ORIG_FASTPOWF_MANGLE
    (rvv_sleef_orig_fastpowf_vpow2i) (e_half, vlen);
  VFLOAT p_rest = RVV_SLEEF_ORIG_FASTPOWF_MANGLE
    (rvv_sleef_orig_fastpowf_vpow2i) (e_rest, vlen);
  return VFMUL_VV (VFMUL_VV (d, p_half, vlen), p_rest, vlen);
}

static inline VFLOAT
RVV_SLEEF_ORIG_FASTPOWF_MANGLE (rvv_sleef_orig_fastpowf_vtruncate)
  (VFLOAT d, size_t vlen)
{
  /* helperrvv.h's vtruncate_vf_vf is RTZ conversion followed by float cast.  */
  VINT truncated = __riscv_vfcvt_rtz_x (d, vlen);
  return __riscv_vfcvt_f (truncated, vlen);
}

static inline VFLOAT
RVV_SLEEF_ORIG_FASTPOWF_MANGLE (rvv_sleef_orig_fastpowf_logk3f)
  (VFLOAT d, size_t vlen)
{
  VFLOAT x, x2, t, m;

  VBOOL small = __riscv_vmflt
    (d, RVV_SLEEF_ORIG_FASTPOWF_FLT_MIN, vlen);
  VFLOAT scaled_small = VFMUL_VF
    (d, RVV_SLEEF_ORIG_FASTPOWF_SCALE_SMALL, vlen);
  d = VMERGE_VVM_F (d, scaled_small, small, vlen);

  VFLOAT scaled = VFMUL_VF
    (d, RVV_SLEEF_ORIG_FASTPOWF_INV_075, vlen);
  VINT e = RVV_SLEEF_ORIG_FASTPOWF_MANGLE
    (rvv_sleef_orig_fastpowf_vilogb2k) (scaled, vlen);
  m = RVV_SLEEF_ORIG_FASTPOWF_MANGLE
    (rvv_sleef_orig_fastpowf_vldexp3) (d, __riscv_vneg (e, vlen), vlen);
  e = VMERGE_VVM_I
    (e, __riscv_vsub (e, 64, vlen), small, vlen);

  x = VFDIV_VV
    (VFSUB_VV (m, VFMV_VF (1.0f, vlen), vlen),
     VFADD_VV (VFMV_VF (1.0f, vlen), m, vlen), vlen);
  x2 = VFMUL_VV (x, x, vlen);

  t = VFMV_VF (0.2392828464508056640625f, vlen);
  t = VFMADD_VV
    (t, x2, VFMV_VF (0.28518211841583251953125f, vlen), vlen);
  t = VFMADD_VV
    (t, x2, VFMV_VF (0.400005877017974853515625f, vlen), vlen);
  t = VFMADD_VV
    (t, x2, VFMV_VF (0.666666686534881591796875f, vlen), vlen);
  t = VFMADD_VV (t, x2, VFMV_VF (2.0f, vlen), vlen);

  VFLOAT e_f = __riscv_vfcvt_f (e, vlen);
  VFLOAT ln2e = VFMUL_VF
    (e_f, RVV_SLEEF_ORIG_FASTPOWF_LN2, vlen);
  return VFMADD_VV (x, t, ln2e, vlen);
}

static inline VFLOAT
RVV_SLEEF_ORIG_FASTPOWF_MANGLE (rvv_sleef_orig_fastpowf_expk3f)
  (VFLOAT d, size_t vlen)
{
  VINT q = RVV_SLEEF_ORIG_FASTPOWF_VFCVT_X_RNE
    (VFMUL_VF (d, RVV_SLEEF_ORIG_FASTPOWF_R_LN2, vlen), vlen);
  VFLOAT q_f = __riscv_vfcvt_f (q, vlen);
  VFLOAT s = VFMADD_VV
    (q_f, VFMV_VF (-RVV_SLEEF_ORIG_FASTPOWF_L2U, vlen), d, vlen);
  s = VFMADD_VV
    (q_f, VFMV_VF (-RVV_SLEEF_ORIG_FASTPOWF_L2L, vlen), s, vlen);

  VFLOAT u = VFMV_VF (0.000198527617612853646278381f, vlen);
  u = VFMADD_VV
    (u, s, VFMV_VF (0.00139304355252534151077271f, vlen), vlen);
  u = VFMADD_VV
    (u, s, VFMV_VF (0.00833336077630519866943359f, vlen), vlen);
  u = VFMADD_VV
    (u, s, VFMV_VF (0.0416664853692054748535156f, vlen), vlen);
  u = VFMADD_VV
    (u, s, VFMV_VF (0.166666671633720397949219f, vlen), vlen);
  u = VFMADD_VV (u, s, VFMV_VF (0.5f, vlen), vlen);

  VFLOAT s2 = VFMUL_VV (s, s, vlen);
  u = VFMADD_VV
    (s2, u, VFADD_VV (s, VFMV_VF (1.0f, vlen), vlen), vlen);
  u = RVV_SLEEF_ORIG_FASTPOWF_MANGLE
    (rvv_sleef_orig_fastpowf_vldexp2) (u, q, vlen);

  VBOOL too_small = __riscv_vmflt
    (d, -104.0f, vlen);
  VUINT keep = __riscv_vmerge
    (VMVU_VX (~0u, vlen), VMVU_VX (0, vlen), too_small, vlen);
  return U_AS_F (__riscv_vand (F_AS_U (u), keep, vlen));
}

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL)
  (VFLOAT x, VFLOAT y, size_t vlen)
{
  /* xfastpowf_u3500 source order: compute the finite result first.  The
     source's vrint_vi2_vf is an explicit RNE conversion inside expk3f;
     surrounding floating-point arithmetic retains the caller's FRM.  */
  VFLOAT log_abs_x = RVV_SLEEF_ORIG_FASTPOWF_MANGLE
    (rvv_sleef_orig_fastpowf_logk3f) (__riscv_vfabs (x, vlen), vlen);
  VFLOAT result = RVV_SLEEF_ORIG_FASTPOWF_MANGLE
    (rvv_sleef_orig_fastpowf_expk3f) (VFMUL_VV (log_abs_x, y, vlen), vlen);

  VFLOAT abs_y = __riscv_vfabs (y, vlen);
  VFLOAT y_trunc = RVV_SLEEF_ORIG_FASTPOWF_MANGLE
    (rvv_sleef_orig_fastpowf_vtruncate) (y, vlen);
  VBOOL y_is_int = __riscv_vmor
    (__riscv_vmfeq (y_trunc, y, vlen),
     __riscv_vmfgt (abs_y, RVV_SLEEF_ORIG_FASTPOWF_INT_BOUND, vlen), vlen);
  VINT y_trunc_i = __riscv_vfcvt_rtz_x (y, vlen);
  VBOOL y_is_odd = __riscv_vmand
    (__riscv_vmseq (__riscv_vand (y_trunc_i, 1, vlen), 1, vlen),
     y_is_int, vlen);
  y_is_odd = __riscv_vmand
    (y_is_odd, __riscv_vmflt
     (abs_y, RVV_SLEEF_ORIG_FASTPOWF_INT_BOUND, vlen), vlen);

  VUINT x_sign = __riscv_vand (F_AS_U (x), 0x80000000u, vlen);
  VBOOL x_negative = __riscv_vmseq (U_AS_I (x_sign), 0x80000000, vlen);
  VBOOL sign_flip = __riscv_vmand (x_negative, y_is_odd, vlen);
  VUINT sign_mask = __riscv_vmerge
    (VMVU_VX (0, vlen), VMVU_VX (0x80000000u, vlen), sign_flip, vlen);
  result = U_AS_F (__riscv_vxor (F_AS_U (result), sign_mask, vlen));

  result = VFMERGE_VFM
    (result, fp_posZero_f, __riscv_vmfeq (x, 0.0f, vlen), vlen);
  result = VFMERGE_VFM
    (result, fp_posOne_f, __riscv_vmfeq (y, 0.0f, vlen), vlen);

  return result;
}

#undef RVV_SLEEF_ORIG_FASTPOWF_MANGLE
#undef RVV_SLEEF_ORIG_FASTPOWF_VFCVT_X_RNE
#undef RVV_SLEEF_ORIG_FASTPOWF_INT_BOUND
#undef RVV_SLEEF_ORIG_FASTPOWF_FLT_MIN
#undef RVV_SLEEF_ORIG_FASTPOWF_SCALE_SMALL
#undef RVV_SLEEF_ORIG_FASTPOWF_INV_075
#undef RVV_SLEEF_ORIG_FASTPOWF_LN2
#undef RVV_SLEEF_ORIG_FASTPOWF_L2L
#undef RVV_SLEEF_ORIG_FASTPOWF_L2U
#undef RVV_SLEEF_ORIG_FASTPOWF_R_LN2
#undef RVV_SLEEF_ORIG_FASTPOWF_CONSTANTS
#undef LMUL
#undef SIMDLEN
