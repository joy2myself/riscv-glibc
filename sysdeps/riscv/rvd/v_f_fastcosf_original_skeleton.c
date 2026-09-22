/* Single-precision source-original SLEEF fastcosf skeleton.

   This is a target-owned translation of SLEEF xfastcosf_u3500.  It keeps the
   source phase reduction, reduced-precision polynomial, quadrant sign
   reconstruction, and late xcosf fallback.  The fallback is instantiated
   from a private, target-owned translation of the pinned xcosf body in this
   same translation unit; it is not a call to clean cosf, xcosf_u1, or a
   SLEEF object/source.

   Clean-only quarantine, small-domain replacement, and fallback ordering are
   intentionally not imported into this comparison implementation.

   SLEEF is compiled with -ffp-contract=off; the source-original object is
   built with the same policy.  */

#include <vecmath_config.h>
#include <v_math.h>

#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

/* xfastcosf_u3500 calls xcosf for lanes outside its fast range.  Import only
   the target-owned rempif/DF helper definitions from the xcosf_u1
   translation; the xcosf_u35 fallback body below is separate because SLEEF's
   xcosf and xcosf_u1 are different source routines.  The embedded cosf
   skeleton selects helper-only mode from V_MATH_FUNC in v_math.h.  */
#include "v_f_cosf_original_skeleton.c"

/* Pinned helperrvv.h defines vrint_vf_vf as an explicit-RNE conversion to an
   integer followed by the source's ordinary integer-to-float conversion.
   Keep that operation local to the conversion itself; source constants and
   ordinary arithmetic must continue to observe the caller's FRM.  */
static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (rvv_sleef_orig_fastcosf_vrint, LMUL)
  (VFLOAT x, size_t vlen)
{
  VINT rounded = RVV_VFCVT_X_F_RNE (x, vlen);
  return __riscv_vfcvt_f (rounded, vlen);
}

/* Target-owned translation of pinned SLEEF xcosf/Sleef_cosf_u35.  The
   source's xcosf_u35 fallback has a different reduction and polynomial from
   xcosf_u1/Sleef_cosf_u10, so only the shared rempif helpers are reused.  */
static VFLOAT
RVV_LIBMVEC_IMPL_NAME (rvv_sleef_orig_fastcosf_xcosf, LMUL)
  (VFLOAT d, size_t vlen)
{
  VINT q;
  VFLOAT u, s, r = d;

  /* xcosf: q = rint(d * (1/pi) - 1/2), then q = 2*q + 1. */
  VFLOAT q_input = d;
#if RVV_MATH_IS_FASTCOSF_CLEAN_EH (V_MATH_FUNC)
  /* The source computes q before its late rempif selection.  Keep that
     source ordering for safe finite lanes, but quarantine the large lanes
     whose int32 conversion is outside the safe range.  The large branch
     below overwrites q and d from rempif, so zero is only a provider input. */
  VBOOL q_safe = __riscv_vmflt
    (__riscv_vfabs (r, vlen), RVV_FASTCOSF_XCOSF_TRIG_RANGE_MAX_F, vlen);
  q_input = VFMERGE_VFM
    (q_input, fp_posZero_f, __riscv_vmnot (q_safe, vlen), vlen);
#endif
  VFLOAT reduced = VFSUB_VF
    (VFMUL_VF (q_input, RVV_SLEEF_TRIG_ONE_OVER_PI, vlen),
     0.5f, vlen);
  q = RVV_VFCVT_X_F_RNE (reduced, vlen);
  q = __riscv_vadd (__riscv_vadd (q, q, vlen), 1, vlen);
  u = __riscv_vfcvt_f (q, vlen);

  /* Three-stage split pi/2 subtraction, as in pinned xcosf. */
  d = VFMADD_VV
    (u, VFMV_VF (RVV_SLEEF_COSF_NEG_PI_A2_HALF, vlen), d, vlen);
  d = VFMADD_VV
    (u, VFMV_VF (RVV_SLEEF_COSF_NEG_PI_B2_HALF, vlen), d, vlen);
  d = VFMADD_VV
    (u, VFMV_VF (RVV_SLEEF_COSF_NEG_PI_C2_HALF, vlen), d, vlen);

  VBOOL g = __riscv_vmflt
    (__riscv_vfabs (r, vlen), RVV_SLEEF_TRIG_RANGE_MAX2,
                    vlen);
  if (__riscv_vcpop (g, vlen) != vlen)
    {
      /* Medium range: use the source's four-part pi split. */
      s = __riscv_vfcvt_f (q, vlen);
      u = VFMADD_VV
        (s, VFMV_VF (RVV_FASTCOSF_XCOSF_NEG_PI_A_HALF, vlen), r,
         vlen);
      u = VFMADD_VV
        (s, VFMV_VF (RVV_FASTCOSF_XCOSF_NEG_PI_B_HALF, vlen), u,
         vlen);
      u = VFMADD_VV
        (s, VFMV_VF (RVV_FASTCOSF_XCOSF_NEG_PI_C_HALF, vlen), u,
         vlen);
      u = VFMADD_VV
        (s, VFMV_VF (RVV_FASTCOSF_XCOSF_NEG_PI_D_HALF, vlen), u,
         vlen);
      d = VMERGE_VVM_F (u, d, g, vlen);

      g = __riscv_vmflt
        (__riscv_vfabs (r, vlen), RVV_FASTCOSF_XCOSF_TRIG_RANGE_MAX_F,
                        vlen);
      if (__riscv_vcpop (g, vlen) != vlen)
        {
          /* Large range: this is the pinned xcosf rempif reconstruction. */
          VFLOAT rem_hi, rem_lo;
          VINT rem_q;
          RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_cosf_rempif, LMUL)
            (r, &rem_hi, &rem_lo, &rem_q, vlen);

          VINT q2 = __riscv_vand (rem_q, 3, vlen);
          VBOOL rem_positive = __riscv_vmfgt
            (rem_hi, fp_posZero_f, vlen);
          VINT q2_add = VMERGE_VXM_I
            (VMVI_VX (7, vlen), 8, rem_positive, vlen);
          q2 = __riscv_vadd
            (__riscv_vadd (q2, q2, vlen), q2_add, vlen);
          q2 = __riscv_vsra (q2, 1, vlen);

          VBOOL even = __riscv_vmseq
            (__riscv_vand (rem_q, 1, vlen), 0, vlen);
          VFLOAT y = VFMERGE_VFM
            (VFMV_VF (-1.0f, vlen), 0.0f, rem_positive, vlen);
          VFLOAT corr_hi = __riscv_vfsgnjx
            (VFMV_VF (RVV_SLEEF_COSF_REMPI_HALF_PI_HI, vlen), y, vlen);
          VFLOAT corr_lo = __riscv_vfsgnjx
            (VFMV_VF (RVV_SLEEF_COSF_REMPI_HALF_PI_LO, vlen), y, vlen);

          VFLOAT alt_hi, alt_lo;
          RVV_SLEEF_ORIG_DFADD2_VF2_VF2
            (rem_hi, rem_lo, corr_hi, corr_lo, alt_hi, alt_lo, vlen);
          VFLOAT t_hi = VMERGE_VVM_F (rem_hi, alt_hi, even, vlen);
          VFLOAT t_lo = VMERGE_VVM_F (rem_lo, alt_lo, even, vlen);

          /* Preserve source's late non-finite propagation.  SLEEF performs
             the add first and only then ORs the all-ones marker into the
             result; putting the marker into t_hi before the add changes the
             NaN payload/rounding path. */
          VBOOL input_inf = __riscv_vmfeq
            (__riscv_vfabs (r, vlen), fp_posInf_f, vlen);
          VBOOL input_nan = __riscv_vmfne (r, r, vlen);
          VBOOL input_nonfinite = __riscv_vmor
            (input_inf, input_nan, vlen);
          VUINT poison = __riscv_vmerge
            (VMVU_VX (0, vlen), VMVU_VX (~0u, vlen), input_nonfinite,
                      vlen);
          /* Source vsel(g, d, u) selects d for the middle-range lanes and
             the rempif reconstruction for the complementary large lanes. */
          VFLOAT rem = VFADD_VV (t_hi, t_lo, vlen);
          rem = U_AS_F (__riscv_vor (F_AS_U (rem), poison, vlen));
          d = __riscv_vmerge (rem, d, g, vlen);
          q = __riscv_vmerge (q2, q, g, vlen);
        }
    }

  /* xcosf polynomial and final quadrant sign reconstruction. */
  s = VFMUL_VV (d, d, vlen);
  VBOOL flip = __riscv_vmseq
    (__riscv_vand (q, 2, vlen), 0, vlen);
  VUINT sign_mask = __riscv_vmerge
    (VMVU_VX (0, vlen), VMVU_VX (0x80000000u, vlen), flip, vlen);
  d = U_AS_F (__riscv_vxor (F_AS_U (d), sign_mask, vlen));

  u = VFMV_VF (2.6083159809786593541503e-06f, vlen);
  u = VFMADD_VV
    (u, s, VFMV_VF (-0.0001981069071916863322258f, vlen), vlen);
  u = VFMADD_VV
    (u, s, VFMV_VF (0.00833307858556509017944336f, vlen), vlen);
  u = VFMADD_VV
    (u, s, VFMV_VF (-0.166666597127914428710938f, vlen), vlen);

  /* Keep SLEEF xcosf's final operation order and non-fused add:
       vadd(vmul(s, vmul(u, d)), d).
     In particular, do not reassociate this as (s*d)*u+d or turn the
     final add into an explicit FMA.  The source-original object is meant
     to preserve this source-level rounding contract.  */
  VFLOAT ud = VFMUL_VV (u, d, vlen);
  VFLOAT sud = VFMUL_VV (s, ud, vlen);
  VFLOAT result = VFADD_VV (sud, d, vlen);

  return result;
}

#if !RVV_MATH_IS_FASTCOSF_HELPERS_ONLY (V_MATH_FUNC)
static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT d, size_t vlen)
{
  VFLOAT t = d;
  VFLOAT s, u;
  VINT q;

  /* Source line 717: `(float)M_1_PI` is an ordinary caller-FRM scalar
     conversion, followed by the vector multiply-accumulate reduction.  */
  float m1_pi_scalar
    = RVV_SLEEF_ORIG_F64_TO_F32 (RVV_FASTCOSF_SOURCE_M1_PI_D);

  /* SLEEF xfastcosf_u3500 reduction:
       s = d * M_1_PI - 0.5; u = rint(s); q = rint(s).  */
  s = VFMADD_VV (d,
                 VFMV_VF (m1_pi_scalar, vlen),
                 VFMV_VF (-0.5f, vlen), vlen);
  float pi_scalar = RVV_SLEEF_ORIG_F64_TO_F32 (RVV_FASTCOSF_SOURCE_PI_D);
  float pi_neg_scalar = -pi_scalar;
  float pi_half_scalar = pi_scalar * 0.5f;
  u = RVV_LIBMVEC_IMPL_NAME
    (rvv_sleef_orig_fastcosf_vrint, LMUL) (s, vlen);
  q = RVV_VFCVT_X_F_RNE (s, vlen);

  /* Source line 720: reduce with the rounded vector value and the ordinary
     single-precision cast of pi.  Keep the subtraction as a separate vector
     operation, matching vsub_vf_vf_vf(d, (float) M_PI * 0.5f), before the
     fused multiply-add. */
  VFLOAT d_reduced = VFSUB_VF (d, pi_half_scalar, vlen);
  d = VFMADD_VV
    (u, VFMV_VF (pi_neg_scalar, vlen), d_reduced, vlen);

  s = VFMUL_VV (d, d, vlen);

  /* SLEEF xfastcosf_u3500 reduced-precision fast polynomial.  These are
     source literals, not coefficients imported from the clean cosf path. */
  u = VFMV_VF
    (RVV_SLEEF_ORIG_F64_TO_F32 (RVV_FASTCOSF_SOURCE_POLY0_D), vlen);
  u = VFMADD_VV
    (u, s, VFMV_VF
      (RVV_SLEEF_ORIG_F64_TO_F32 (RVV_FASTCOSF_SOURCE_POLY1_D), vlen), vlen);
  u = VFMADD_VV
    (u, s, VFMV_VF
      (RVV_SLEEF_ORIG_F64_TO_F32 (RVV_FASTCOSF_SOURCE_POLY2_D), vlen), vlen);
  u = VFMADD_VV (VFMUL_VV (s, d, vlen), u, d, vlen);

  /* Source flips the sign when q is even.  */
  VBOOL even = __riscv_vmseq (__riscv_vand (q, 1, vlen), 0, vlen);
  VUINT sign_mask = __riscv_vmerge
    (VMVU_VX (0, vlen), VMVU_VX (0x80000000u, vlen), even, vlen);
  u = U_AS_F (__riscv_vxor (F_AS_U (u), sign_mask, vlen));

  /* Source selects xcosf only after the fast kernel has been evaluated.  */
  VBOOL fast_range = __riscv_vmflt
    (__riscv_vfabs (t, vlen), RVV_FASTCOSF_FAST_RANGE_F, vlen);
  if (__riscv_vcpop (fast_range, vlen) != vlen)
    {
      VFLOAT fallback = RVV_LIBMVEC_IMPL_NAME
        (rvv_sleef_orig_fastcosf_xcosf, LMUL) (t, vlen);
      u = VMERGE_VVM_F (fallback, u, fast_range, vlen);
    }

  return u;
}
#endif

#if !RVV_MATH_IS_FASTCOSF_HELPERS_ONLY (V_MATH_FUNC)
#undef LMUL
#undef SIMDLEN
#endif
