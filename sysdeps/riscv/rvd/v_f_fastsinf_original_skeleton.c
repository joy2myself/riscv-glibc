/* Single-precision source-original SLEEF fastsinf skeleton.

   Target-owned translation of SLEEF xfastsinf_u3500.  It keeps the source
   1/pi reduction, reduced-precision polynomial, quadrant sign operation and
   late xsinf fallback.  The fallback is instantiated from the target-owned
   source-original sinf skeleton in this translation unit, not from clean
   sinf and not from a SLEEF object or source file.

   Clean-only quarantine, tiny-domain replacement and fallback reordering are
   intentionally not imported into this comparison implementation.  The
   object is built with SLEEF's -ffp-contract=off policy.  */

#include <vecmath_config.h>
#include <v_math.h>

#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

/* xfastsinf_u3500 calls xsinf after the fast kernel.  Keep that source
   fallback target-owned and private to this translation unit.  The ordinary
   sinf source-original skeleton below supplies rempif and its target-owned
   helpers, but xsinf itself must be translated here: it is the u35 entry,
   whereas v_f_sinf_original_skeleton.c is the separate xsinf_u1/u10
   comparator.  */
#define RVV_SLEEF_ORIG_KEEP_LMUL
#undef V_MATH_FUNC
#define V_MATH_FUNC rvv_sleef_orig_fastsinf_sinf
#include "v_f_sinf_original_skeleton.c"
#undef V_MATH_FUNC
#define V_MATH_FUNC fastsinf_original
#undef RVV_SLEEF_ORIG_KEEP_LMUL

#define RVV_SLEEF_ORIG_FASTSINF_REMPIF \
  RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_sinf_rempif, LMUL)
#define RVV_SLEEF_ORIG_FASTSINF_VFCVT_X_RNE(value, vlen) \
  RVV_VFCVT_X_F_RNE ((value), (vlen))

#ifndef RVV_SLEEF_ORIG_FASTSINF_U35_CONSTANTS
# define RVV_SLEEF_ORIG_FASTSINF_U35_CONSTANTS
# define RVV_SLEEF_ORIG_FASTSINF_M1_PI 0.31830988618379067154
# define RVV_SLEEF_ORIG_FASTSINF_PI_A2 3.1414794921875f
# define RVV_SLEEF_ORIG_FASTSINF_PI_B2 0.00011315941810607910156f
# define RVV_SLEEF_ORIG_FASTSINF_PI_C2 1.9841872589410058936e-09f
# define RVV_SLEEF_ORIG_FASTSINF_PI_A 3.140625f
# define RVV_SLEEF_ORIG_FASTSINF_PI_B 0.0009670257568359375f
# define RVV_SLEEF_ORIG_FASTSINF_PI_C 6.2771141529083251953e-07f
# define RVV_SLEEF_ORIG_FASTSINF_PI_D 1.2154201256553420762e-10f
# define RVV_SLEEF_ORIG_FASTSINF_TRIGRANGEMAX2 125.0f
# define RVV_SLEEF_ORIG_FASTSINF_TRIGRANGEMAX 39000.0f
#endif

/* Target-owned translation of SLEEF xsinf (Sleef_sinf_u35), exactly the
   fallback selected by xfastsinf_u3500.  It deliberately remains separate
   from the xsinf_u1 translation used by sinf_original.  */
static VFLOAT
RVV_LIBMVEC_IMPL_NAME (rvv_sleef_orig_fastsinf_u35, LMUL)
  (VFLOAT d, size_t vlen)
{
  VFLOAT r = d;
  VFLOAT u, s;
  VINT q;

  u = VFMUL_VF (d, RVV_SLEEF_ORIG_FASTSINF_M1_PI, vlen);
  q = RVV_SLEEF_ORIG_FASTSINF_VFCVT_X_RNE (u, vlen);
  u = __riscv_vfcvt_f (q, vlen);

  d = VFMADD_VV
    (u, VFMV_VF (-RVV_SLEEF_ORIG_FASTSINF_PI_A2, vlen), d, vlen);
  d = VFMADD_VV
    (u, VFMV_VF (-RVV_SLEEF_ORIG_FASTSINF_PI_B2, vlen), d, vlen);
  d = VFMADD_VV
    (u, VFMV_VF (-RVV_SLEEF_ORIG_FASTSINF_PI_C2, vlen), d, vlen);

  VBOOL small2 = __riscv_vmflt
    (__riscv_vfabs (r, vlen), RVV_SLEEF_ORIG_FASTSINF_TRIGRANGEMAX2,
                    vlen);
  if (__riscv_vcpop (small2, vlen) != vlen)
    {
      s = u;
      u = VFMADD_VV
        (s, VFMV_VF (-RVV_SLEEF_ORIG_FASTSINF_PI_A, vlen), r, vlen);
      u = VFMADD_VV
        (s, VFMV_VF (-RVV_SLEEF_ORIG_FASTSINF_PI_B, vlen), u, vlen);
      u = VFMADD_VV
        (s, VFMV_VF (-RVV_SLEEF_ORIG_FASTSINF_PI_C, vlen), u, vlen);
      u = VFMADD_VV
        (s, VFMV_VF (-RVV_SLEEF_ORIG_FASTSINF_PI_D, vlen), u, vlen);
      d = VMERGE_VVM_F (u, d, small2, vlen);

      VBOOL small = __riscv_vmflt
        (__riscv_vfabs (r, vlen), RVV_SLEEF_ORIG_FASTSINF_TRIGRANGEMAX,
                        vlen);
      if (__riscv_vcpop (small, vlen) != vlen)
        {
          VFLOAT rem_hi, rem_lo;
          VINT rem_q;
          RVV_SLEEF_ORIG_FASTSINF_REMPIF
            (r, &rem_hi, &rem_lo, &rem_q, vlen);

          VINT q2 = __riscv_vand (rem_q, 3, vlen);
          VBOOL rem_positive = __riscv_vmfgt
            (rem_hi, fp_posZero_f, vlen);
          VINT q2_add = VMERGE_VXM_I
            (VMVI_VX (1, vlen), 2, rem_positive, vlen);
          q2 = __riscv_vadd (__riscv_vadd (q2, q2, vlen), q2_add, vlen);
          q2 = __riscv_vsra (q2, 2, vlen);

          VBOOL rem_q_odd = __riscv_vmseq
            (__riscv_vand (rem_q, 1, vlen), 1, vlen);
          VFLOAT corr_hi = __riscv_vfsgnjx
            (VFMV_VF (-3.1415927410125732422f * 0.5f, vlen), rem_hi,
             vlen);
          VFLOAT corr_lo = __riscv_vfsgnjx
            (VFMV_VF (-8.7422776573475857731e-08f * -0.5f, vlen),
             rem_hi, vlen);
          VFLOAT alt_hi, alt_lo;
          RVV_SLEEF_ORIG_DFADD2_VF2_VF2
            (rem_hi, rem_lo, corr_hi, corr_lo, alt_hi, alt_lo, vlen);
          VFLOAT selected_hi = VMERGE_VVM_F
            (rem_hi, alt_hi, rem_q_odd, vlen);
          VFLOAT selected_lo = VMERGE_VVM_F
            (rem_lo, alt_lo, rem_q_odd, vlen);

          /* This is SLEEF's source-native late Inf/NaN bit operation.  */
          VBOOL input_inf = __riscv_vmfeq
            (__riscv_vfabs (r, vlen), fp_posInf_f, vlen);
          VBOOL input_nan = __riscv_vmfne (r, r, vlen);
          VBOOL input_nonfinite = __riscv_vmor
            (input_inf, input_nan, vlen);
          VUINT poison = __riscv_vmerge
            (VMVU_VX (0, vlen), VMVU_VX (~0u, vlen), input_nonfinite,
             vlen);
          VFLOAT selected = VFADD_VV (selected_hi, selected_lo, vlen);
          selected = U_AS_F
            (__riscv_vor (F_AS_U (selected), poison, vlen));

          VBOOL very_large = __riscv_vmnot (small, vlen);
          q = __riscv_vmerge (q, q2, very_large, vlen);
          d = __riscv_vmerge (d, selected, very_large, vlen);
        }
    }

  s = VFMUL_VV (d, d, vlen);
  VBOOL q_odd = __riscv_vmseq
    (__riscv_vand (q, 1, vlen), 1, vlen);
  VFLOAT sign_flip = VFMERGE_VFM
    (VFMV_VF (0.0f, vlen), -0.0f, q_odd, vlen);
  d = __riscv_vfsgnjx (d, sign_flip, vlen);

  u = VFMV_VF (2.6083159809786593541503e-06f, vlen);
  u = VFMADD_VV
    (u, s, VFMV_VF (-0.0001981069071916863322258f, vlen), vlen);
  u = VFMADD_VV
    (u, s, VFMV_VF (0.00833307858556509017944336f, vlen), vlen);
  u = VFMADD_VV
    (u, s, VFMV_VF (-0.166666597127914428710938f, vlen), vlen);
  u = VFADD_VV (VFMUL_VV (s, VFMUL_VV (u, d, vlen), vlen), d, vlen);

  /* Preserve xsinf's source-native negative-zero endpoint.  */
  VBOOL r_negzero = __riscv_vmseq
    (U_AS_I (F_AS_U (r)), 0x80000000, vlen);
  return VMERGE_VVM_F (u, r, r_negzero, vlen);
}

#define RVV_SLEEF_ORIG_FASTSINF_FALLBACK \
  RVV_LIBMVEC_IMPL_NAME (rvv_sleef_orig_fastsinf_u35, LMUL)

#ifndef RVV_SLEEF_ORIG_FASTSINF_CONSTANTS
# define RVV_SLEEF_ORIG_FASTSINF_CONSTANTS
# define RVV_SLEEF_ORIG_FASTSINF_PI 3.14159265358979323846
#endif

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT d, size_t vlen)
{
  VFLOAT t = d;
  VFLOAT u, s;
  VINT q;

  /* Source xfastsinf_u3500: u = rint(d/pi), q = rint(d/pi).  These are
     explicit RNE conversions; surrounding floating-point arithmetic keeps
     the caller's FRM.  */
  s = VFMUL_VF (d, RVV_SLEEF_ORIG_FASTSINF_M1_PI, vlen);
  u = __riscv_vfcvt_f
    (RVV_SLEEF_ORIG_FASTSINF_VFCVT_X_RNE (s, vlen), vlen);
  q = RVV_SLEEF_ORIG_FASTSINF_VFCVT_X_RNE (s, vlen);

  /* Source spells this as -(float)M_PI: round the positive double to float
     under the caller's FRM, then apply the sign bit.  */
  VFLOAT pi = VFMV_VF (RVV_SLEEF_ORIG_FASTSINF_PI, vlen);
  pi = VFNEG_V (pi, vlen);
  d = VFMADD_VV (u, pi, d, vlen);
  s = VFMUL_VV (d, d, vlen);

  /* SLEEF reduced-precision fast polynomial.  */
  u = VFMV_VF (-0.1881748176e-3, vlen);
  u = VFMADD_VV (u, s, VFMV_VF (+0.8323502727e-2, vlen), vlen);
  u = VFMADD_VV (u, s, VFMV_VF (-0.1666651368e+0, vlen), vlen);
  u = VFMADD_VV (VFMUL_VV (s, d, vlen), u, d, vlen);

  /* Source flips the sign when q is odd.  */
  VBOOL odd = __riscv_vmseq (__riscv_vand (q, 1, vlen), 1, vlen);
  VUINT sign_mask = __riscv_vmerge
    (VMVU_VX (0, vlen), VMVU_VX (0x80000000u, vlen), odd, vlen);
  u = U_AS_F (__riscv_vxor (F_AS_U (u), sign_mask, vlen));

  /* Source calls xsinf only after the fast kernel and selects late.  */
  VBOOL fast_range = __riscv_vmflt
    (__riscv_vfabs (t, vlen), 30.0f, vlen);
  if (__riscv_vcpop (fast_range, vlen) != vlen)
    {
      VFLOAT fallback = RVV_SLEEF_ORIG_FASTSINF_FALLBACK (t, vlen);
      u = VMERGE_VVM_F (fallback, u, fast_range, vlen);
    }

  return u;
}

#undef RVV_SLEEF_ORIG_FASTSINF_FALLBACK
#undef RVV_SLEEF_ORIG_FASTSINF_VFCVT_X_RNE
#undef RVV_SLEEF_ORIG_FASTSINF_REMPIF
#undef RVV_SLEEF_ORIG_FASTSINF_U35_CONSTANTS
#undef RVV_SLEEF_ORIG_FASTSINF_PI
#undef RVV_SLEEF_ORIG_FASTSINF_M1_PI
#undef RVV_SLEEF_ORIG_FASTSINF_CONSTANTS
#undef LMUL
#undef SIMDLEN
