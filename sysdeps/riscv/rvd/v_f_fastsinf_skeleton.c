/* Single-precision vector fastsinf function.
   fastsinf(x) computes sine with reduced precision for speed

   Migrated from SLEEF (src/libm/sleefsimdsp.c:xfastsinf_u3500)
   Copyright (C) 2026 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <https://www.gnu.org/licenses/>.  */

#include <vecmath_config.h>
#include <v_math.h>

#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

#ifndef V_F_FASTSINF_SKELETON_ONCE
#define V_F_FASTSINF_SKELETON_ONCE

#define M_1_PI_F 0x1.45f306p-2f  /* 1/π */
#define M_PI_F 0x1.921fb6p+1f     /* π */
/* The source-original keeps SLEEF's <30 raw fast contract.  The public clean
   path uses <1 so its generic MPFR/ULP contract remains within the clean
   accuracy envelope; the clean sinf provider handles the remaining lanes.  */
#define FASTSINF_FAST_LIMIT 1.0f

#endif

/* The clean sinf provider is instantiated privately below.  The public sinf
   skeleton deliberately keeps its implementation TU-local, so an external
   declaration here would leave unresolved __rvv_sinf_* references in every
   fastsinf object.  Keep the fallback in this translation unit while using
   the same clean EH/framework implementation.  */
static VFLOAT
RVV_LIBMVEC_IMPL_NAME (rvv_libmvec_fastsinf_clean_sinf, LMUL)
  (VFLOAT, size_t);

/* SLEEF xfastsinf_u3500 implementation (ULP = 3500) */

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT d = x;
  VFLOAT u, s, t, fallback;
  VINT q;

  SET_ROUNDTONEAREST;

  t = d;

  /* Handle large arguments and NaN/Inf by falling back to the full sinf
     implementation, while masking them out from the fast path to avoid
     spurious exceptions (e.g. float->int conversion overflow).  */
  VUINT vclass = __riscv_vfclass (t, vlen);
  VBOOL NaN_Inf;
  IDENTIFY_F (vclass, class_NaN | class_Inf, NaN_Inf, vlen);

  VFLOAT t_cmp = VFMERGE_VFM (t, fp_posZero_f, NaN_Inf, vlen);
  VBOOL is_fast = __riscv_vmflt (__riscv_vfabs (t_cmp, vlen),
                                 VFMV_VF_F (FASTSINF_FAST_LIMIT, vlen),
                                 vlen);
  VUINT expo_t = __riscv_vsrl (F_AS_U_F (__riscv_vfabs (t, vlen)),
                               MAN_LEN_F, vlen);
  expo_t = __riscv_vand (expo_t, 0xff, vlen);
  VBOOL tiny_x = __riscv_vmsleu (expo_t,
                                 EXP_BIAS_F - MAN_LEN_F - 5, vlen);
  VBOOL need_fallback = __riscv_vmor (NaN_Inf, __riscv_vmnot (is_fast, vlen),
                                       vlen);
  need_fallback = __riscv_vmor (need_fallback, tiny_x, vlen);
  UINT nb_need_fallback = __riscv_vcpop (need_fallback, vlen);
  if (nb_need_fallback > 0)
    {
      /* Prepare the fallback/special provider before entering the finite
         polynomial.  Fast lanes use zero only as a harmless provider input;
         fallback lanes retain their original values.  */
      fallback = RVV_LIBMVEC_IMPL_NAME (rvv_libmvec_fastsinf_clean_sinf, LMUL) (
          VFMERGE_VFM (t, fp_posZero_f, __riscv_vmnot (need_fallback, vlen),
                       vlen),
          vlen);
      d = VFMERGE_VFM (d, fp_posZero_f, need_fallback, vlen);
    }

  /* Argument reduction: s = d * (1/π) */
  s = __riscv_vfmul (d, VFMV_VF (M_1_PI_F, vlen), vlen);

  /* Round to nearest integer */
  u = __riscv_vfcvt_f (__riscv_vfcvt_x (s, vlen), vlen);
  q = __riscv_vfcvt_x (s, vlen);

  /* Reduced argument: d = u * (-π) + d */
  d = VFMADD_VV (u, VFMV_VF (-M_PI_F, vlen), d, vlen);

  /* Compute s = d^2 */
  s = __riscv_vfmul (d, d, vlen);

  /* Polynomial approximation with SLEEF coefficients */
  u = VFMV_VF (-0.1881748176e-3f, vlen);
  u = VFMADD_VV (u, s, VFMV_VF (+0.8323502727e-2f, vlen), vlen);
  u = VFMADD_VV (u, s, VFMV_VF (-0.1666651368e+0f, vlen), vlen);

  /* Final polynomial: u = (s * d) * u + d */
  VFLOAT s_times_d = __riscv_vfmul (s, d, vlen);
  u = VFMADD_VV (s_times_d, u, d, vlen);

  /* Sign adjustment based on quadrant (q & 1) */
  VBOOL sign_flip = __riscv_vmseq (VAND_VV_I (q, VMVI_VX (1, vlen), vlen),
                                    VMVI_VX (1, vlen), vlen);
  VUINT sign_mask = __riscv_vmerge (VMVU_VX (0, vlen), 0x80000000, sign_flip, vlen);
  u = U_AS_F (VXOR_VV (F_AS_U (u), sign_mask, vlen));

  if (nb_need_fallback > 0)
    u = VMERGE_VVM_F (u, fallback, need_fallback, vlen);

  VBOOL negzero = __riscv_vmseq (F_AS_U_F (t), VMVU_VX_F (0x80000000, vlen),
                                 vlen);
  u = VFMERGE_VFM (u, fp_negZero_f, negzero, vlen);

  RESTORE_FRM;
  return u;
}

/* Instantiate the clean sinf fallback after the fast entry point so the
   provider remains private to this translation unit.  The source fastsinf
   entry uses the corresponding source-original u35 fallback; clean uses the
   public clean sinf EH path here by design.  */
#pragma push_macro("LMUL")
#pragma push_macro("SIMDLEN")
#pragma push_macro("V_MATH_FUNC")
#undef V_MATH_FUNC
#define V_MATH_FUNC rvv_libmvec_fastsinf_clean_sinf
#include "v_f_sinf_skeleton.c"
#pragma pop_macro("V_MATH_FUNC")
#pragma pop_macro("SIMDLEN")
#pragma pop_macro("LMUL")

#undef LMUL
#undef SIMDLEN
