/* Single-precision vector log1pf function.
   log1pf(x) computes log(1 + x) accurately for small x

   Based on the pinned SLEEF xlog1pf implementation.

   Copyright (C) 2026 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <https://www.gnu.org/licenses/>.  */

#include <float.h>
#include <vecmath_config.h>
#include <v_math.h>

/* The finite path, DF operation order, polynomial, and large-input fallback
   are shared with the source-original comparison entry.  The importer-name
   predicate in v_math.h makes this include contribute helpers only; the
   public clean entry below owns the VECLIBM-style EH envelope.  */
#include "v_f_log1p_original_skeleton.c"

#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT d = x;
  VFLOAT r_special;
  VBOOL special_args;

  SET_ROUNDTONEAREST;

  /* Classify before doing arithmetic on caller NaNs, infinities, domain
     endpoints, or signed zero.  This is the clean VECLIBM EH envelope; the
     source-shaped finite kernel is called only after sanitization.  */
  VUINT vclass = __riscv_vfclass (d, vlen);
  VBOOL x_nan, x_snan, x_pos_inf, x_zero;
  IDENTIFY_F (vclass, class_NaN, x_nan, vlen);
  IDENTIFY_F (vclass, class_sNaN, x_snan, vlen);
  IDENTIFY_F (vclass, class_posInf, x_pos_inf, vlen);
  IDENTIFY_F (vclass, class_Zero, x_zero, vlen);

  VBOOL cmp_mask = __riscv_vmor (x_nan, x_pos_inf, vlen);
  VFLOAT d_cmp = __riscv_vfmerge (d, fp_posZero_f, cmp_mask, vlen);
  VBOOL x_lt_neg1 = __riscv_vmflt (d_cmp, fp_negOne_f, vlen);
  VBOOL x_eq_neg1 = __riscv_vmfeq (d_cmp, fp_negOne_f, vlen);

  special_args = __riscv_vmor (x_nan, x_pos_inf, vlen);
  special_args = __riscv_vmor (special_args, x_lt_neg1, vlen);
  special_args = __riscv_vmor (special_args, x_eq_neg1, vlen);
  special_args = __riscv_vmor (special_args, x_zero, vlen);

  r_special = VFMV_VF_F (0.0f, vlen);
  r_special = __riscv_vfmerge
    (r_special, fp_qNaN_f, __riscv_vmor (x_nan, x_lt_neg1, vlen), vlen);
  r_special = __riscv_vfmerge (r_special, fp_negInf_f, x_eq_neg1, vlen);
  r_special = __riscv_vfmerge (r_special, fp_posInf_f, x_pos_inf, vlen);
  r_special = __riscv_vmerge (r_special, d, x_zero, vlen);

  /* Match the VECLIBM exception side effects while keeping the raw lanes out
     of the finite kernel.  log1p(-1) is a pole (divide by zero), and x < -1
     plus signaling NaN inputs require invalid.  The masked reciprocal and
     masked sNaN addition are used only as EH triggers; their prepared values
     are merged back into r_special before the input quarantine.  */
  VFLOAT pole_zero = __riscv_vfsgnj
    (VFMV_VF_F (fp_posZero_f, vlen), d, vlen);
  VFLOAT pole_effect = __riscv_vfrec7 (x_eq_neg1, pole_zero, vlen);
  r_special = __riscv_vmerge (r_special, pole_effect, x_eq_neg1, vlen);

  VBOOL invalid_args = __riscv_vmor (x_lt_neg1, x_snan, vlen);
  VFLOAT invalid_seed = __riscv_vfmerge
    (VFMV_VF_F (fp_posZero_f, vlen), fp_sNaN_f, invalid_args, vlen);
  VFLOAT invalid_effect = __riscv_vfadd
    (invalid_args, invalid_seed, invalid_seed, vlen);
  r_special = __riscv_vmerge (r_special, invalid_effect, invalid_args, vlen);

  VBOOL finite_args = __riscv_vmnot (special_args, vlen);
  UINT nb_finite_args = __riscv_vcpop (finite_args, vlen);
  if (nb_finite_args == 0)
    {
      RESTORE_FRM;
      return r_special;
    }

  /* Prevent raw special values from reaching the source finite path.  The
     final merge below restores the scalar-compatible special result and
     preserves the source input bits for both signed-zero lanes.  */
  d = __riscv_vfmerge (d, fp_posZero_f, special_args, vlen);

  VFLOAT r;
  if (nb_finite_args == vlen)
    r = RVV_SLEEF_ORIG_MANGLE
      (rvv_sleef_orig_log1pf_kernel_m, LMUL) (d, vlen);
  else
    r = RVV_SLEEF_ORIG_MANGLE
      (rvv_sleef_orig_log1pf_kernel_masked_m, LMUL)
      (d, finite_args, vlen);
  r = __riscv_vmerge (r, r_special, special_args, vlen);

  /* The source-shaped DF reduction can collapse -FLT_TRUE_MIN to +0 on RVV.
     The public clean routine must retain the correctly rounded finite endpoint
     value.  Keep this narrow target adapter outside the source-original entry;
     it does not change the source reduction, polynomial, DF order, or
     fallback identity.  */
  VFLOAT x_abs = __riscv_vfabs (x, vlen);
  x_abs = __riscv_vfmerge (x_abs, fp_posInf_f, special_args, vlen);
  VBOOL tiny_finite = __riscv_vmflt (x_abs, 0x1.0p-24f, vlen);
  r = VMERGE_VVM_F (r, x, __riscv_vmandn (tiny_finite, special_args, vlen),
                    vlen);

  RESTORE_FRM;
  return r;
}

#undef LMUL
#undef SIMDLEN
