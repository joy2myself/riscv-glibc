/* Single-precision vector sinpif function.
   sinpif(x) computes sin(pi * x)

   The normal path follows SLEEF xsinpif_u05/sinpifk trig-pi reduction and
   polynomial reconstruction.  The special-case handling is adapted for glibc
   libmvec.

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

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT x, size_t vlen)
{
  VFLOAT vx_orig = x;
  VFLOAT vx, ax, vy, vy_special;
  VBOOL special_args;
  VUINT expo_x;

  SET_ROUNDTONEAREST;

  vx = vx_orig;

  /* Special cases: NaN/Inf, tiny |x| (sinpif rounds to pi*x), and the SLEEF
     large-finite zeroing region |x| > 1e7f.  */
  VUINT vclass = __riscv_vfclass (vx, vlen);
  VBOOL NaN_Inf;
  IDENTIFY_F (vclass, class_NaN | class_Inf, NaN_Inf, vlen);

  ax = __riscv_vfabs (vx, vlen);
  expo_x = __riscv_vsrl (F_AS_U_F (ax), MAN_LEN_F, vlen);
  expo_x = __riscv_vand (expo_x, 0xff, vlen);
  VBOOL tiny_x = __riscv_vmsleu
    (expo_x, RVV_LIBMVEC_TRIGPI_TINY_EXPONENT, vlen);
  VBOOL x_large = __riscv_vmfgt
    (ax, RVV_SLEEF_TRIGPI_RANGE_MAX_F, vlen);
  special_args = __riscv_vmor (NaN_Inf, __riscv_vmor (tiny_x, x_large, vlen),
                               vlen);

  UINT nb_special_args = __riscv_vcpop (special_args, vlen);
  if (nb_special_args > 0)
    {
      /* Default special-case result for large |x|: signed zero.  */
      vy_special = __riscv_vfsgnjx (VFMV_VF_F (fp_posZero_f, vlen), vx_orig,
                                    vlen);

      /* For tiny |x|, sin(pi*x) ~= pi*x.  */
      VFLOAT tiny_input
          = VFMERGE_VFM (vx_orig, fp_posZero_f, __riscv_vmnot (tiny_x, vlen),
                         vlen);
      VFLOAT sin_tiny = __riscv_vfmul
        (tiny_input, RVV_PI_F, vlen);
      sin_tiny = __riscv_vfmacc (sin_tiny, tiny_input,
                                 VFMV_VF_F (RVV_PI_LO_F, vlen),
                                 vlen);
      vy_special = VMERGE_VVM_F (vy_special, sin_tiny, tiny_x, vlen);
      VBOOL x_zero = __riscv_vmfeq (vx_orig, fp_posZero_f, vlen);
      vy_special = VMERGE_VVM_F (vy_special, vx_orig, x_zero, vlen);

      /* Substitute Inf with sNaN and propagate NaNs via x+x.  */
      VBOOL inf_mask;
      IDENTIFY_F (vclass, class_Inf, inf_mask, vlen);
      VFLOAT vy_NaN_Inf = __riscv_vfmerge (vx_orig, fp_sNaN_f, inf_mask,
                                            vlen);
      vy_NaN_Inf = __riscv_vfadd (NaN_Inf, vy_NaN_Inf, vy_NaN_Inf, vlen);
      vy_special = VMERGE_VVM_F (vy_special, vy_NaN_Inf, NaN_Inf, vlen);

      /* Mask special lanes out from the main path to avoid spurious
         exceptions (e.g. NaN->int conversion).  */
      vx = VFMERGE_VFM (vx, fp_posZero_f, special_args, vlen);
      ax = VFMERGE_VFM (ax, fp_posZero_f, special_args, vlen);
    }

  /* SLEEF sinpifk reduction: u = 4*x, q = trunc(u), followed by the
     source signed correction and even-quotient reduction.  Special lanes
     were sanitized above; keeping this finite graph source-shaped avoids
     turning the clean EH adaptation into a second algorithm.  */
  VFLOAT u4 = __riscv_vfmul (vx, RVV_SLEEF_TRIG_REMPI_FOUR, vlen);
  VINT n = __riscv_vfcvt_rtz_x (u4, vlen);
  VUINT n_bits = I_AS_U_F (n);
  VUINT n_sign = __riscv_vsrl (n_bits, BIT_WIDTH_F - 1, vlen);
  VUINT n_correction = __riscv_vxor
    (n_sign, VMVU_VX_F (1, vlen), vlen);
  n = __riscv_vadd (n, U_AS_I_F (n_correction), vlen);
  n = __riscv_vand (n, VMVI_VX_F (~1, vlen), vlen);
  VFLOAT n_flt = __riscv_vfcvt_f (n, vlen);
  VFLOAT t_arg = __riscv_vfsub (u4, n_flt, vlen);
  VFLOAT s = __riscv_vfmul (t_arg, t_arg, vlen);
  VFLOAT s2_hi, s2_lo;
  DFMUL_VF2_VF_VF (t_arg, t_arg, s2_hi, s2_lo, vlen);

  VBOOL sin_use_c
      = __riscv_vmseq (__riscv_vand (I_AS_U_F (n), 2, vlen), 2, vlen);

  VFLOAT p = VFMV_VF_F (RVV_SLEEF_TRIGPI_POLY_COS0_F, vlen);
  p = VFMERGE_VFM (p, RVV_SLEEF_TRIGPI_POLY_SIN0_F, sin_use_c, vlen);
  VFLOAT c = VFMV_VF_F (RVV_SLEEF_TRIGPI_POLY_COS1_F, vlen);
  c = VFMERGE_VFM (c, RVV_SLEEF_TRIGPI_POLY_SIN1_F, sin_use_c, vlen);
  p = __riscv_vfmadd (p, s, c, vlen);
  c = VFMV_VF_F (RVV_SLEEF_TRIGPI_POLY_COS2_F, vlen);
  c = VFMERGE_VFM (c, RVV_SLEEF_TRIGPI_POLY_SIN2_F, sin_use_c, vlen);
  p = __riscv_vfmadd (p, s, c, vlen);

  VFLOAT x_hi, x_lo;
  VFLOAT c_hi = VFMV_VF_F (RVV_SLEEF_TRIGPI_DF_FIRST_COS_HI_F, vlen);
  VFLOAT c_lo = VFMV_VF_F (RVV_SLEEF_TRIGPI_DF_FIRST_COS_LO_F, vlen);
  c_hi = VFMERGE_VFM
    (c_hi, RVV_SLEEF_TRIGPI_DF_FIRST_SIN_HI_F, sin_use_c, vlen);
  c_lo = VFMERGE_VFM
    (c_lo, RVV_SLEEF_TRIGPI_DF_FIRST_SIN_LO_F, sin_use_c, vlen);
  DFADD2_VF2_VF_VF2 (__riscv_vfmul (p, s, vlen), c_hi, c_lo, x_hi, x_lo,
                     vlen);

  VFLOAT tmp_hi, tmp_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2 (s2_hi, s2_lo, x_hi, x_lo, tmp_hi, tmp_lo, vlen);
  c_hi = VFMV_VF_F (RVV_SLEEF_TRIGPI_DF_SECOND_COS_HI_F, vlen);
  c_lo = VFMV_VF_F (RVV_SLEEF_TRIGPI_DF_SECOND_COS_LO_F, vlen);
  c_hi = VFMERGE_VFM
    (c_hi, RVV_SLEEF_TRIGPI_DF_SECOND_SIN_HI_F, sin_use_c, vlen);
  c_lo = VFMERGE_VFM
    (c_lo, RVV_SLEEF_TRIGPI_DF_SECOND_SIN_LO_F, sin_use_c, vlen);
  DFADD2_VF2_VF2_VF2 (tmp_hi, tmp_lo, c_hi, c_lo, x_hi, x_lo, vlen);

  VFLOAT mul_hi = t_arg;
  VFLOAT mul_lo = VFMV_VF_F (fp_posZero_f, vlen);
  mul_hi = VMERGE_VVM_F (mul_hi, s2_hi, sin_use_c, vlen);
  mul_lo = VMERGE_VVM_F (mul_lo, s2_lo, sin_use_c, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2 (x_hi, x_lo, mul_hi, mul_lo, tmp_hi, tmp_lo, vlen);
  x_hi = tmp_hi;
  x_lo = tmp_lo;

  VFLOAT add_hi, add_lo;
  DFADD2_VF2_VF (x_hi, x_lo, VFMV_VF_F (fp_posOne_f, vlen), add_hi, add_lo,
                 vlen);
  x_hi = VMERGE_VVM_F (x_hi, add_hi, sin_use_c, vlen);
  x_lo = VMERGE_VVM_F (x_lo, add_lo, sin_use_c, vlen);
  vy = __riscv_vfadd (x_hi, x_lo, vlen);

  /* SLEEF sinpifk flips both double-float components before collapse.  */
  VUINT sign_bits
      = __riscv_vsll (__riscv_vand (I_AS_U_F (n), 4, vlen), BIT_WIDTH_F - 3,
                      vlen);
  x_hi = U_AS_F_F (__riscv_vxor (F_AS_U_F (x_hi), sign_bits, vlen));
  x_lo = U_AS_F_F (__riscv_vxor (F_AS_U_F (x_lo), sign_bits, vlen));
  vy = __riscv_vfadd (x_hi, x_lo, vlen);

  /* Preserve clean's exact odd signed-zero endpoint contract.  */
  VBOOL sin_zero = __riscv_vmand (
      __riscv_vmseq (F_AS_U_F (t_arg), 0, vlen),
      __riscv_vmnot (sin_use_c, vlen), vlen);
  VFLOAT signed_zero = __riscv_vfsgnj
    (VFMV_VF_F (fp_posZero_f, vlen), vx_orig, vlen);
  vy = __riscv_vmerge (vy, signed_zero, sin_zero, vlen);

  if (nb_special_args > 0)
    vy = VMERGE_VVM_F (vy, vy_special, special_args, vlen);

  RESTORE_FRM;
  return vy;
}

#undef LMUL
#undef SIMDLEN
