/* Single-precision vector sincospif function.
   sincospif(x) computes both sin(pi * x) and cos(pi * x).

   The normal path follows SLEEF xsincospif_u05 trig-pi reduction and
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

/* This version computes both sin(pi*x) and cos(pi*x) simultaneously.  */

static void
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT_F vx, VFLOAT_F *vsin_out,
                                           VFLOAT_F *vcos_out, size_t vlen)
{
  VFLOAT_F vx_orig, vx_core, ax;
  VFLOAT_F vsin, vcos, vsin_special, vcos_special;
  VBOOL_F special_args;
  VUINT_F expo_x;

  SET_ROUNDTONEAREST;

  vx_orig = vx;
  vx_core = vx_orig;

  /* Special cases: NaN/Inf, tiny |x| (sincospif rounds to (pi*x, 1)), and
     |x| >= 2^24 (float can represent only even integers, so sincospif(x) is
     (+/-0, 1)).  */
  VUINT_F vclass = __riscv_vfclass (vx_orig, vlen);
  VBOOL_F NaN_Inf;
  IDENTIFY_F (vclass, class_NaN | class_Inf, NaN_Inf, vlen);

  ax = __riscv_vfabs (vx_orig, vlen);
  expo_x = __riscv_vsrl (F_AS_U_F (ax), MAN_LEN_F, vlen);
  expo_x = __riscv_vand (expo_x, 0xff, vlen);
  VBOOL_F tiny_x = __riscv_vmsleu (expo_x, EXP_BIAS_F - MAN_LEN_F - 5, vlen);
  VBOOL_F x_large = __riscv_vmsgeu (expo_x, EXP_BIAS_F + 24, vlen);
  special_args = __riscv_vmor (NaN_Inf, __riscv_vmor (tiny_x, x_large, vlen),
                               vlen);

  UINT nb_special_args = __riscv_vcpop (special_args, vlen);
  if (nb_special_args > 0)
    {
      /* Default special-case results for large |x|: (+/-0, 1).  */
      vsin_special = __riscv_vfsgnj (VFMV_VF_F (fp_posZero_f, vlen), vx_orig,
                                     vlen);
      vcos_special = VFMV_VF_F (fp_posOne_f, vlen);

      /* For tiny |x|, sin(pi*x) ~= pi*x and cos(pi*x) rounds to 1. Use 1 + x
         so nonzero x raises inexact.  */
      VFLOAT_F tiny_input
          = VFMERGE_VFM (vx_orig, fp_posZero_f, __riscv_vmnot (tiny_x, vlen),
                         vlen);
      VFLOAT_F sin_tiny = __riscv_vfmul
        (tiny_input, RVV_PI_F, vlen);
      sin_tiny = __riscv_vfmacc
        (sin_tiny, tiny_input, VFMV_VF_F (RVV_PI_LO_F, vlen),
         vlen);
      vsin_special = VMERGE_VVM_F (vsin_special, sin_tiny, tiny_x, vlen);
      VBOOL_F x_zero = __riscv_vmfeq (vx_orig, fp_posZero_f, vlen);
      vsin_special = VMERGE_VVM_F (vsin_special, vx_orig, x_zero, vlen);

      VFLOAT_F cos_tiny = __riscv_vfadd (tiny_x, vx_orig, fp_posOne_f, vlen);
      vcos_special = VMERGE_VVM_F (vcos_special, cos_tiny, tiny_x, vlen);

      /* Substitute Inf with sNaN and propagate NaNs via x+x.  */
      VBOOL_F inf_mask;
      IDENTIFY_F (vclass, class_Inf, inf_mask, vlen);
      VFLOAT_F v_NaN_Inf = __riscv_vfmerge (vx_orig, fp_sNaN_f, inf_mask, vlen);
      v_NaN_Inf = __riscv_vfadd (NaN_Inf, v_NaN_Inf, v_NaN_Inf, vlen);
      vsin_special = VMERGE_VVM_F (vsin_special, v_NaN_Inf, NaN_Inf, vlen);
      vcos_special = VMERGE_VVM_F (vcos_special, v_NaN_Inf, NaN_Inf, vlen);

      /* Mask special lanes out from the main path to avoid spurious
         exceptions (e.g. NaN->int conversion).  */
      vx_core = VFMERGE_VFM (vx_core, fp_posZero_f, special_args, vlen);
    }

  /* SLEEF xsincospif_u05 signed reduction: q = trunc(4*x), followed by
     q += ((q >> 31) ^ 1) and q &= ~1.  Special lanes were sanitized above;
     keep this finite graph source-shaped and retain clean endpoint repair
     below as a separate framework contract.  */
  VFLOAT_F u4 = __riscv_vfmul (vx_core, 4.0f, vlen);
  VINT_F q = __riscv_vfcvt_rtz_x (u4, vlen);
  VUINT_F q_bits = I_AS_U_F (q);
  VUINT_F q_sign = __riscv_vsrl (q_bits, 31, vlen);
  VUINT_F q_correction = __riscv_vxor
    (q_sign, VMVU_VX_F (1, vlen), vlen);
  q = __riscv_vadd (q, U_AS_I_F (q_correction), vlen);
  q = __riscv_vand (q, VMVI_VX_F (~1, vlen), vlen);
  VFLOAT_F q_flt = __riscv_vfcvt_f (q, vlen);
  VFLOAT_F t_arg = __riscv_vfsub (u4, q_flt, vlen);
  VFLOAT_F s = __riscv_vfmul (t_arg, t_arg, vlen);
  VFLOAT_F s2_hi, s2_lo;
  RVV_SLEEF_ORIG_DFMUL_VF_VF (t_arg, t_arg, s2_hi, s2_lo, vlen);

  VFLOAT_F p = VFMV_VF_F (0.3093842054e-6f, vlen);
  p = __riscv_vfmadd (p, s, VFMV_VF_F (-0.3657307388e-4f, vlen), vlen);
  p = __riscv_vfmadd (p, s, VFMV_VF_F (0.2490393585e-2f, vlen), vlen);
  VFLOAT_F x_hi, x_lo, tmp_hi, tmp_lo;
  RVV_SLEEF_ORIG_DFADD2_VF_VF2 (__riscv_vfmul (p, s, vlen),
                     VFMV_VF_F (-0.080745510756969451904f, vlen),
                     VFMV_VF_F (-1.3373665339076936258e-09f, vlen), x_hi,
                     x_lo, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (s2_hi, s2_lo, x_hi, x_lo, tmp_hi, tmp_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF2 (tmp_hi, tmp_lo,
                      VFMV_VF_F (0.78539818525314331055f, vlen),
                      VFMV_VF_F (-2.1857338617566484855e-08f, vlen), x_hi,
                      x_lo, vlen);
  DFMUL_VF2_VF2_VF (x_hi, x_lo, t_arg, tmp_hi, tmp_lo, vlen);
  VFLOAT_F rx = __riscv_vfadd (tmp_hi, tmp_lo, vlen);

  p = VFMV_VF_F (-0.2430611801e-7f, vlen);
  p = __riscv_vfmadd (p, s, VFMV_VF_F (0.3590577080e-5f, vlen), vlen);
  p = __riscv_vfmadd (p, s, VFMV_VF_F (-0.3259917721e-3f, vlen), vlen);
  RVV_SLEEF_ORIG_DFADD2_VF_VF2 (__riscv_vfmul (p, s, vlen),
                     VFMV_VF_F (0.015854343771934509277f, vlen),
                     VFMV_VF_F (4.4940051354032242811e-10f, vlen), x_hi, x_lo,
                     vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (s2_hi, s2_lo, x_hi, x_lo, tmp_hi, tmp_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF2 (tmp_hi, tmp_lo,
                      VFMV_VF_F (-0.30842512845993041992f, vlen),
                      VFMV_VF_F (-9.0728339030733922277e-09f, vlen), x_hi,
                      x_lo, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (x_hi, x_lo, s2_hi, s2_lo, tmp_hi, tmp_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF (tmp_hi, tmp_lo, VFMV_VF_F (fp_posOne_f, vlen), x_hi, x_lo,
                 vlen);
  VFLOAT_F ry = __riscv_vfadd (x_hi, x_lo, vlen);

  VBOOL_F use_rx_for_sin
      = __riscv_vmseq (__riscv_vand (I_AS_U_F (q), 2, vlen), 0, vlen);
  vsin = __riscv_vmerge (ry, rx, use_rx_for_sin, vlen);
  vcos = __riscv_vmerge (rx, ry, use_rx_for_sin, vlen);

  VUINT_F sin_sign
      = __riscv_vsll (__riscv_vand (I_AS_U_F (q), 4, vlen), BIT_WIDTH_F - 3,
                      vlen);
  VUINT_F cos_sign
      = __riscv_vsll (__riscv_vand (I_AS_U_F (__riscv_vadd (q, 2, vlen)), 4,
                                    vlen),
                      BIT_WIDTH_F - 3, vlen);
  vsin = __riscv_vfsgnjx (vsin, U_AS_F_F (sin_sign), vlen);
  vcos = __riscv_vfsgnjx (vcos, U_AS_F_F (cos_sign), vlen);

  /* Preserve exact zeros at integer and half-integer points.  */
  VBOOL_F t_zero = __riscv_vmseq (F_AS_U_F (t_arg), 0, vlen);
  VBOOL_F sinX_zero = __riscv_vmand (t_zero, use_rx_for_sin, vlen);
  VBOOL_F cosX_zero
      = __riscv_vmand (t_zero, __riscv_vmnot (use_rx_for_sin, vlen), vlen);
  VFLOAT_F signed_zero = __riscv_vfsgnj
    (VFMV_VF_F (fp_posZero_f, vlen), vx_orig, vlen);
  vsin = __riscv_vmerge (vsin, signed_zero, sinX_zero, vlen);
  vcos = __riscv_vmerge (vcos, VFMV_VF_F (fp_posZero_f, vlen), cosX_zero, vlen);

  if (nb_special_args > 0)
    {
      vsin = VMERGE_VVM_F (vsin, vsin_special, special_args, vlen);
      vcos = VMERGE_VVM_F (vcos, vcos_special, special_args, vlen);
    }

  *vsin_out = vsin;
  *vcos_out = vcos;

  RESTORE_FRM;
}

#undef LMUL
