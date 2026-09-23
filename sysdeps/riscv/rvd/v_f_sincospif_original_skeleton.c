/* Single-precision source-original SLEEF sincospif skeleton.

   This is a target-owned translation of SLEEF xsincospif_u05.  The source
   reduction, even-quotient correction, polynomial constants, double-float
   reconstruction, quadrant sign selection, negative-zero repair, large-range
   bit clearing, and infinity bit propagation are kept in source order.

   SLEEF source files are evidence only; no SLEEF source file is included,
   compiled, or linked here.  The only adaptations are the clean libmvec
   ret2 importer ABI, RVV vector types/LMUL/VL, and target-owned DF helpers.

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

static void
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT_F vx,
                                           VFLOAT_F *vsin_out,
                                           VFLOAT_F *vcos_out, size_t vlen)
{
  /* Direct source-shaped xsincospif_u05 translation.  Clean's
     tiny/large/non-finite quarantine and endpoint repair belong to the clean
     EH implementation, not to this source-original comparator.  */
  VFLOAT_F u = __riscv_vfmul (vx, 4.0f, vlen);
  VINT_F q = __riscv_vfcvt_rtz_x (u, vlen);

  /* SLEEF sinpifk's signed quotient correction:
       q = (q + ((q >> 31) ^ 1)) & ~1.  */
  VUINT_F q_u = I_AS_U_F (q);
  VUINT_F q_sign = __riscv_vsrl (q_u, 31, vlen);
  VUINT_F q_correction = __riscv_vxor
    (q_sign, VMVU_VX_F (1, vlen), vlen);
  q = __riscv_vadd (q, U_AS_I_F (q_correction), vlen);
  q = __riscv_vand (q, VMVI_VX_F (~1, vlen), vlen);

  VFLOAT_F q_f = __riscv_vfcvt_f (q, vlen);
  VFLOAT_F s = __riscv_vfsub (u, q_f, vlen);
  VFLOAT_F t = s;
  s = __riscv_vfmul (s, s, vlen);

  VFLOAT_F s2_hi, s2_lo;
  RVV_SLEEF_ORIG_DFMUL_VF_VF (t, t, s2_hi, s2_lo, vlen);

  /* Source first kernel: sine reconstruction.  SLEEF evaluates this
     complete kernel before evaluating the independent cosine kernel; it
     does not select coefficient families per lane.  */
  VFLOAT_F p = VFMV_VF_F (0.3093842054e-6f, vlen);
  p = __riscv_vfmadd
    (p, s, VFMV_VF_F (-0.3657307388e-4f, vlen), vlen);
  p = __riscv_vfmadd
    (p, s, VFMV_VF_F (+0.2490393585e-2f, vlen), vlen);

  VFLOAT_F x_hi, x_lo;
  RVV_SLEEF_ORIG_DFADD2_VF_VF2
     (__riscv_vfmul (p, s, vlen),
     VFMV_VF_F (-0.080745510756969451904f, vlen),
     VFMV_VF_F (-1.3373665339076936258e-09f, vlen),
     x_hi, x_lo, vlen);

  VFLOAT_F tmp_hi, tmp_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (s2_hi, s2_lo, x_hi, x_lo, tmp_hi, tmp_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF2
    (tmp_hi, tmp_lo,
     VFMV_VF_F (+0.78539818525314331055f, vlen),
     VFMV_VF_F (-2.1857338617566484855e-08f, vlen),
     x_hi, x_lo, vlen);

  /* Source first result: x *= t, then collapse and repair -0 before the
     quadrant pair selection.  */
  DFMUL_VF2_VF2_VF (x_hi, x_lo, t, tmp_hi, tmp_lo, vlen);
  VFLOAT_F rx = __riscv_vfadd (tmp_hi, tmp_lo, vlen);
  VBOOL_F negative_zero = __riscv_vmseq
    (F_AS_U_F (vx), VMVU_VX_F (0x80000000u, vlen), vlen);
  rx = __riscv_vmerge (rx, VFMV_VF_F (-0.0f, vlen), negative_zero, vlen);

  /* Source second kernel: cosine reconstruction, in its original order.  */
  p = VFMV_VF_F (-0.2430611801e-7f, vlen);
  p = __riscv_vfmadd
    (p, s, VFMV_VF_F (+0.3590577080e-5f, vlen), vlen);
  p = __riscv_vfmadd
    (p, s, VFMV_VF_F (-0.3259917721e-3f, vlen), vlen);

  RVV_SLEEF_ORIG_DFADD2_VF_VF2
    (__riscv_vfmul (p, s, vlen),
     VFMV_VF_F (+0.015854343771934509277f, vlen),
     VFMV_VF_F (+4.4940051354032242811e-10f, vlen),
     x_hi, x_lo, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (s2_hi, s2_lo, x_hi, x_lo, tmp_hi, tmp_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF2
    (tmp_hi, tmp_lo,
     VFMV_VF_F (-0.30842512845993041992f, vlen),
     VFMV_VF_F (-9.0728339030733922277e-09f, vlen),
     x_hi, x_lo, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (x_hi, x_lo, s2_hi, s2_lo, tmp_hi, tmp_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF
    (tmp_hi, tmp_lo, VFMV_VF_F (1.0f, vlen), x_hi, x_lo, vlen);
  VFLOAT_F ry = __riscv_vfadd (x_hi, x_lo, vlen);

  /* Select the pair.  */
  VBOOL_F use_rx_for_sin = __riscv_vmseq
    (__riscv_vand (q, 2, vlen), 0, vlen);
  VFLOAT_F vsin = __riscv_vmerge (ry, rx, use_rx_for_sin, vlen);
  VFLOAT_F vcos = __riscv_vmerge (rx, ry, use_rx_for_sin, vlen);

  /* Source xsincospif_u05 applies the q&4 and (q+2)&4 sign flips to the
     returned pair after quadrant selection.  */
  VUINT_F q_plus2 = I_AS_U_F (__riscv_vadd (q, 2, vlen));
  VBOOL_F flip_sin = __riscv_vmseq
    (__riscv_vand (q, 4, vlen), 4, vlen);
  VBOOL_F flip_cos = __riscv_vmseq
    (__riscv_vand (q_plus2, 4, vlen), 4, vlen);
  VUINT_F sin_bits = F_AS_U_F (vsin);
  VUINT_F cos_bits = F_AS_U_F (vcos);
  VUINT_F sign_bits = VMVU_VX_F (0x80000000u, vlen);
  sin_bits = __riscv_vxor
    (sin_bits, __riscv_vmerge (VMVU_VX_F (0, vlen), sign_bits,
                               flip_sin, vlen), vlen);
  cos_bits = __riscv_vxor
    (cos_bits, __riscv_vmerge (VMVU_VX_F (0, vlen), sign_bits,
                               flip_cos, vlen), vlen);

  vsin = U_AS_F_F (sin_bits);
  vcos = U_AS_F_F (cos_bits);

  /* Match the source's final clear-and-OR sequence exactly.  For finite
     |x| > 1e7 this clears both result bit patterns; for infinity the later
     bitwise OR turns both results into the source all-ones NaN pattern.
     NaN lanes are not selected by visinf and therefore retain the NaN
     produced by the source arithmetic.  */
  VBOOL_F large = __riscv_vmfgt
    (__riscv_vfabs (vx, vlen), RVV_SLEEF_TRIGPI_RANGE_MAX_F, vlen);
  VUINT_F zero_bits = VMVU_VX_F (0, vlen);
  sin_bits = __riscv_vmerge (F_AS_U_F (vsin), zero_bits, large, vlen);
  cos_bits = __riscv_vmerge (F_AS_U_F (vcos), zero_bits, large, vlen);

  VBOOL_F input_inf = __riscv_vmfeq
    (__riscv_vfabs (vx, vlen), fp_posInf_f, vlen);
  VUINT_F inf_bits = __riscv_vmerge
    (zero_bits, VMVU_VX_F (~0u, vlen), input_inf, vlen);
  sin_bits = __riscv_vor (sin_bits, inf_bits, vlen);
  cos_bits = __riscv_vor (cos_bits, inf_bits, vlen);

  vsin = U_AS_F_F (sin_bits);
  vcos = U_AS_F_F (cos_bits);

  *vsin_out = vsin;
  *vcos_out = vcos;
}

#undef LMUL
