/* Copyright (C) 2026 Free Software Foundation, Inc.
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

#ifndef _VECMATH_CONFIG_H
#define _VECMATH_CONFIG_H

#include <stdint.h>

__BEGIN_DECLS

#define __VCFG_PASTE4_BASE(A, B, C, D) A##B##C##D
#define __VCFG_PASTE4(A, B, C, D) __VCFG_PASTE4_BASE (A, B, C, D)
#define __VCFG_PASTE5_BASE(A, B, C, D, E) A##B##C##D##E
#define __VCFG_PASTE5(A, B, C, D, E) __VCFG_PASTE5_BASE (A, B, C, D, E)

#define __MAKE_VMERGE_VVM(T) \
  __VCFG_PASTE5 (__riscv_vmerge_vvm_, T, BIT_WIDTH, m, LMUL)
#define __MAKE_VMERGE_VXM(T) \
  __VCFG_PASTE5 (__riscv_vmerge_vxm_, T, BIT_WIDTH, m, LMUL)
#define __MAKE_VFMERGE_VFM \
  __VCFG_PASTE4 (__riscv_vfmerge_vfm_f, BIT_WIDTH, m, LMUL)
#define __MAKE_VFMUL_VF_M \
  __VCFG_PASTE5 (__riscv_vfmul_vf_f, BIT_WIDTH, m, LMUL, _m)
#define __MAKE_VFREC7_M \
  __VCFG_PASTE5 (__riscv_vfrec7_v_f, BIT_WIDTH, m, LMUL, _m)
#define __MAKE_VFADD_VV_M \
  __VCFG_PASTE5 (__riscv_vfadd_vv_f, BIT_WIDTH, m, LMUL, _m)
#define __MAKE_VFADD_VF_M \
  __VCFG_PASTE5 (__riscv_vfadd_vf_f, BIT_WIDTH, m, LMUL, _m)
#define __MAKE_VFRSUB_VF_M \
  __VCFG_PASTE5 (__riscv_vfrsub_vf_f, BIT_WIDTH, m, LMUL, _m)
#define __MAKE_VFSUB_VV_M \
  __VCFG_PASTE5 (__riscv_vfsub_vv_f, BIT_WIDTH, m, LMUL, _m)
#define __MAKE_VFSUB_VF_M \
  __VCFG_PASTE5 (__riscv_vfsub_vf_f, BIT_WIDTH, m, LMUL, _m)
#define __MAKE_VFSUB_VV \
  __VCFG_PASTE4 (__riscv_vfsub_vv_f, BIT_WIDTH, m, LMUL)
#define __MAKE_VFSUB_VF \
  __VCFG_PASTE4 (__riscv_vfsub_vf_f, BIT_WIDTH, m, LMUL)
#define __MAKE_VFADD_VF \
  __VCFG_PASTE4 (__riscv_vfadd_vf_f, BIT_WIDTH, m, LMUL)
#define __MAKE_VFMADD_VV \
  __VCFG_PASTE4 (__riscv_vfmadd_vv_f, BIT_WIDTH, m, LMUL)
#define __MAKE_VFMACC_VV \
  __VCFG_PASTE4 (__riscv_vfmacc_vv_f, BIT_WIDTH, m, LMUL)
#define __MAKE_VFMUL_VV \
  __VCFG_PASTE4 (__riscv_vfmul_vv_f, BIT_WIDTH, m, LMUL)
#define __MAKE_VFMUL_VV_M \
  __VCFG_PASTE5 (__riscv_vfmul_vv_f, BIT_WIDTH, m, LMUL, _m)
#define __MAKE_VFMSUB_VV_M \
  __VCFG_PASTE5 (__riscv_vfmsub_vv_f, BIT_WIDTH, m, LMUL, _m)
#define __MAKE_VFMSUB_VF_M \
  __VCFG_PASTE5 (__riscv_vfmsub_vf_f, BIT_WIDTH, m, LMUL, _m)
#define __MAKE_VFMADD_VV_M \
  __VCFG_PASTE5 (__riscv_vfmadd_vv_f, BIT_WIDTH, m, LMUL, _m)
#define __MAKE_VFNMSUB_VV_M \
  __VCFG_PASTE5 (__riscv_vfnmsub_vv_f, BIT_WIDTH, m, LMUL, _m)
#define __MAKE_VFMACC_VV_M \
  __VCFG_PASTE5 (__riscv_vfmacc_vv_f, BIT_WIDTH, m, LMUL, _m)
#define __MAKE_VFMACC_VF_M \
  __VCFG_PASTE5 (__riscv_vfmacc_vf_f, BIT_WIDTH, m, LMUL, _m)
#define __MAKE_VFDIV_VV_M \
  __VCFG_PASTE5 (__riscv_vfdiv_vv_f, BIT_WIDTH, m, LMUL, _m)
#define __MAKE_VFNMSAC_VV_M \
  __VCFG_PASTE5 (__riscv_vfnmsac_vv_f, BIT_WIDTH, m, LMUL, _m)
#define __MAKE_VFRDIV_VF_M \
  __VCFG_PASTE5 (__riscv_vfrdiv_vf_f, BIT_WIDTH, m, LMUL, _m)
#define __MAKE_VAND_VV_I \
  __VCFG_PASTE4 (__riscv_vand_vv_i, BIT_WIDTH, m, LMUL)
#define __MAKE_VXOR_VV \
  __VCFG_PASTE4 (__riscv_vxor_vv_u, BIT_WIDTH, m, LMUL)

#define VMERGE_VVM_F(op1, op2, mask, vlen) \
  __MAKE_VMERGE_VVM (f) ((op1), (op2), (mask), (vlen))
#define VMERGE_VVM_I(op1, op2, mask, vlen) \
  __MAKE_VMERGE_VVM (i) ((op1), (op2), (mask), (vlen))
#define VMERGE_VXM_I(op1, scalar, mask, vlen) \
  __MAKE_VMERGE_VXM (i) ((op1), (scalar), (mask), (vlen))
#define VFMERGE_VFM(op1, scalar, mask, vlen) \
  __MAKE_VFMERGE_VFM ((op1), (scalar), (mask), (vlen))
#define VFMUL_VF_M(mask, op1, scalar, vlen) \
  __MAKE_VFMUL_VF_M ((mask), (op1), (scalar), (vlen))
#define VFREC7_M(mask, op1, vlen) \
  __MAKE_VFREC7_M ((mask), (op1), (vlen))
#define VFADD_VV_M(mask, op1, op2, vlen) \
  __MAKE_VFADD_VV_M ((mask), (op1), (op2), (vlen))
#define VFADD_VF_M(mask, op1, scalar, vlen) \
  __MAKE_VFADD_VF_M ((mask), (op1), (scalar), (vlen))
#define VFRSUB_VF_M(mask, op1, scalar, vlen) \
  __MAKE_VFRSUB_VF_M ((mask), (op1), (scalar), (vlen))
#define VFSUB_VV_M(mask, op1, op2, vlen) \
  __MAKE_VFSUB_VV_M ((mask), (op1), (op2), (vlen))
#define VFSUB_VF_M(mask, op1, scalar, vlen) \
  __MAKE_VFSUB_VF_M ((mask), (op1), (scalar), (vlen))
#define VFSUB_VV(op1, op2, vlen) \
  __MAKE_VFSUB_VV ((op1), (op2), (vlen))
#define VFSUB_VF(op1, scalar, vlen) \
  __MAKE_VFSUB_VF ((op1), (scalar), (vlen))
#define VFADD_VF(op1, scalar, vlen) \
  __MAKE_VFADD_VF ((op1), (scalar), (vlen))
#define VFMADD_VV(op1, op2, op3, vlen) \
  __MAKE_VFMADD_VV ((op1), (op2), (op3), (vlen))
#define VFMACC_VV(op1, op2, op3, vlen) \
  __MAKE_VFMACC_VV ((op1), (op2), (op3), (vlen))
#define VFMUL_VV(op1, op2, vlen) \
  __MAKE_VFMUL_VV ((op1), (op2), (vlen))
#define VFMUL_VV_M(mask, op1, op2, vlen) \
  __MAKE_VFMUL_VV_M ((mask), (op1), (op2), (vlen))
#define VFMSUB_VV_M(mask, op1, op2, op3, vlen) \
  __MAKE_VFMSUB_VV_M ((mask), (op1), (op2), (op3), (vlen))
#define VFMSUB_VF_M(mask, op1, scalar, op3, vlen) \
  __MAKE_VFMSUB_VF_M ((mask), (op1), (scalar), (op3), (vlen))
#define VFMADD_VV_M(mask, op1, op2, op3, vlen) \
  __MAKE_VFMADD_VV_M ((mask), (op1), (op2), (op3), (vlen))
#define VFNMSUB_VV_M(mask, op1, op2, op3, vlen) \
  __MAKE_VFNMSUB_VV_M ((mask), (op1), (op2), (op3), (vlen))
#define VFMACC_VV_M(mask, op1, op2, op3, vlen) \
  __MAKE_VFMACC_VV_M ((mask), (op1), (op2), (op3), (vlen))
#define VFMACC_VF_M(mask, op1, scalar, op2, vlen) \
  __MAKE_VFMACC_VF_M ((mask), (op1), (scalar), (op2), (vlen))
#define VFDIV_VV_M(mask, op1, op2, vlen) \
  __MAKE_VFDIV_VV_M ((mask), (op1), (op2), (vlen))
#define VFNMSAC_VV_M(mask, op1, op2, op3, vlen) \
  __MAKE_VFNMSAC_VV_M ((mask), (op1), (op2), (op3), (vlen))
#define VFRDIV_VF_M(mask, op1, scalar, vlen) \
  __MAKE_VFRDIV_VF_M ((mask), (op1), (scalar), (vlen))
#define VAND_VV_I(op1, op2, vlen) \
  __MAKE_VAND_VV_I ((op1), (op2), (vlen))
#define VXOR_VV(op1, op2, vlen) \
  __MAKE_VXOR_VV ((op1), (op2), (vlen))

#define RVV_VLUXEI32_F(base, index, vlen) \
  __VCFG_PASTE4 (__riscv_vluxei32_v_f, BIT_WIDTH, m, LMUL) \
    ((base), (index), (vlen))

#define VMOR(op1, op2, vlen) \
  __riscv_vmor ((op1), (op2), (vlen))

#define RVV_VFCVT_X_F_RNE(value, vlen) \
  __VCFG_PASTE5 (__riscv_vfcvt_x_f_v_i, BIT_WIDTH, m, LMUL, _rm) \
    ((value), __RISCV_FRM_RNE, (vlen))

/* Explicitly masked float/integer conversions used by the double atan2
   fixed-point core.  Keep the intrinsic spelling in the common framework so
   the routine skeleton contains only its operation graph.  */
#define __MAKE_VFCVT_X_F_M \
  __VCFG_PASTE5 (__riscv_vfcvt_x_f_v_i, BIT_WIDTH, m, LMUL, _m)
#define __MAKE_VFCVT_F_X_M \
  __VCFG_PASTE5 (__riscv_vfcvt_f_x_v_f, BIT_WIDTH, m, LMUL, _m)
#define VFCVT_X_F_M(mask, value, vlen) \
  __MAKE_VFCVT_X_F_M ((mask), (value), (vlen))
#define VFCVT_F_X_M(mask, value, vlen) \
  __MAKE_VFCVT_F_X_M ((mask), (value), (vlen))

#ifndef FE_TONEAREST
#define FE_TONEAREST 0x000
#endif

#define read_frm()                                                            \
  ({                                                                          \
    unsigned long __value;                                                    \
    __asm__ __volatile__ ("frrm %0" : "=r"(__value)::"memory");               \
    __value;                                                                  \
  })

#define write_frm(value)                                                      \
  ({                                                                          \
    unsigned long __value;                                                    \
    __asm__ __volatile__ ("fsrm %0, %1"                                       \
                          : "=r"(__value)                                     \
                          : "r"(value)                                        \
                          : "memory");                                        \
    __value;                                                                  \
  })

#define SET_ROUNDTONEAREST                                                    \
  int __original_frm = read_frm ();                                           \
  if (__original_frm != FE_TONEAREST)                                         \
  {                                                                           \
    write_frm (FE_TONEAREST);                                                 \
  }

#define RESTORE_FRM                                                           \
  do                                                                          \
  {                                                                           \
    if (__original_frm != FE_TONEAREST)                                       \
    {                                                                         \
      write_frm (__original_frm);                                             \
    }                                                                         \
  }                                                                           \
  while (0)

#define PSTEP(coeff_j, x, poly, vlen)                                         \
  __riscv_vfmadd ((poly), (x), VFMV_VF ((coeff_j), (vlen)), (vlen))

#define PSTEP_I(coeff_j, x, poly, vlen)                                       \
  __riscv_vsadd (__riscv_vsmul ((poly), (x), 1, (vlen)), (coeff_j), (vlen))

#define PSTEP_I_SLL(coeff_j, x, k, poly, vlen)                                \
  __riscv_vsadd (                                                             \
      __riscv_vsll (__riscv_vsmul ((poly), (x), 1, (vlen)), (k), (vlen)),     \
      (coeff_j), (vlen))

#define PSTEP_I_SRA(coeff_j, x, k, poly, vlen)                                \
  __riscv_vsadd (                                                             \
      __riscv_vsra (__riscv_vsmul ((poly), (x), 1, (vlen)), (k), (vlen)),     \
      (coeff_j), (vlen))

#define PSTEPN_I(coeff_j, x, poly, vlen)                                      \
  __riscv_vrsub (__riscv_vsmul ((poly), (x), 1, (vlen)), (coeff_j), (vlen))

#define PSTEP_I_HI(coeff_j, x, poly, vlen)                                    \
  __riscv_vadd (__riscv_vmulh ((poly), (x), (vlen)), (coeff_j), (vlen))

#define PSTEP_I_ab(pick_a, coeff_a, coeff_b, x, poly, vlen)                   \
  __riscv_vsadd (__riscv_vsmul ((poly), (x), 1, (vlen)),                      \
                 VMERGE_VXM_I (VMVI_VX ((coeff_b), (vlen)), (coeff_a),       \
                               (pick_a), (vlen)),                             \
                 (vlen))

#define IDENTIFY(vclass, stencil, identity_mask, vlen)                        \
  identity_mask = __riscv_vmsgtu (__riscv_vand ((vclass), (stencil), (vlen)), \
                                  0, (vlen))

#define FCLIP(x, min, max, vlen)                                              \
  __riscv_vfmax (__riscv_vfmin ((x), (max), (vlen)), (min), (vlen))

#define VSRL_I_AS_U(x, nbits, vlen)                                           \
  U_AS_I (__riscv_vsrl (I_AS_U ((x)), (nbits), (vlen)))

#define FAST2SUM(x, y, sum_hi, sum_lo, vlen)                                  \
  do                                                                          \
    {                                                                         \
      (sum_hi) = __riscv_vfadd ((x), (y), (vlen));                            \
      (sum_lo) = __riscv_vfsub ((x), (sum_hi), (vlen));                       \
      (sum_lo) = __riscv_vfadd ((sum_lo), (y), (vlen));                       \
    }                                                                         \
  while (0)

#define KNUTH2SUM(x, y, sum_hi, sum_lo, vlen)                                 \
  do                                                                          \
    {                                                                         \
      (sum_hi) = __riscv_vfadd ((x), (y), (vlen));                            \
      VFLOAT _x_hat = __riscv_vfsub ((sum_hi), (y), (vlen));                  \
      (sum_lo) = __riscv_vfadd (                                              \
          __riscv_vfsub ((x), _x_hat, (vlen)),                                \
          __riscv_vfsub ((y), __riscv_vfsub ((sum_hi), _x_hat, (vlen)),       \
                         (vlen)),                                             \
          (vlen));                                                            \
    }                                                                         \
  while (0)

#define FIX2FLT(x, scale, y_hi, y_lo, vlen)                                   \
  do                                                                          \
    {                                                                         \
      (y_hi) = __riscv_vfcvt_f ((x), (vlen));                                 \
      (y_lo) = __riscv_vfcvt_f (                                              \
          __riscv_vsub ((x), __riscv_vfcvt_x ((y_hi), (vlen)), (vlen)),       \
          (vlen));                                                            \
      (y_hi) = __riscv_vfmul ((y_hi), (scale), (vlen));                       \
      (y_lo) = __riscv_vfmul ((y_lo), (scale), (vlen));                       \
    }                                                                         \
  while (0)

#define FLT2FIX(x_hi, x_lo, scale, y, vlen)                                   \
  do                                                                          \
    {                                                                         \
      (y) = __riscv_vfcvt_x (__riscv_vfmul ((x_hi), (scale), (vlen)),         \
                             (vlen));                                         \
      (y) = __riscv_vadd (                                                    \
          (y),                                                                \
          __riscv_vfcvt_x (__riscv_vfmul ((x_lo), (scale), (vlen)), (vlen)),  \
          (vlen));                                                            \
    }                                                                         \
  while (0)

#define SQR_X2(x_hi, x_lo, prod_hi, prod_lo, vlen)                            \
  do                                                                          \
    {                                                                         \
      (prod_hi) = __riscv_vfmul ((x_hi), (x_hi), (vlen));                     \
      (prod_lo) = __riscv_vfmsub ((x_hi), (x_hi), (prod_hi), (vlen));         \
      (prod_lo) = __riscv_vfmacc ((prod_lo), (x_hi), (x_lo), (vlen));         \
      (prod_lo) = __riscv_vfmacc ((prod_lo), (x_lo), (x_hi), (vlen));         \
    }                                                                         \
  while (0)

#define FAST_LDEXP(num, exp, vlen)                                            \
  do                                                                          \
    {                                                                         \
      VINT _n1 = __riscv_vsra ((exp), 1, (vlen));                             \
      VINT _n2 = __riscv_vsub ((exp), _n1, (vlen));                           \
      _n1 = __riscv_vsll (_n1, MAN_LEN, (vlen));                              \
      (num) = I_AS_F (__riscv_vadd (F_AS_I ((num)), _n1, (vlen)));            \
      _n2 = __riscv_vadd (_n2, EXP_BIAS, (vlen));                             \
      _n2 = __riscv_vsll (_n2, MAN_LEN, (vlen));                              \
      (num) = __riscv_vfmul ((num), I_AS_F (_n2), (vlen));                    \
    }                                                                         \
  while (0)

/* Define the table for table-driven implementations.  Scalar constants live
   in v_math.h.  */
#define V_LOG_TABLE_BITS 7
extern const struct v_log_data
{
  const int64_t table[1 << V_LOG_TABLE_BITS];
} __v_log_data attribute_hidden;

/* Shared constants and tables for the double-precision exp2m1 implementation.  */
#define V_EXP2M1_TABLE_BITS 7
#define V_EXP2M1_N (1 << V_EXP2M1_TABLE_BITS)
extern const struct v_exp2m1_data
{
  const uint64_t exp2_tab[2 * V_EXP2M1_N];
  const uint64_t scalem1_tab[88];
} __v_exp2m1_data attribute_hidden;
/* Indexed coefficient tables for the double-precision lgamma implementation.  */
extern const struct v_lgamma_data
{
  const double neg_bnd_to_2_c[26][2];
  const double neg_3_to_bnd_c[27][2];
} __v_lgamma_data attribute_hidden;

/* Lookup tables for the double-precision exp10m1 implementation.  */
extern const struct v_exp10m1_data
{
  const uint64_t exp_data[128];
  const uint64_t exp_tail[128];
  const uint64_t scalem1[88];
} __v_exp10m1_data attribute_hidden;

/* Shared SLEEF rempif table for the single-precision trigonometric paths.  */
#define V_SLEEF_REMPITABSP_LEN 416
extern const struct v_sleef_trig_data
{
  const float rempitabsp[V_SLEEF_REMPITABSP_LEN];
} __v_sleef_trig_data attribute_hidden;

/* Large-argument reduction table for the double-precision cos
   implementation.  */
extern const struct v_cos_data
{
  const double two_over_pi[28];
} __v_cos_data attribute_hidden;

#define PSTEP_ab(pick_a, coeff_a, coeff_b, x, poly, vlen)                     \
  __riscv_vfmadd ((poly), (x),                                                \
                  VFMERGE_VFM (VFMV_VF ((coeff_b), (vlen)), (coeff_a),        \
                               (pick_a), (vlen)),                             \
                  (vlen))

#define FAST2SUM(x, y, sum_hi, sum_lo, vlen)                                  \
  do                                                                          \
    {                                                                         \
      (sum_hi) = __riscv_vfadd ((x), (y), (vlen));                            \
      (sum_lo) = __riscv_vfsub ((x), (sum_hi), (vlen));                       \
      (sum_lo) = __riscv_vfadd ((sum_lo), (y), (vlen));                       \
    }                                                                         \
  while (0)

#define POS2SUM(x, y, sum_hi, sum_lo, vlen)                                   \
  do                                                                          \
    {                                                                         \
      VFLOAT _first = __riscv_vfmax ((x), (y), (vlen));                       \
      VFLOAT _second = __riscv_vfmin ((x), (y), (vlen));                      \
      (sum_hi) = __riscv_vfadd ((x), (y), (vlen));                            \
      (sum_lo) = __riscv_vfadd (                                              \
          __riscv_vfsub (_first, (sum_hi), (vlen)), _second, (vlen));         \
    }                                                                         \
  while (0)

#define PROD_X1Y1(x, y, prod_hi, prod_lo, vlen)                               \
  do                                                                          \
    {                                                                         \
      (prod_hi) = __riscv_vfmul ((x), (y), (vlen));                           \
      (prod_lo) = __riscv_vfmsub ((x), (y), (prod_hi), (vlen));               \
    }                                                                         \
  while (0)

#define PROD_X2Y2(x_hi, x_lo, y_hi, y_lo, prod_hi, prod_lo, vlen)             \
  do                                                                          \
    {                                                                         \
      (prod_hi) = __riscv_vfmul ((x_hi), (y_hi), (vlen));                     \
      (prod_lo) = __riscv_vfmsub ((x_hi), (y_hi), (prod_hi), (vlen));         \
      (prod_lo) = __riscv_vfmacc ((prod_lo), (x_hi), (y_lo), (vlen));         \
      (prod_lo) = __riscv_vfmacc ((prod_lo), (x_lo), (y_hi), (vlen));         \
    }                                                                         \
  while (0)

#define DIV_N1D2(numer, denom, delta_d, q_hi, q_lo, vlen)                     \
  do                                                                          \
    {                                                                         \
      (q_hi) = __riscv_vfdiv ((numer), (denom), (vlen));                      \
      (q_lo) = __riscv_vfnmsub ((q_hi), (denom), (numer), (vlen));            \
      (q_lo) = __riscv_vfnmsac ((q_lo), (q_hi), (delta_d), (vlen));           \
      (q_lo) = __riscv_vfmul ((q_lo), __riscv_vfrec7 ((denom), (vlen)),       \
                              (vlen));                                         \
    }                                                                         \
  while (0)

#define DIV_N2D2(numer, delta_n, denom, delta_d, q, vlen)                     \
  do                                                                          \
    {                                                                         \
      VFLOAT _tmp;                                                            \
      (q) = __riscv_vfdiv ((numer), (denom), (vlen));                         \
      _tmp = __riscv_vfnmsub ((q), (denom), (numer), (vlen));                 \
      _tmp = __riscv_vfnmsac (_tmp, (q), (delta_d), (vlen));                  \
      _tmp = __riscv_vfadd (_tmp, (delta_n), (vlen));                         \
      _tmp = __riscv_vfmul (_tmp, __riscv_vfrec7 ((denom), (vlen)), (vlen));  \
      (q) = __riscv_vfadd ((q), _tmp, (vlen));                                \
    }                                                                         \
  while (0)

#define DIV2_N2D2(numer, delta_n, denom, delta_d, q_hi, q_lo, vlen)          \
  do                                                                          \
    {                                                                         \
      VFLOAT _q;                                                              \
      (q_hi) = __riscv_vfdiv ((numer), (denom), (vlen));                      \
      _q = __riscv_vfnmsub ((q_hi), (denom), (numer), (vlen));                \
      _q = __riscv_vfnmsac (_q, (q_hi), (delta_d), (vlen));                   \
      _q = __riscv_vfadd (_q, (delta_n), (vlen));                             \
      (q_lo) = __riscv_vfmul (_q, __riscv_vfrec7 ((denom), (vlen)), (vlen));  \
    }                                                                         \
  while (0)

#define __MAKE_VFADD_VV \
  __VCFG_PASTE4 (__riscv_vfadd_vv_f, BIT_WIDTH, m, LMUL)
#define __MAKE_VFSUB_VV \
  __VCFG_PASTE4 (__riscv_vfsub_vv_f, BIT_WIDTH, m, LMUL)
#define __MAKE_VFMUL_VV \
  __VCFG_PASTE4 (__riscv_vfmul_vv_f, BIT_WIDTH, m, LMUL)
#define __MAKE_VFDIV_VV \
  __VCFG_PASTE4 (__riscv_vfdiv_vv_f, BIT_WIDTH, m, LMUL)
#define __MAKE_VFMADD_VV \
  __VCFG_PASTE4 (__riscv_vfmadd_vv_f, BIT_WIDTH, m, LMUL)
#define __MAKE_VFRSUB_VF \
  __VCFG_PASTE4 (__riscv_vfrsub_vf_f, BIT_WIDTH, m, LMUL)
#define __MAKE_VFMUL_VF \
  __VCFG_PASTE4 (__riscv_vfmul_vf_f, BIT_WIDTH, m, LMUL)
#define __MAKE_VFMSUB_VV \
  __VCFG_PASTE4 (__riscv_vfmsub_vv_f, BIT_WIDTH, m, LMUL)
#define __MAKE_VFMACC_VV \
  __VCFG_PASTE4 (__riscv_vfmacc_vv_f, BIT_WIDTH, m, LMUL)
#define __MAKE_VFNMSUB_VV \
  __VCFG_PASTE4 (__riscv_vfnmsub_vv_f, BIT_WIDTH, m, LMUL)
#define __MAKE_VFNMSAC_VV \
  __VCFG_PASTE4 (__riscv_vfnmsac_vv_f, BIT_WIDTH, m, LMUL)
#define __MAKE_VFNEG_V \
  __VCFG_PASTE4 (__riscv_vfneg_v_f, BIT_WIDTH, m, LMUL)

#define VFADD_VV(op1, op2, vlen) \
  __MAKE_VFADD_VV ((op1), (op2), (vlen))
#define VFRSUB_VF(op1, scalar, vlen) \
  __MAKE_VFRSUB_VF ((op1), (scalar), (vlen))
#define VFSUB_VV(op1, op2, vlen) \
  __MAKE_VFSUB_VV ((op1), (op2), (vlen))
#define VFMUL_VV(op1, op2, vlen) \
  __MAKE_VFMUL_VV ((op1), (op2), (vlen))
#define VFMUL_VF(op1, scalar, vlen) \
  __MAKE_VFMUL_VF ((op1), (scalar), (vlen))
#define VFDIV_VV(op1, op2, vlen) \
  __MAKE_VFDIV_VV ((op1), (op2), (vlen))
#define VFMADD_VV(op1, op2, op3, vlen) \
  __MAKE_VFMADD_VV ((op1), (op2), (op3), (vlen))
#define VFMACC_VV(op1, op2, op3, vlen) \
  __MAKE_VFMACC_VV ((op1), (op2), (op3), (vlen))
#define VFMSUB_VV(op1, op2, op3, vlen) \
  __MAKE_VFMSUB_VV ((op1), (op2), (op3), (vlen))
#define VFNMSUB_VV(op1, op2, op3, vlen) \
  __MAKE_VFNMSUB_VV ((op1), (op2), (op3), (vlen))
#define VFNMSAC_VV(op1, op2, op3, vlen) \
  __MAKE_VFNMSAC_VV ((op1), (op2), (op3), (vlen))
#define VFNEG_V(op1, vlen) \
  __MAKE_VFNEG_V ((op1), (vlen))

__END_DECLS

#endif
