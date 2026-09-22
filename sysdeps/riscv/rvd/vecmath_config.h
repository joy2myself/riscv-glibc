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

#define VMERGE_VVM_F(op1, op2, mask, vlen) \
  __MAKE_VMERGE_VVM (f) ((op1), (op2), (mask), (vlen))
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

#define IDENTIFY(vclass, stencil, identity_mask, vlen)                        \
  identity_mask = __riscv_vmsgtu (__riscv_vand ((vclass), (stencil), (vlen)), \
                                  0, (vlen))

/* Define the table for table-driven implementations. */
#define V_LOG_TABLE_BITS 7
extern const struct v_log_data
{
  const int64_t table[1 << V_LOG_TABLE_BITS];
} __v_log_data attribute_hidden;

__END_DECLS

#endif
