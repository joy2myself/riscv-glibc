/* RISC-V RVV based strcmp.
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

#include <sysdep.h>
#include <sys/asm.h>

#ifndef STRCMP
# define STRCMP strcmp
#endif

#define result a0

#define str1 a0
#define str2 a1

#define ivl a2
#define temp1 a3
#define temp2 a4

#define vstr1 v0
#define vstr2 v8
#define vmask1 v16
#define vmask2 v17

ENTRY (STRCMP)
.option push
.option arch, +v
    li t2, 16
.Lcmp_scalar_prefix:
    lbu t0, 0(a0)
    lbu t1, 0(a1)
    bne t0, t1, .Lcmp_first
    beqz t0, .Lcmp_first
    addi a0, a0, 1
    addi a1, a1, 1
    addi t2, t2, -1
    bnez t2, .Lcmp_scalar_prefix
.Lcmp_loop:
    vsetvli a2, zero, e8, m4, ta, ma
    vle8ff.v v0, (a0)
    vle8ff.v v8, (a1)
    vmseq.vi v16, v0, 0
    vmsne.vv v17, v0, v8
    vmor.mm v16, v16, v17
    vfirst.m a3, v16
    bgez a3, .Lcmp_found
    csrr a2, vl
    add a0, a0, a2
    add a1, a1, a2
    j .Lcmp_loop
.Lcmp_found:
    add a0, a0, a3
    add a1, a1, a3
    lbu t0, 0(a0)
    lbu t1, 0(a1)
.Lcmp_first:
    sub a0, t0, t1
    ret
.option pop
END (STRCMP)
libc_hidden_builtin_def (strcmp)
