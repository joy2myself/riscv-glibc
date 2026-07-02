/* RISC-V RVV based memrchr.
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

#ifndef MEMRCHR
# define MEMRCHR __memrchr
#endif

#define src a0
#define ch a1
#define num a2
#define vl a3
#define first_idx a4
#define idx a5

#define vmask v0
#define vdata v8
#define vidx v16
#define vmax v24

ENTRY (MEMRCHR)
.option push
.option arch, +v
    add src, src, num
L(loop):
    beqz num, L(nohit)

    vsetvli vl, num, e8, m4, ta, ma
    sub src, src, vl
    vle8.v vdata, (src)
    sub num, num, vl

    vmseq.vx vmask, vdata, ch
    vfirst.m first_idx, vmask
    bltz first_idx, L(loop)
L(hit):
    vsetvli zero, vl, e16, m8, ta, ma
    vid.v vidx

    /* reduce max index over hit lanes */
    vmv.s.x vmax, x0
    vredmaxu.vs vmax, vidx, vmax, vmask.t
    vmv.x.s idx, vmax

    add a0, src, idx
    ret
L(nohit):
    li a0, 0
    ret
.option pop
END (MEMRCHR)
libc_hidden_def (MEMRCHR)
weak_alias (MEMRCHR, memrchr)
libc_hidden_builtin_def (memrchr)
