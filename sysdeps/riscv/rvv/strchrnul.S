/* RISC-V RVV based strchrnul.
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

#ifndef STRCHRNUL
# define STRCHRNUL __strchrnul
#endif

#define str a0
#define ch a1
#define end_offset a2
#define ch_offset a3
#define temp1 a4
#define temp2 a5
#define cur_vl a6
#define ivl t0

#define ELEM_LMUL_SETTING m1
#define vstr v0
#define vmask_end v8
#define vmask_ch v9

ENTRY (STRCHRNUL)
.option push
.option arch, +v
    andi a1, a1, 255
    lbu t1, 0(a0)
    beq t1, a1, .Lscan_return
    beqz t1, .Lscan_end
.Lscan_loop:
    vsetvli t0, zero, e8, m4, ta, ma
    vle8ff.v v0, (a0)
    vmseq.vi v8, v0, 0
    vmseq.vx v9, v0, a1
    vmor.mm v8, v8, v9
    vfirst.m a2, v8
    bgez a2, .Lscan_found
    csrr a6, vl
    add a0, a0, a6
    j .Lscan_loop
.Lscan_found:
    add a0, a0, a2
    lbu t1, 0(a0)
    beq t1, a1, .Lscan_return
.Lscan_end:

.Lscan_return:
    ret
.option pop
END (STRCHRNUL)

#ifndef STRCHRNUL_NO_ALIAS
libc_hidden_def (STRCHRNUL)
weak_alias (STRCHRNUL, strchrnul)
#endif
