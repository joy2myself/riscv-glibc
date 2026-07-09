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
L(strchrnul_loop):
	vsetvli	ivl, zero, e8, ELEM_LMUL_SETTING, ta, ma
	vle8ff.v	vstr, (str)
	vmseq.vi	vmask_end, vstr, 0
	vmseq.vx	vmask_ch, vstr, ch
	vfirst.m	end_offset, vmask_end
	vfirst.m	ch_offset, vmask_ch
	sltz	temp1, ch_offset
	sltu	temp2, end_offset, ch_offset
	or	temp1, temp1, temp2
	beqz	temp1, L(found_ch)
	bltz	end_offset, L(advance)
	add	str, str, end_offset
	ret
L(advance):
	csrr	cur_vl, vl
	add	str, str, cur_vl
	j	L(strchrnul_loop)
L(found_ch):
	add	str, str, ch_offset
	ret
.option pop
END (STRCHRNUL)

#ifndef STRCHRNUL_NO_ALIAS
libc_hidden_def (STRCHRNUL)
weak_alias (STRCHRNUL, strchrnul)
#endif
