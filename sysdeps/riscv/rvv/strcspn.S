/* RISC-V RVV based strcspn.
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

#ifndef STRCSPN
# define STRCSPN strcspn
#endif

#define result a0

#define str a0
#define reject a1

#define start a2
#define cur a3
#define reject_ptr a4
#define ch a5
#define ivl a6
#define hit_idx a7

#ifndef ELEM_LMUL_SETTING
# define ELEM_LMUL_SETTING m4
#endif
#define vdata v0
#define vhit v8
#define vtmp v16

ENTRY (STRCSPN)
.option push
.option arch, +v
	mv	start, str
	mv	cur, str

L(loop):
	vsetvli	ivl, zero, e8, ELEM_LMUL_SETTING, ta, ma
	vle8ff.v	vdata, (cur)
	csrr	ivl, vl

	vmseq.vi	vhit, vdata, 0
	mv	reject_ptr, reject

L(reject_loop):
	lbu	ch, 0(reject_ptr)
	beqz	ch, L(find)
	vmseq.vx	vtmp, vdata, ch
	vmor.mm	vhit, vhit, vtmp
	addi	reject_ptr, reject_ptr, 1
	j	L(reject_loop)

L(find):
	vfirst.m	hit_idx, vhit
	bgez	hit_idx, L(found)
	add	cur, cur, ivl
	j	L(loop)

L(found):
	add	result, cur, hit_idx
	sub	result, result, start
	ret
.option pop
END (STRCSPN)

#ifndef STRCSPN_NO_ALIAS
libc_hidden_builtin_def (strcspn)
#endif
