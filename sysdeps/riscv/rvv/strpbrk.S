/* RISC-V RVV based strpbrk.
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

#ifndef STRPBRK
# define STRPBRK strpbrk
#endif

#define result a0

#define str a0
#define accept a1

#define cur a2
#define accept_ptr a3
#define ch a4
#define ivl a5
#define hit_idx a6
#define hit_ptr a7

#ifndef ELEM_LMUL_SETTING
# define ELEM_LMUL_SETTING m4
#endif
#define vdata v0
#define vhit v8
#define vtmp v16

ENTRY (STRPBRK)
.option push
.option arch, +v
	mv	cur, str

L(loop):
	vsetvli	ivl, zero, e8, ELEM_LMUL_SETTING, ta, ma
	vle8ff.v	vdata, (cur)
	csrr	ivl, vl

	vmseq.vi	vhit, vdata, 0
	mv	accept_ptr, accept

L(accept_loop):
	lbu	ch, 0(accept_ptr)
	beqz	ch, L(find)
	vmseq.vx	vtmp, vdata, ch
	vmor.mm	vhit, vhit, vtmp
	addi	accept_ptr, accept_ptr, 1
	j	L(accept_loop)

L(find):
	vfirst.m	hit_idx, vhit
	bgez	hit_idx, L(found)
	add	cur, cur, ivl
	j	L(loop)

L(found):
	add	hit_ptr, cur, hit_idx
	lbu	ch, 0(hit_ptr)
	beqz	ch, L(null)
	mv	result, hit_ptr
	ret

L(null):
	li	result, 0
	ret
.option pop
END (STRPBRK)

#ifndef STRPBRK_NO_ALIAS
libc_hidden_builtin_def (strpbrk)
#endif
