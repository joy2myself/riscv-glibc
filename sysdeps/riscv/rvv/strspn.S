/* RISC-V RVV based strspn.
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

#ifndef STRSPN
# define STRSPN strspn
#endif

#define result a0

#define str a0
#define accept a1

#define start a2
#define cur a3
#define accept_ptr a4
#define ch a5
#define ivl a6
#define miss_idx a7

#ifndef ELEM_LMUL_SETTING
# define ELEM_LMUL_SETTING m4
#endif
#define vdata v0
#define vmatch v8
#define vtmp v16
#define vmiss v24

ENTRY (STRSPN)
.option push
.option arch, +v
	mv	start, str
	mv	cur, str

L(loop):
	vsetvli	ivl, zero, e8, ELEM_LMUL_SETTING, ta, ma
	vle8ff.v	vdata, (cur)
	csrr	ivl, vl

	vmxor.mm	vmatch, vmatch, vmatch
	mv	accept_ptr, accept

L(accept_loop):
	lbu	ch, 0(accept_ptr)
	beqz	ch, L(find)
	vmseq.vx	vtmp, vdata, ch
	vmor.mm	vmatch, vmatch, vtmp
	addi	accept_ptr, accept_ptr, 1
	j	L(accept_loop)

L(find):
	vmnot.m	vmiss, vmatch
	vfirst.m	miss_idx, vmiss
	bgez	miss_idx, L(found)
	add	cur, cur, ivl
	j	L(loop)

L(found):
	add	result, cur, miss_idx
	sub	result, result, start
	ret
.option pop
END (STRSPN)

#ifndef STRSPN_NO_ALIAS
libc_hidden_builtin_def (strspn)
#endif
