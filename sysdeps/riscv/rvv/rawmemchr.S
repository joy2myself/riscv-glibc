/* RISC-V RVV based rawmemchr.
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

#ifndef RAWMEMCHR
# define RAWMEMCHR __rawmemchr
#endif

#define result a0

#define src a0
#define value a1

#define ivl a3
#define temp a4

#define ELEM_LMUL_SETTING m8
#define vdata v0
#define vmask v8

ENTRY (RAWMEMCHR)
.option push
.option arch, +v

L(loop):
	vsetvli	ivl, zero, e8, ELEM_LMUL_SETTING, ta, ma
	vle8ff.v	vdata, (src)
	vmseq.vx	vmask, vdata, value
	vfirst.m	temp, vmask
	bgez	temp, L(found)
	csrr	ivl, vl
	add	src, src, ivl
	j	L(loop)

L(found):
	add	result, src, temp
	ret
.option pop
END (RAWMEMCHR)

#ifndef RAWMEMCHR_NO_ALIAS
libc_hidden_def (__rawmemchr)
weak_alias (__rawmemchr, rawmemchr)
#endif
