/* RISC-V RVV based mempcpy.
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

#ifndef MEMPCPY
# define MEMPCPY __mempcpy
#endif

#define dst a0
#define src a1
#define num a2

#define ivl a3

#ifndef ELEM_LMUL_SETTING
# define ELEM_LMUL_SETTING m4
#endif
#define vdata v0

ENTRY (MEMPCPY)
.option push
.option arch, +v

L(loop):
	vsetvli	ivl, num, e8, ELEM_LMUL_SETTING, ta, ma
	vle8.v	vdata, (src)
	sub	num, num, ivl
	add	src, src, ivl
	vse8.v	vdata, (dst)
	add	dst, dst, ivl
	bnez	num, L(loop)

	ret
.option pop
END (MEMPCPY)

#ifndef MEMPCPY_NO_ALIAS
libc_hidden_def (__mempcpy)
weak_alias (__mempcpy, mempcpy)
libc_hidden_builtin_def (mempcpy)
#endif
