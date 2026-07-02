/* RISC-V RVV based memmove.
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

#ifndef MEMMOVE
# define MEMMOVE memmove
#endif

#define dst a0
#define src a1
#define num a2

#define delta t0
#define ivl a3
#define dst_ptr a4
#define src_backward_ptr a5
#define dst_backward_ptr a6
#define vdata v0

ENTRY (MEMMOVE)
.option push
.option arch, +v
	sub	delta, dst, src
	beqz	delta, L(done)
	mv	dst_ptr, dst
	bltu	delta, num, L(prepare_backward)

L(forward_copy_loop):
	vsetvli	ivl, num, e8, m8, ta, ma
	vle8.v	vdata, (src)
	sub	num, num, ivl
	add	src, src, ivl
	vse8.v	vdata, (dst_ptr)
	add	dst_ptr, dst_ptr, ivl
	bnez	num, L(forward_copy_loop)
	j	L(done)

L(prepare_backward):
	add	src_backward_ptr, src, num
	add	dst_backward_ptr, dst, num

L(backward_copy_loop):
	vsetvli	ivl, num, e8, m8, ta, ma
	sub	src_backward_ptr, src_backward_ptr, ivl
	vle8.v	vdata, (src_backward_ptr)
	sub	num, num, ivl
	sub	dst_backward_ptr, dst_backward_ptr, ivl
	vse8.v	vdata, (dst_backward_ptr)
	bnez	num, L(backward_copy_loop)

L(done):
	ret

.option pop
END (MEMMOVE)
libc_hidden_builtin_def (memmove)
