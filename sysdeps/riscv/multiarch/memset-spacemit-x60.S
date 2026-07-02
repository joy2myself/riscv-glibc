/* SpaceMiT X60 tuned memset.
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

#include <profile-level.h>

#if PROFILE_SHOULD_BUILD (RVA23)
#include <sysdep.h>
#include <sys/asm.h>

#define dst		a0
#define value		a1
#define num		a2

#define ivl		a3
#define dst_ptr		a4
#define fill		a5
#define fill2		a6
#define tmp		a7
#define tmp2		t0
#define vdata		v0

	.attribute unaligned_access, 1
ENTRY (__memset_spacemit_x60)
	beqz	num, L(ret)

	li	tmp, 15
	bleu	num, tmp, L(copy_1_15)

	li	tmp, 31
	bleu	num, tmp, L(vector_m2_once)

	li	tmp, 63
	bleu	num, tmp, L(copy_32_63_scalar)

	li	tmp, 255
	bleu	num, tmp, L(vector_m4_loop)

	li	tmp, 2048
	bleu	num, tmp, L(vector_m4_fixed)

	j	L(vector_m8_fixed)

L(copy_1_15):
	li	tmp, 3
	bleu	num, tmp, L(copy_1_3_raw)

.option push
.option arch,+v
L(vector_m1_once):
	vsetvli	ivl, num, e8, m1, ta, ma
	vmv.v.x	vdata, value
	vse8.v	vdata, (dst)
	.option pop
	ret

.option push
.option arch,+v
L(vector_m2_once):
	vsetvli	ivl, num, e8, m2, ta, ma
	vmv.v.x	vdata, value
	vse8.v	vdata, (dst)
	.option pop
	ret

L(vector_m4_loop):
	mv	dst_ptr, dst
.option push
.option arch,+v
	vsetvli	ivl, num, e8, m4, ta, ma
	vmv.v.x	vdata, value
L(loop_m4):
	vse8.v	vdata, (dst_ptr)
	sub	num, num, ivl
	add	dst_ptr, dst_ptr, ivl
	beqz	num, L(vector_ret)
	vsetvli	ivl, num, e8, m4, ta, ma
	j	L(loop_m4)

L(vector_m4_fixed):
	mv	dst_ptr, dst
	vsetvli	ivl, zero, e8, m4, ta, ma
	vmv.v.x	vdata, value
L(loop_m4_fixed):
	vse8.v	vdata, (dst_ptr)
	sub	num, num, ivl
	add	dst_ptr, dst_ptr, ivl
	bgtu	num, ivl, L(loop_m4_fixed)
	beqz	num, L(vector_ret)
	beq	num, ivl, L(vector_m4_store)
	vsetvli	ivl, num, e8, m4, ta, ma
L(vector_m4_store):
	vse8.v	vdata, (dst_ptr)
	j	L(vector_ret)

L(vector_m8_fixed):
	mv	dst_ptr, dst
	vsetvli	ivl, zero, e8, m8, ta, ma
	vmv.v.x	vdata, value
L(loop_m8_fixed):
	vse8.v	vdata, (dst_ptr)
	sub	num, num, ivl
	add	dst_ptr, dst_ptr, ivl
	bgtu	num, ivl, L(loop_m8_fixed)
	beqz	num, L(vector_ret)
	beq	num, ivl, L(vector_m8_store)
	vsetvli	ivl, num, e8, m8, ta, ma
L(vector_m8_store):
	vse8.v	vdata, (dst_ptr)
L(vector_ret):
	.option pop
	ret

L(copy_32_63_scalar):
	andi	fill, value, 0xff
	slli	fill2, fill, 8
	or	fill, fill, fill2
	slli	fill2, fill, 16
	or	fill, fill, fill2
	slli	fill2, fill, 32
	or	fill, fill, fill2

	add	dst_ptr, dst, num
	sd	fill, 0(dst)
	sd	fill, 8(dst)
	sd	fill, 16(dst)
	sd	fill, 24(dst)
	sd	fill, -32(dst_ptr)
	sd	fill, -24(dst_ptr)
	sd	fill, -16(dst_ptr)
	sd	fill, -8(dst_ptr)
	ret

L(copy_1_3_raw):
	li	tmp, 1
	beq	num, tmp, L(copy_1_raw)
	li	tmp, 2
	beq	num, tmp, L(copy_2_raw)
	sb	value, 0(dst)
	sb	value, 1(dst)
	sb	value, 2(dst)
	ret

L(copy_2_raw):
	sb	value, 0(dst)
	sb	value, 1(dst)
	ret

L(copy_1_raw):
	sb	value, 0(dst)
L(ret):
	ret
END (__memset_spacemit_x60)
#endif
