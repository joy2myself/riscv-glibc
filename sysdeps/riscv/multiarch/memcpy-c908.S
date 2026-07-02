/* C908 tuned memcpy using scalar/pair small copies and size-selected RVV.
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
#define src		a1
#define num		a2
#define ivl		a3
#define dst_ptr		a4
#define tmp		a5
#define saved_vtype	a6
#define byte_tmp	a7
#define tmp2		t0
#define tmp3		t1
#define vdata		v0

	.attribute unaligned_access, 1
ENTRY (__memcpy_c908)
.option push
.option arch,+v,+xtheadmempair
	beqz	num, L(ret)
	mv	dst_ptr, dst

	li	tmp, 32
	bleu	num, tmp, L(copy_1_32)

	li	tmp, 64
	bleu	num, tmp, L(copy_33_64_pair)

	li	tmp, 128
	bleu	num, tmp, L(copy_65_128_m8)

	li	tmp, 512
	bltu	num, tmp, L(setup_m4)

	/* Large copies use m8, but first align the store stream to 16 bytes.  */
	andi	tmp, dst_ptr, 15
	beqz	tmp, L(setup_m8)
	li	byte_tmp, 16
	sub	tmp, byte_tmp, tmp

L(align_loop):
	lbu	byte_tmp, 0(src)
	sb	byte_tmp, 0(dst_ptr)
	addi	src, src, 1
	addi	dst_ptr, dst_ptr, 1
	addi	num, num, -1
	addi	tmp, tmp, -1
	bnez	tmp, L(align_loop)

L(setup_m8):
	vsetvli	ivl, zero, e8, m8, ta, ma
	csrr	saved_vtype, vtype
	j	L(copy_loop)

L(setup_m4):
	vsetvli	ivl, zero, e8, m4, ta, ma
	csrr	saved_vtype, vtype

L(copy_loop):
	vsetvl	ivl, num, saved_vtype
	vle8.v	vdata, (src)
	vse8.v	vdata, (dst_ptr)
	sub	num, num, ivl
	add	src, src, ivl
	add	dst_ptr, dst_ptr, ivl
	bnez	num, L(copy_loop)
L(ret):
	ret

L(copy_65_128_m8):
	vsetvli	ivl, num, e8, m8, ta, ma
	vle8.v	vdata, (src)
	vse8.v	vdata, (dst)
	ret

L(copy_33_64_pair):
	add	ivl, src, num
	add	dst_ptr, dst, num
	addi	ivl, ivl, -32
	addi	dst_ptr, dst_ptr, -32
	th.ldd	tmp, saved_vtype, (src), 0, 4
	th.ldd	byte_tmp, tmp2, (src), 1, 4
	th.sdd	tmp, saved_vtype, (dst), 0, 4
	th.sdd	byte_tmp, tmp2, (dst), 1, 4
	th.ldd	tmp, saved_vtype, (ivl), 0, 4
	th.ldd	byte_tmp, tmp2, (ivl), 1, 4
	th.sdd	tmp, saved_vtype, (dst_ptr), 0, 4
	th.sdd	byte_tmp, tmp2, (dst_ptr), 1, 4
	ret

L(copy_1_32):
	li	tmp, 8
	bleu	num, tmp, L(copy_1_8)
	li	tmp, 16
	bleu	num, tmp, L(copy_9_16)

	ld	tmp, 0(src)
	ld	ivl, 8(src)
	sd	tmp, 0(dst)
	sd	ivl, 8(dst)
	add	ivl, src, num
	add	dst_ptr, dst, num
	ld	tmp, -16(ivl)
	ld	ivl, -8(ivl)
	sd	tmp, -16(dst_ptr)
	sd	ivl, -8(dst_ptr)
	ret

L(copy_9_16):
	ld	tmp, 0(src)
	add	ivl, src, num
	add	dst_ptr, dst, num
	ld	ivl, -8(ivl)
	sd	tmp, 0(dst)
	sd	ivl, -8(dst_ptr)
	ret

L(copy_1_8):
	li	tmp, 4
	bltu	num, tmp, L(copy_1_3)
	lw	tmp, 0(src)
	add	ivl, src, num
	add	dst_ptr, dst, num
	lw	ivl, -4(ivl)
	sw	tmp, 0(dst)
	sw	ivl, -4(dst_ptr)
	ret

L(copy_1_3):
	li	tmp, 1
	beq	num, tmp, L(copy_1)
	li	tmp, 2
	beq	num, tmp, L(copy_2)
	lhu	tmp, 0(src)
	lbu	ivl, 2(src)
	sh	tmp, 0(dst)
	sb	ivl, 2(dst)
	ret

L(copy_2):
	lhu	tmp, 0(src)
	sh	tmp, 0(dst)
	ret

L(copy_1):
	lbu	tmp, 0(src)
	sb	tmp, 0(dst)
	ret
.option pop
END (__memcpy_c908)
#endif
