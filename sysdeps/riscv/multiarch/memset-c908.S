/* C908 tuned memset with XThead pair front paths.
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
#define tmp3		t1
#define vdata		v0

#define ELEM_LMUL_SETTING m8

	.attribute unaligned_access, 1
ENTRY (__memset_c908)
	beqz	num, L(ret)

	li	tmp, 3
	bleu	num, tmp, L(copy_1_3_raw)

	li	tmp, 256
	bgtu	num, tmp, L(vector_prepare)

L(scalar_entry):
	li	tmp, 8
	bleu	num, tmp, L(copy_1_8_prep)

	/* Build the repeated-byte word once for the mid/long paths.  */
	andi	fill, value, 0xff
	slli	fill2, fill, 8
	or	fill, fill, fill2
	slli	fill2, fill, 16
	or	fill, fill, fill2
	slli	fill2, fill, 32
	or	fill, fill, fill2
	mv	fill2, fill

	li	tmp, 32
	bleu	num, tmp, L(copy_9_32)

	li	tmp, 64
	bleu	num, tmp, L(copy_33_64_pair)

	li	tmp, 128
	bleu	num, tmp, L(copy_65_128_pair)

	j	L(copy_129_256_pair)

L(vector_prepare):
	mv	dst_ptr, dst

.option push
.option arch,+v
L(vector_fixed_entry):
	vsetvli	ivl, zero, e8, ELEM_LMUL_SETTING, ta, ma
	vmv.v.x	vdata, value
	li	tmp2, 128
L(loop_fixed):
	vse8.v	vdata, (dst_ptr)
	addi	num, num, -128
	addi	dst_ptr, dst_ptr, 128
	bgtu	num, tmp2, L(loop_fixed)
	beqz	num, L(vector_fixed_ret)
	beq	num, tmp2, L(vector_fixed_store)
	vsetvli	ivl, num, e8, ELEM_LMUL_SETTING, ta, ma
L(vector_fixed_store):
	vse8.v	vdata, (dst_ptr)
L(vector_fixed_ret):
	.option pop
L(ret):
	ret

L(copy_1_8):
	li	tmp, 4
	bltu	num, tmp, L(copy_1_3)

L(copy_4_8):
	sw	fill, 0(dst)
	add	tmp2, dst, num
	sw	fill, -4(tmp2)
	ret

L(copy_1_3):
	li	tmp, 1
	beq	num, tmp, L(copy_1)
	li	tmp, 2
	beq	num, tmp, L(copy_2)
	sh	fill, 0(dst)
	sb	fill, 2(dst)
	ret

L(copy_1_8_prep):
	andi	fill, value, 0xff
	slli	fill2, fill, 8
	or	fill, fill, fill2
	slli	fill2, fill, 16
	or	fill, fill, fill2
	j	L(copy_1_8)

L(copy_2):
	sh	fill, 0(dst)
	ret

L(copy_1):
	sb	fill, 0(dst)
	ret

L(copy_9_32):
	li	tmp, 16
	bleu	num, tmp, L(copy_9_16)

	li	tmp, 32
	beq	num, tmp, L(copy_32_pair)

L(copy_17_32):
	sd	fill, 0(dst)
	sd	fill, 8(dst)
	add	tmp2, dst, num
	sd	fill, -16(tmp2)
	sd	fill, -8(tmp2)
	ret

L(copy_32_pair):
.option push
.option arch,+xtheadmempair
	th.sdd	fill, fill2, (dst), 0, 4
	th.sdd	fill, fill2, (dst), 1, 4
.option pop
	ret

L(copy_33_64_pair):
	add	tmp2, dst, num
	addi	tmp2, tmp2, -32
.option push
.option arch,+xtheadmempair
	th.sdd	fill, fill2, (dst), 0, 4
	th.sdd	fill, fill2, (dst), 1, 4
	th.sdd	fill, fill2, (tmp2), 0, 4
	th.sdd	fill, fill2, (tmp2), 1, 4
.option pop
	ret

L(copy_65_128_pair):
	add	tmp2, dst, num
	addi	tmp2, tmp2, -64
.option push
.option arch,+xtheadmempair
	th.sdd	fill, fill2, (dst), 0, 4
	th.sdd	fill, fill2, (dst), 1, 4
	th.sdd	fill, fill2, (dst), 2, 4
	th.sdd	fill, fill2, (dst), 3, 4
	th.sdd	fill, fill2, (tmp2), 0, 4
	th.sdd	fill, fill2, (tmp2), 1, 4
	th.sdd	fill, fill2, (tmp2), 2, 4
	th.sdd	fill, fill2, (tmp2), 3, 4
.option pop
	ret

L(copy_129_256_pair):
	add	tmp2, dst, num
	addi	tmp2, tmp2, -128
	addi	tmp3, dst, 64
.option push
.option arch,+xtheadmempair
	th.sdd	fill, fill2, (dst), 0, 4
	th.sdd	fill, fill2, (dst), 1, 4
	th.sdd	fill, fill2, (dst), 2, 4
	th.sdd	fill, fill2, (dst), 3, 4
	th.sdd	fill, fill2, (tmp3), 0, 4
	th.sdd	fill, fill2, (tmp3), 1, 4
	th.sdd	fill, fill2, (tmp3), 2, 4
	th.sdd	fill, fill2, (tmp3), 3, 4
	addi	tmp3, tmp2, 64
	th.sdd	fill, fill2, (tmp2), 0, 4
	th.sdd	fill, fill2, (tmp2), 1, 4
	th.sdd	fill, fill2, (tmp2), 2, 4
	th.sdd	fill, fill2, (tmp2), 3, 4
	th.sdd	fill, fill2, (tmp3), 0, 4
	th.sdd	fill, fill2, (tmp3), 1, 4
	th.sdd	fill, fill2, (tmp3), 2, 4
	th.sdd	fill, fill2, (tmp3), 3, 4
.option pop
	ret

L(copy_9_16):
	sd	fill, 0(dst)
	add	tmp2, dst, num
	sd	fill, -8(tmp2)
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
	ret
END (__memset_c908)
#endif
