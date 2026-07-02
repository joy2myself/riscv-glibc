/* SpaceMiT X60 tuned memmove implementation with size-selected RVV dispatch.
   Copyright (C) 2026 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public License as
   published by the Free Software Foundation; either version 2.1 of the
   License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with the GNU C Library; if not, see <https://www.gnu.org/licenses/>.
   */

#include <profile-level.h>

#if PROFILE_SHOULD_BUILD (RVA23)
# include <sysdep.h>
# include <sys/asm.h>

# define MEMMOVE __memmove_spacemit_x60
# define dst a0
# define src a1
# define num a2
# define delta t0
# define ivl a3
# define dst_ptr a4
# define src_end a5
# define dst_end a6
# define tmp0 a7
# define tmp1 t1
# define tmp2 t2
# define tmp3 t3
# define tmp4 t4
# define tmp5 t5
# define tmp6 t6
# define saved_vtype t4
# define vdata v0

	.attribute unaligned_access, 1
ENTRY (MEMMOVE)
.option push
.option arch, +v
	beqz	num, L(ret)
	sub	delta, dst, src
	beqz	delta, L(ret)

	/* Small copies are overlap-safe because all source words are read first.  */
	li	tmp0, 64
	bleu	num, tmp0, L(copy_1_64)

	/* Choose the vector regime once, then share one delta dispatch.  */
	li	tmp0, 128
	bltu	num, tmp0, L(setup_m1)
	li	tmp0, 512
	bltu	num, tmp0, L(setup_m4)

	mv	dst_ptr, dst
	andi	tmp0, dst_ptr, 15
	beqz	tmp0, L(setup_m8)

	/* Only non-overlap/forward copies may consume an aligned scalar prologue.  */
	bltu	delta, num, L(setup_m8)
	li	tmp1, 16
	sub	tmp0, tmp1, tmp0
L(align_forward_loop):
	lbu	tmp1, 0(src)
	sb	tmp1, 0(dst_ptr)
	addi	src, src, 1
	addi	dst_ptr, dst_ptr, 1
	addi	num, num, -1
	addi	tmp0, tmp0, -1
	bnez	tmp0, L(align_forward_loop)
	j	L(setup_m8_ready)

L(setup_m8):
	mv	dst_ptr, dst
L(setup_m8_ready):
	vsetvli	ivl, zero, e8, m8, ta, ma
	csrr	saved_vtype, vtype
	j	L(dispatch)

L(setup_m4):
	mv	dst_ptr, dst
	vsetvli	ivl, zero, e8, m4, ta, ma
	csrr	saved_vtype, vtype
	j	L(dispatch)

L(setup_m1):
	mv	dst_ptr, dst
	vsetvli	ivl, zero, e8, m1, ta, ma
	csrr	saved_vtype, vtype

L(dispatch):
	bltu	delta, num, L(backward_loop_prepare)

L(forward_loop):
	vsetvl	ivl, num, saved_vtype
	vle8.v	vdata, (src)
	vse8.v	vdata, (dst_ptr)
	sub	num, num, ivl
	add	src, src, ivl
	add	dst_ptr, dst_ptr, ivl
	bnez	num, L(forward_loop)
	ret

L(backward_loop_prepare):
	add	src_end, src, num
	add	dst_end, dst, num

L(backward_loop):
	vsetvl	ivl, num, saved_vtype
	sub	src_end, src_end, ivl
	vle8.v	vdata, (src_end)
	sub	dst_end, dst_end, ivl
	vse8.v	vdata, (dst_end)
	sub	num, num, ivl
	bnez	num, L(backward_loop)
	ret

L(copy_1_64):
	li	tmp0, 8
	bleu	num, tmp0, L(copy_1_8)
	li	tmp0, 16
	bleu	num, tmp0, L(copy_9_16)
	li	tmp0, 32
	bleu	num, tmp0, L(copy_17_32)

	add	tmp5, src, num
	add	tmp6, dst, num
	addi	tmp5, tmp5, -32
	addi	tmp6, tmp6, -32
	ld	ivl, 0(src)
	ld	src_end, 8(src)
	ld	dst_end, 16(src)
	ld	tmp0, 24(src)
	ld	tmp1, 0(tmp5)
	ld	tmp2, 8(tmp5)
	ld	tmp3, 16(tmp5)
	ld	tmp4, 24(tmp5)
	sd	ivl, 0(dst)
	sd	src_end, 8(dst)
	sd	dst_end, 16(dst)
	sd	tmp0, 24(dst)
	sd	tmp1, 0(tmp6)
	sd	tmp2, 8(tmp6)
	sd	tmp3, 16(tmp6)
	sd	tmp4, 24(tmp6)
	ret

L(copy_17_32):
	ld	ivl, 0(src)
	ld	src_end, 8(src)
	add	tmp5, src, num
	add	tmp6, dst, num
	ld	dst_end, -16(tmp5)
	ld	tmp0, -8(tmp5)
	sd	ivl, 0(dst)
	sd	src_end, 8(dst)
	sd	dst_end, -16(tmp6)
	sd	tmp0, -8(tmp6)
	ret

L(copy_9_16):
	ld	ivl, 0(src)
	add	tmp5, src, num
	add	tmp6, dst, num
	ld	tmp0, -8(tmp5)
	sd	ivl, 0(dst)
	sd	tmp0, -8(tmp6)
	ret

L(copy_1_8):
	li	tmp0, 4
	bltu	num, tmp0, L(copy_1_3)
	lw	ivl, 0(src)
	add	tmp5, src, num
	add	tmp6, dst, num
	lw	tmp0, -4(tmp5)
	sw	ivl, 0(dst)
	sw	tmp0, -4(tmp6)
	ret

L(copy_1_3):
	li	tmp0, 1
	beq	num, tmp0, L(copy_1)
	li	tmp0, 2
	beq	num, tmp0, L(copy_2)
	lhu	ivl, 0(src)
	lbu	tmp0, 2(src)
	sh	ivl, 0(dst)
	sb	tmp0, 2(dst)
	ret

L(copy_2):
	lhu	ivl, 0(src)
	sh	ivl, 0(dst)
	ret

L(copy_1):
	lbu	ivl, 0(src)
	sb	ivl, 0(dst)
L(ret):
	ret
.option pop
END (MEMMOVE)
#endif
