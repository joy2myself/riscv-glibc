/* C908 tuned strcpy.
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
#include <profile-level.h>

#if PROFILE_SHOULD_BUILD (RVA23)
# define pDst		a0
# define pSrc		a1
# define pDstPtr		a2
# define iVL		a3
# define iVlMask		a4
# define iActPos		a5
# define iSrcAlign	a6
# define vMask		v0
# define vStr		v8

ENTRY (__strcpy_c908)
.option push
.option arch, +v,+zbb
	mv	pDstPtr, pDst

	vsetvli	iVL, zero, e8, m1, ta, ma
	vle8ff.v vStr, (pSrc)

	vmseq.vx vMask, vStr, zero
	vfirst.m iActPos, vMask
	bgez	iActPos, L(m1_exit)

	vse8.v	vStr, (pDstPtr)
	csrr	iVL, vl
	add	pSrc, pSrc, iVL
	add	pDstPtr, pDstPtr, iVL

	vsetvli	iVL, zero, e8, m4, ta, ma
	vle8ff.v vStr, (pSrc)

	vmseq.vx vMask, vStr, zero
	vfirst.m iActPos, vMask
	bgez	iActPos, L(m4_exit)

	addi	iVlMask, iVL, -1
	vse8.v	vStr, (pDstPtr)
	vsetvli	zero, zero, e8, m4, ta, ma
	and	iSrcAlign, pSrc, iVlMask
	add	pSrc, pSrc, iVL
	add	pDstPtr, pDstPtr, iVL
	andn	pSrc, pSrc, iVlMask
	sub	pDstPtr, pDstPtr, iSrcAlign

	.p2align 3
L(m4_loop):
	vle8.v	vStr, (pSrc)

	vmseq.vx vMask, vStr, zero
	vfirst.m iActPos, vMask
	bltz	iActPos, L(m4_store)

L(m4_exit):
	addi	iActPos, iActPos, 1
	vsetvli	x0, iActPos, e8, m4, ta, ma
	vse8.v	vStr, (pDstPtr)
	ret

L(m4_store):
	vse8.v	vStr, (pDstPtr)
	add	pSrc, pSrc, iVL
	add	pDstPtr, pDstPtr, iVL
	j	L(m4_loop)

L(m1_exit):
	addi	iActPos, iActPos, 1
	vsetvli	x0, iActPos, e8, m1, ta, ma
	vse8.v	vStr, (pDstPtr)
	ret
.option pop
END (__strcpy_c908)

#endif
