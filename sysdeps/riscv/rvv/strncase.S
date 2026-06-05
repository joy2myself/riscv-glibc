/* RISC-V RVV based ASCII strncasecmp core.
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

#ifndef STRNCASECMP
# define STRNCASECMP __strncasecmp
# define STRNCASECMP_CHECK_LOCALE 1
# define STRNCASECMP_PUBLIC 1
#endif

#ifdef STRNCASECMP_CHECK_LOCALE
#include <locale-defines.h>
#endif

#ifndef STRNCASECMP_GENERIC
# define STRNCASECMP_GENERIC __strncasecmp_generic
#endif

#define result a0

#define str1 a0
#define str2 a1
#define num a2

#define p1 a3
#define p2 a4
#define ivl a5
#define idx a6
#define lower_bound a0
#define upper_bound a1
#define case_add a7

#ifndef ELEM_LMUL_SETTING
# define ELEM_LMUL_SETTING m4
#endif
#define vstr1 v8
#define vstr2 v16
#define vtmp v24
#define vmask v0

ENTRY (STRNCASECMP)
.option push
.option arch, +v
	beqz	num, L(equal)
	beq	str1, str2, L(equal)

#ifdef STRNCASECMP_CHECK_LOCALE
	la.tls.ie	t4, __libc_tsd_LOCALE
	add	t4, t4, tp
	REG_L	t4, 0(t4)
	li	t5, LOCALE_T___LOCALES+LC_CTYPE*SZREG
	add	t4, t4, t5
	REG_L	t4, 0(t4)
	li	t5, LOCALE_DATA_VALUES+_NL_CTYPE_NONASCII_CASE*SIZEOF_VALUES
	add	t4, t4, t5
	lw	t4, 0(t4)
	bnez	t4, L(nonascii)
#endif
	mv	p1, str1
	mv	p2, str2
	li	lower_bound, 'A'
	li	upper_bound, 'Z' + 1
	li	case_add, 'a' - 'A'

L(loop):
	vsetvli	ivl, num, e8, ELEM_LMUL_SETTING, ta, mu
	vle8ff.v	vstr1, (p1)
	vle8ff.v	vstr2, (p2)
	csrr	ivl, vl

	vmsgeu.vx	vmask, vstr1, lower_bound
	vmsltu.vx	vtmp, vstr1, upper_bound
	vmand.mm	vmask, vmask, vtmp
	vadd.vx	vstr1, vstr1, case_add, vmask.t

	vmsgeu.vx	vmask, vstr2, lower_bound
	vmsltu.vx	vtmp, vstr2, upper_bound
	vmand.mm	vmask, vmask, vtmp
	vadd.vx	vstr2, vstr2, case_add, vmask.t

	vmseq.vi	vmask, vstr1, 0
	vmsne.vv	vtmp, vstr1, vstr2
	vmor.mm	vmask, vmask, vtmp
	vfirst.m	idx, vmask
	bgez	idx, L(found)

	sub	num, num, ivl
	beqz	num, L(equal)
	add	p1, p1, ivl
	add	p2, p2, ivl
	j	L(loop)

L(found):
	vslidedown.vx	vtmp, vstr1, idx
	vmv.x.s	lower_bound, vtmp
	andi	lower_bound, lower_bound, 255
	vslidedown.vx	vtmp, vstr2, idx
	vmv.x.s	upper_bound, vtmp
	andi	upper_bound, upper_bound, 255
	sub	result, lower_bound, upper_bound
	ret

L(equal):
	li	result, 0
	ret

#ifdef STRNCASECMP_CHECK_LOCALE
L(nonascii):
	tail	STRNCASECMP_GENERIC
#endif
.option pop
END (STRNCASECMP)

libc_hidden_def (STRNCASECMP)
#ifdef STRNCASECMP_PUBLIC
weak_alias (STRNCASECMP, strncasecmp)
#endif
