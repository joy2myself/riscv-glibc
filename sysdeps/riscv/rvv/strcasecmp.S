/* RISC-V RVV based ASCII strcasecmp core.
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

#ifndef STRCASECMP
# define STRCASECMP __strcasecmp
# define STRCASECMP_CHECK_LOCALE 1
# define STRCASECMP_PUBLIC 1
#endif

#ifdef STRCASECMP_CHECK_LOCALE
#include <locale-defines.h>
#endif

#ifndef STRCASECMP_GENERIC
# define STRCASECMP_GENERIC __strcasecmp_generic
#endif

#define result a0

#define str1 a0
#define str2 a1

#define p1 a2
#define p2 a3
#define ivl a4
#define idx a5
#define lower_bound a6
#define upper_bound a7
#define case_add a0

#ifndef ELEM_LMUL_SETTING
# define ELEM_LMUL_SETTING m4
#endif
#define vstr1 v8
#define vstr2 v16
#define vtmp v24
#define vmask v0

ENTRY (STRCASECMP)
.option push
.option arch, +v
	beq	str1, str2, L(equal)

#ifdef STRCASECMP_CHECK_LOCALE
	la.tls.ie	t3, __libc_tsd_LOCALE
	add	t3, t3, tp
	REG_L	t3, 0(t3)
	li	t4, LOCALE_T___LOCALES+LC_CTYPE*SZREG
	add	t3, t3, t4
	REG_L	t3, 0(t3)
	li	t4, LOCALE_DATA_VALUES+_NL_CTYPE_NONASCII_CASE*SIZEOF_VALUES
	add	t3, t3, t4
	lw	t3, 0(t3)
	bnez	t3, L(nonascii)
#endif
	mv	p1, str1
	mv	p2, str2
	li	lower_bound, 'A'
	li	upper_bound, 'Z' + 1
	li	case_add, 'a' - 'A'

L(loop):
	vsetvli	ivl, zero, e8, ELEM_LMUL_SETTING, ta, mu
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

#ifdef STRCASECMP_CHECK_LOCALE
L(nonascii):
	tail	STRCASECMP_GENERIC
#endif
.option pop
END (STRCASECMP)

libc_hidden_def (STRCASECMP)
#ifdef STRCASECMP_PUBLIC
weak_alias (STRCASECMP, strcasecmp)
#endif
