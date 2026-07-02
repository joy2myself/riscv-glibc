/* RISC-V RVV based strstr.
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

#ifndef STRSTR
# define STRSTR strstr
#endif

#define result a0
#define haystack a0
#define needle a1

#define haystack_base a2
#define needle_base a3
#define scan a4
#define needle_len a5
#define haystack_len a6
#define ivl a7

#define offset t0
#define end_scan t1
#define candidate t2
#define cmp_haystack t3
#define cmp_needle t4
#define remaining t5
#define first_ch t6
#define last_ch a1

#define ELEM_LMUL_SETTING m8
#define vhaystack v0
#define vneedle v8
#define vmask v16
#define vmask_last v24

ENTRY (STRSTR)
.option push
.option arch, +v
    mv haystack_base, haystack
    mv needle_base, needle

    mv scan, needle_base
L(strlen_needle_loop):
    vsetvli ivl, zero, e8, ELEM_LMUL_SETTING, ta, ma
    vle8ff.v vhaystack, (scan)
    csrr ivl, vl
    vmseq.vi vmask, vhaystack, 0
    vfirst.m offset, vmask
    bgez offset, L(strlen_needle_found)
    add scan, scan, ivl
    j L(strlen_needle_loop)

L(strlen_needle_found):
    add scan, scan, offset
    sub needle_len, scan, needle_base
    beqz needle_len, L(return_haystack)

    mv scan, haystack_base
L(strlen_haystack_loop):
    vsetvli ivl, zero, e8, ELEM_LMUL_SETTING, ta, ma
    vle8ff.v vhaystack, (scan)
    csrr ivl, vl
    vmseq.vi vmask, vhaystack, 0
    vfirst.m offset, vmask
    bgez offset, L(strlen_haystack_found)
    add scan, scan, ivl
    j L(strlen_haystack_loop)

L(strlen_haystack_found):
    add scan, scan, offset
    sub haystack_len, scan, haystack_base
    bltu haystack_len, needle_len, L(not_found)

    lbu first_ch, 0(needle_base)
    add remaining, needle_base, needle_len
    addi remaining, remaining, -1
    lbu last_ch, 0(remaining)

    add end_scan, haystack_base, haystack_len
    sub end_scan, end_scan, needle_len
    addi end_scan, end_scan, 1
    mv scan, haystack_base

L(scan_candidates):
    sub remaining, end_scan, scan
    blez remaining, L(not_found)

    vsetvli ivl, remaining, e8, ELEM_LMUL_SETTING, ta, ma
    vle8.v vhaystack, (scan)
    vmseq.vx vmask, vhaystack, first_ch

    add cmp_haystack, scan, needle_len
    addi cmp_haystack, cmp_haystack, -1
    vle8.v vneedle, (cmp_haystack)
    vmseq.vx vmask_last, vneedle, last_ch
    vmand.mm vmask, vmask, vmask_last

    vfirst.m offset, vmask
    bgez offset, L(verify_candidate)

    add scan, scan, ivl
    j L(scan_candidates)

L(verify_candidate):
    add candidate, scan, offset
    mv cmp_haystack, candidate
    mv cmp_needle, needle_base
    mv remaining, needle_len

L(compare_loop):
    vsetvli ivl, remaining, e8, ELEM_LMUL_SETTING, ta, ma
    vle8.v vhaystack, (cmp_haystack)
    vle8.v vneedle, (cmp_needle)
    vmsne.vv vmask, vhaystack, vneedle
    vfirst.m offset, vmask
    bgez offset, L(candidate_mismatch)

    sub remaining, remaining, ivl
    add cmp_haystack, cmp_haystack, ivl
    add cmp_needle, cmp_needle, ivl
    bnez remaining, L(compare_loop)

    mv result, candidate
    ret

L(candidate_mismatch):
    addi scan, candidate, 1
    j L(scan_candidates)

L(return_haystack):
    mv result, haystack_base
    ret

L(not_found):
    li result, 0
    ret
.option pop
END (STRSTR)

#ifndef STRSTR_NO_ALIAS
libc_hidden_builtin_def (strstr)
#endif
