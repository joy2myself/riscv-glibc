/* RISC-V RVV based strnlen.
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

#ifndef STRNLEN
# define STRNLEN __strnlen
#endif

#define str a0
#define copy_str a2
#define ret_value a0
#define max_len a1
#define cur_vl a3
#define end_offset a4

#define ELEM_LMUL_SETTING m1
#define vstr v0
#define vmask_end v8

ENTRY (STRNLEN)
.option push
.option arch, +v
    mv copy_str, str
    mv ret_value, max_len
L(strnlen_vector_loop):
    beqz max_len, L(end_strnlen_vector_loop)
    vsetvli zero, max_len, e8, ELEM_LMUL_SETTING, ta, ma
    vle8ff.v vstr, (copy_str)
    vmseq.vi vmask_end, vstr, 0
    vfirst.m end_offset, vmask_end /* first occurence of \0 */
    csrr cur_vl, vl
    add copy_str, copy_str, cur_vl
    sub max_len, max_len, cur_vl
    bltz end_offset, L(strnlen_vector_loop)
    add max_len, max_len, cur_vl
    sub ret_value, ret_value, max_len
    add ret_value, ret_value, end_offset
L(end_strnlen_vector_loop):
    ret
.option pop
END (STRNLEN)
weak_alias (STRNLEN, strnlen)
libc_hidden_def (STRNLEN)
libc_hidden_def (strnlen)
