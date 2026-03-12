/* Symbol redirection for loader/static initialization code.
   Copyright (C) 2025-2026 Free Software Foundation, Inc.
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

#ifndef _DL_IFUNC_GENERIC_H
#define _DL_IFUNC_GENERIC_H

# include <profile-level.h>

# ifndef SHARED
#  if MINIMUM_RISCV_PROFILE_LEVEL >= 23
asm ("memset = __memset_vector");
asm ("memcpy = __memcpy_vector");
asm ("strlen = __strlen_vector");
#  elif MINIMUM_RISCV_PROFILE_LEVEL >= 20
asm ("memset = __memset_generic");
asm ("memcpy = __memcpy_noalignment");
asm ("strlen = __strlen_generic");
#  else
asm ("memset = __memset_generic");
asm ("memcpy = __memcpy_generic");
asm ("strlen = __strlen_generic");
#  endif

# endif

#endif
