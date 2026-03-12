/* Multiple versions of memset.
   All versions must be listed in ifunc-impl-list.c.
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

#if IS_IN (libc)
/* Redefine memset so that the compiler won't complain about the type
   mismatch with the IFUNC selector in strong_alias, below.  */
# undef memset
# define memset __redirect_memset
# include <string.h>
# include <riscv-ifunc.h>
# include <profile-ifunc-macros.h>

extern __typeof (__redirect_memset) __libc_memset;

extern __typeof (__redirect_memset) __memset_generic attribute_hidden;
extern __typeof (__redirect_memset) __memset_vector attribute_hidden;

static inline __typeof (__redirect_memset) *
select_memset_ifunc (void)
{
  INIT_ARCH ();

  if (RISCV_PROFILE_COND (HAS_VECTOR (), V))
    return __memset_vector;

  return __memset_generic;
}

riscv_libc_ifunc (__libc_memset, select_memset_ifunc);

# undef memset
strong_alias (__libc_memset, memset);
# ifdef SHARED
__hidden_ver1 (memset, __GI_memset, __redirect_memset)
  __attribute__ ((visibility ("hidden"))) __attribute_copy__ (memset);
# endif
#else
# include <string/memset.c>
#endif
