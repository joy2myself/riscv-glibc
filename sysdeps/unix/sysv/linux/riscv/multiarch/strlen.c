/* Multiple versions of strlen.
   All versions must be listed in ifunc-impl-list.c.
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

#if IS_IN (libc)
/* Redefine strlen so that the compiler won't complain about the type
   mismatch with the IFUNC selector in strong_alias, below.  */
# undef strlen
# define strlen __redirect_strlen
# include <string.h>
# include <profile-ifunc-macros.h>
# include <riscv-ifunc.h>

extern __typeof (__redirect_strlen) __libc_strlen;

extern __typeof (__redirect_strlen) __strlen_generic attribute_hidden;
extern __typeof (__redirect_strlen) __strlen_vector attribute_hidden;

static inline __typeof (__redirect_strlen) *
select_strlen_ifunc (void)
{
  INIT_ARCH ();

  if (RISCV_PROFILE_COND (HAS_VECTOR (), V))
    return __strlen_vector;
  return __strlen_generic;
}

riscv_libc_ifunc (__libc_strlen, select_strlen_ifunc);

# undef strlen
strong_alias (__libc_strlen, strlen);
# ifdef SHARED
__hidden_ver1 (strlen, __GI_strlen, __redirect_strlen)
  __attribute__ ((visibility ("hidden"))) __attribute_copy__ (strlen);
# endif
#else
# include <string/strlen.c>
#endif
