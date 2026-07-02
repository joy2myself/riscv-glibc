/* Multiple versions of strstr.
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
/* Redefine strstr so that the compiler won't complain about the type
   mismatch with the IFUNC selector in strong_alias, below.  */
# undef strstr
# define strstr __redirect_strstr
# include <stdint.h>
# include <string.h>
# undef strstr

extern __typeof (__redirect_strstr) __strstr_generic attribute_hidden;
extern __typeof (__redirect_strstr) __strstr_vector attribute_hidden;
extern __typeof (__redirect_strstr) __libc_strstr;

# include <profile-ifunc-macros.h>
# include <riscv-ifunc.h>

static inline __typeof (__redirect_strstr) *
select_strstr_ifunc (uint64_t dl_hwcap, __riscv_hwprobe_t hwprobe_func)
{
  (void) dl_hwcap;
  (void) hwprobe_func;
  INIT_ARCH ();

  if (RISCV_PROFILE_COND (HAS_VECTOR (), V))
    return __strstr_vector;

  return __strstr_generic;
}

riscv_libc_ifunc (__libc_strstr, select_strstr_ifunc);
# undef strstr
strong_alias (__libc_strstr, strstr)
#else
# include <string/strstr.c>
#endif
