/* Multiple versions of strchrnul.
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
/* Redefine strchrnul so that the compiler won't complain about the type
   mismatch with the IFUNC selector in weak_alias, below.  */
# define strchrnul __redirect_strchrnul
# define __strchrnul __redirect___strchrnul
# include <stdint.h>
# include <string.h>
# undef strchrnul
# undef __strchrnul
# include <ifunc-init.h>
# include <riscv-ifunc.h>
# include <sys/hwprobe.h>

extern __typeof (__redirect_strchrnul) __strchrnul;
extern __typeof (__redirect_strchrnul) __strchrnul_generic attribute_hidden;
extern __typeof (__redirect_strchrnul) __strchrnul_vector attribute_hidden;

static inline __typeof (__redirect_strchrnul) *
select_strchrnul_ifunc (uint64_t dl_hwcap, __riscv_hwprobe_t hwprobe_func)
{
  if (riscv_hwprobe_has_vector (hwprobe_func))
    return __strchrnul_vector;
  return __strchrnul_generic;
}

riscv_libc_ifunc (__strchrnul, select_strchrnul_ifunc);
weak_alias (__strchrnul, strchrnul)
# ifdef SHARED
__hidden_ver1 (__strchrnul, __GI___strchrnul, __redirect___strchrnul)
  __attribute__ ((visibility ("hidden"))) __attribute_copy__ (__strchrnul);
# endif
#else
# include <string/strchrnul.c>
#endif
