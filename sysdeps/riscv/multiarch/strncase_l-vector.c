/* RISC-V RVV based strncasecmp_l implementation.
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

#include <string.h>

#if IS_IN (libc)
# include <sysdeps/riscv/multiarch/strcasecmp-rvv.h>

extern int __strncasecmp_ascii_vector (const char *, const char *, size_t)
  attribute_hidden;
extern int __strncasecmp_l_generic (const char *, const char *, size_t,
				    locale_t) attribute_hidden;

int __strncasecmp_l_vector (const char *, const char *, size_t, locale_t)
  attribute_hidden;

int
__strncasecmp_l_vector (const char *s1, const char *s2, size_t n,
			locale_t loc)
{
  if (riscv_casecmp_ascii_locale (loc))
    return __strncasecmp_ascii_vector (s1, s2, n);
  return __strncasecmp_l_generic (s1, s2, n, loc);
}
#endif
