/* Generic rawmemchr implementation for RISC-V multiarch.
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
extern __typeof (memchr) __memchr_generic attribute_hidden;
extern __typeof (strlen) __strlen_generic attribute_hidden;

void *
__rawmemchr_generic (const void *s, int c)
{
  if ((unsigned char) c != '\0')
    return __memchr_generic (s, c, (size_t) -1);
  return (char *) s + __strlen_generic (s);
}
#endif
