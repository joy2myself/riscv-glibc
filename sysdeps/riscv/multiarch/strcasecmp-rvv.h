/* Helpers for RISC-V RVV strcasecmp family implementations.
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

#include <locale.h>
#include <locale/localeinfo.h>
#include <stdbool.h>

static inline bool
riscv_casecmp_ascii_locale (locale_t loc)
{
  if (loc == LC_GLOBAL_LOCALE)
    loc = _NL_CURRENT_LOCALE;

  return _nl_lookup_word (loc, LC_CTYPE, _NL_CTYPE_NONASCII_CASE) == 0;
}
