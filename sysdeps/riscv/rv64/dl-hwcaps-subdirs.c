/* Architecture-specific glibc-hwcaps subdirectories.  RISC-V version.
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

#include <dl-hwcaps.h>
#include <ldsodefs.h>
#include <cpu-features.h>

const char _dl_hwcaps_subdirs[] = "rva23u64:rva22u64:rva20u64";
enum { subdirs_count = 3 }; /* Number of components in _dl_hwcaps_subdirs.  */

uint32_t
_dl_hwcaps_subdirs_active (void)
{
  int active = 0;

  /* Test in reverse preference order.  */
  /* v1: RVA20U64.  */
  unsigned int profile_level = GLRO(dl_riscv_cpu_features).profile_level;
  if (profile_level < 20)
    return _dl_hwcaps_subdirs_build_bitmask (subdirs_count, active);
  ++active;

  /* v2: RVA22U64.  */
  if (profile_level < 22)
    return _dl_hwcaps_subdirs_build_bitmask (subdirs_count, active);
  ++active;

  /* v3: RVA23U64.  */
  if (profile_level < 23)
    return _dl_hwcaps_subdirs_build_bitmask (subdirs_count, active);
  ++active;
  return _dl_hwcaps_subdirs_build_bitmask (subdirs_count, active);
}
