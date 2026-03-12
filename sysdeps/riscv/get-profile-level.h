/* Get RISC-V profile level.
   This file is part of the GNU C Library.
   Copyright (C) 2026 Free Software Foundation, Inc.

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

#include <profile-detection.h>

static unsigned int
get_profile_level (void)
{
  struct riscv_hwprobe pair;
  uint64_t cpuperf;
  uint64_t extensions;
  unsigned int zicbop_block_size = 0;
  unsigned int zicboz_block_size = 0;

  pair.key = RISCV_HWPROBE_KEY_CPUPERF_0;
  if (__riscv_hwprobe (&pair, 1, 0, NULL, 0) != 0)
    return 0;
  cpuperf = pair.value;

  pair.key = RISCV_HWPROBE_KEY_IMA_EXT_0;
  if (__riscv_hwprobe (&pair, 1, 0, NULL, 0) != 0)
    return 0;
  extensions = pair.value;

  if (extensions & RISCV_HWPROBE_EXT_ZICBOP)
    zicbop_block_size =
      riscv_hwprobe_get_block_size (RISCV_HWPROBE_KEY_ZICBOP_BLOCK_SIZE);

  if (extensions & RISCV_HWPROBE_EXT_ZICBOZ)
    zicboz_block_size =
      riscv_hwprobe_get_block_size (RISCV_HWPROBE_KEY_ZICBOZ_BLOCK_SIZE);

  return riscv_get_profile_level_from_hwprobe (extensions, cpuperf,
					       zicbop_block_size,
					       zicboz_block_size);
}
