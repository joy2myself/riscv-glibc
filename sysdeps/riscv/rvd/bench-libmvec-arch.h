/* Runtime architecture check for libmvec benchtests.  RISC-V version.
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

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/hwprobe.h>

#define INIT_ARCH()                                                           \
  do                                                                          \
  {                                                                           \
    if (!supported ())                                                        \
      return 77;                                                              \
  }                                                                           \
  while (0)

static bool
supported (void)
{
  struct riscv_hwprobe pair = { .key = RISCV_HWPROBE_KEY_IMA_EXT_0 };
  if (__riscv_hwprobe (&pair, 1, 0, NULL, 0) != 0)
  {
    printf ("__riscv_hwprobe failed.\n");
    return false;
  }

  if ((pair.value & RISCV_HWPROBE_IMA_V) == 0)
  {
    printf ("RVV not supported.\n");
    return false;
  }

  return true;
}
