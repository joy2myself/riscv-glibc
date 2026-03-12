/* glibc-hwcaps subdirectory test.  RISC-V version.
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

#include <stdio.h>
#include <support/check.h>
#include <sys/param.h>
#include <elf.h>
#include <get-profile-level.h>

extern int marker2 (void);
extern int marker3 (void);
extern int marker4 (void);

/* Return the RISC-V profile level, 1 for the baseline.  */
static int
compute_level (void)
{
  unsigned int profile_level = get_profile_level ();

  if (profile_level < 20)
    return 1;

  /* v2: RVA20U64.  */
  if (profile_level < 22)
    return 2;

  /* v3: RVA22U64.  */
  if (profile_level < 23)
    return 3;

  /* v4: RVA23U64.  */
  return 4;
}

static int
do_test (void)
{
  int level = compute_level ();
  printf ("info: detected RISC-V micro-architecture level: %d\n", level);
  TEST_COMPARE (marker2 (), MIN (level, 2));
  TEST_COMPARE (marker3 (), MIN (level, 3));
  TEST_COMPARE (marker4 (), MIN (level, 4));
  return 0;
}

#include <support/test-driver.c>
