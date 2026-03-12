/* RISC-V CPU feature tuning.
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
#include <stdint.h>
#include <unistd.h>
#include <elf/dl-tunables.h>
#include <string.h>
#include <cpu-features.h>
#include <ldsodefs.h>
#include <sys/auxv.h>
#include <dl-tunables-parse.h>

struct riscv_hwcap_tunable
{
  const char *name;
  size_t len;
  uint64_t bit;
};

#define RISCV_HWCAP_ENTRY(name, bit) \
  { #name, sizeof (#name) - 1, (uint64_t) (bit) },
#define RISCV_HWCAP_ALIAS(str, bit) \
  { str, sizeof (str) - 1, (uint64_t) (bit) },

static const struct riscv_hwcap_tunable riscv_hwcaps[] =
  {
#include <cpu-hwcaps-list.h>
  };

static inline void
riscv_apply_tunable_hwcaps (const struct tunable_str_comma_t *t)
{
  for (size_t i = 0; i < sizeof (riscv_hwcaps) / sizeof (riscv_hwcaps[0]);
       ++i)
    if (tunable_str_comma_strcmp (t, riscv_hwcaps[i].name,
				  riscv_hwcaps[i].len))
      {
	if (t->disable)
	  GLRO(dl_riscv_cpu_features).extensions &= ~riscv_hwcaps[i].bit;
	else
	  GLRO(dl_riscv_cpu_features).extensions |= riscv_hwcaps[i].bit;
	break;
      }
}

attribute_hidden void
TUNABLE_CALLBACK (set_hwcaps) (tunable_val_t *valp)
{
  /* The current IFUNC selection is based on microbenchmarks in glibc.
     It should give the best performance for most workloads.  But other
     choices may have better performance for a particular workload or on
     the hardware which wasn't available when the selection was made.
     The environment variable:

     GLIBC_TUNABLES=glibc.cpu.hwcaps=-xxx,yyy,-zzz,....

     can be used to enable CPU/ARCH feature yyy, disable CPU/ARCH feature
     xxx and zzz, where the feature name is case-sensitive and has to
     match the ones defined below.  It can be used by glibc developers
     to tune for a new processor or override the IFUNC selection to
     improve performance for a particular workload.

     NOTE: the IFUNC selection may change over time.  Please check all
     multiarch implementations when experimenting.  */

  struct tunable_str_comma_state_t ts;
  tunable_str_comma_init (&ts, valp);

  struct tunable_str_comma_t n;
  while (tunable_str_comma_next (&ts, &n))
    if (n.len > 0)
      riscv_apply_tunable_hwcaps (&n);

  /* Ensure that the user has not enabled any unsupported features.
     We can only disable features that exist, not enable features
     that the hardware doesn't support.  */
  struct riscv_hwprobe pair;
  pair.key = RISCV_HWPROBE_KEY_IMA_EXT_0;
  if (__riscv_hwprobe (&pair, 1, 0, NULL, 0) == 0)
    GLRO(dl_riscv_cpu_features).extensions &= pair.value;
}
