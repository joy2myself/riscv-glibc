/* RISC-V startup definitions.
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

#ifndef _RISCV_LIBC_START_H
#define _RISCV_LIBC_START_H

#include <sysdeps/generic/libc-start.h>

#ifndef SHARED

/* RISC-V cpu feature initialization uses __riscv_hwprobe, which consults
   GOT-based vDSO state. Delay initialization until after static PIE self
   relocation, but keep it before IFUNC IREL resolution so resolvers still
   see the final CPU feature state.  */
#undef ARCH_INIT_CPU_FEATURES
#define ARCH_INIT_CPU_FEATURES()

#undef ARCH_SETUP_IREL
#define ARCH_SETUP_IREL()						      \
  do									      \
    {									      \
      init_cpu_features (&_dl_riscv_cpu_features);			      \
      apply_irel ();							      \
    }									      \
  while (0)

#endif /* !SHARED */

#endif /* _RISCV_LIBC_START_H */
