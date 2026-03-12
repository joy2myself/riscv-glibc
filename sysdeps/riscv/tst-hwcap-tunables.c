/* Test RISC-V CPU feature tunables.
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

#include <array_length.h>
#include <getopt.h>
#include <spawn.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <cpu-features.h>
#include <profile-detection.h>
#include <support/check.h>
#include <support/support.h>
#include <support/xunistd.h>
#include <sys/hwprobe.h>
#include <sys/wait.h>

/* Nonzero if the program gets called via exec.  */
#define CMDLINE_OPTIONS \
  { "restart", no_argument, &restart, 1 }, \
  { "disable", no_argument, &disable, 1 },

static int restart;
static int disable;

struct riscv_hwcap_name
{
  const char *name;
  uint64_t bit;
  bool alias;
};

#define RISCV_HWCAP_ENTRY(name, bit) { #name, bit, false },
#define RISCV_HWCAP_ALIAS(str, bit) { str, bit, true },
static const struct riscv_hwcap_name riscv_hwcaps[] =
  {
#include <cpu-hwcaps-list.h>
  };
#undef RISCV_HWCAP_ENTRY
#undef RISCV_HWCAP_ALIAS

static bool
hwcap_active (uint64_t extensions, uint64_t bit)
{
  return (extensions & bit) == bit;
}

static int
spawn_runner (pid_t *pid, char *const *argv)
{
  if (strchr (argv[0], '/') != NULL)
    return posix_spawn (pid, argv[0], NULL, NULL, argv, environ);

  return posix_spawnp (pid, argv[0], NULL, NULL, argv, environ);
}

static uint64_t
get_hwprobe_extensions (void)
{
  struct riscv_hwprobe pair = { .key = RISCV_HWPROBE_KEY_IMA_EXT_0 };

  TEST_COMPARE (__riscv_hwprobe (&pair, 1, 0, NULL, 0), 0);
  TEST_COMPARE (pair.key, RISCV_HWPROBE_KEY_IMA_EXT_0);
  return pair.value;
}

static const struct riscv_hwcap_name *
find_hwcap_name (const char *name)
{
  for (size_t i = 0; i < array_length (riscv_hwcaps); ++i)
    if (strcmp (name, riscv_hwcaps[i].name) == 0)
      return &riscv_hwcaps[i];

  FAIL_EXIT1 ("unknown hwcap name: %s\n", name);
  __builtin_unreachable ();
}

static void
print_active_hwcaps (const char *label, uint64_t extensions)
{
  printf ("%s:\n", label);
  for (size_t i = 0; i < array_length (riscv_hwcaps); ++i)
    if ((!riscv_hwcaps[i].alias
	 || strcmp (riscv_hwcaps[i].name, "vector") != 0)
	&& hwcap_active (extensions, riscv_hwcaps[i].bit))
      printf ("  %s\n", riscv_hwcaps[i].name);
}

_Noreturn static void
handle_restart (int argc, char *argv[])
{
  TEST_VERIFY_EXIT (argc == 1);

  const struct riscv_hwcap_name *hwcap = find_hwcap_name (argv[0]);
  const struct cpu_features *cpu_features = __get_cpu_features ();
  bool supported = hwcap_active (get_hwprobe_extensions (), hwcap->bit);
  bool expected = disable ? false : supported;

  TEST_COMPARE (hwcap_active (cpu_features->extensions, hwcap->bit), expected);
  _exit (EXIT_SUCCESS);
}

static void
run_test (char *const *spargs, const struct riscv_hwcap_name *hwcap,
	  uint64_t hwprobe_extensions, bool test_disable)
{
  size_t spargs_len = 0;
  while (spargs[spargs_len] != NULL)
    ++spargs_len;

  char **argv = xmalloc ((spargs_len + 3) * sizeof (*argv));
  size_t i = 0;
  for (; i < spargs_len; ++i)
    argv[i] = spargs[i];

  if (test_disable)
    argv[i++] = (char *) "--disable";

  argv[i++] = (char *) hwcap->name;
  argv[i] = NULL;

  bool supported = hwcap_active (hwprobe_extensions, hwcap->bit);
  printf ("info: checking filter %s%s (expect bit %s)\n",
	  test_disable ? "-" : "", hwcap->name,
	  test_disable ? "disabled" : supported ? "enabled" : "masked");

  char *saved_tunable = getenv ("GLIBC_TUNABLES");
  if (saved_tunable != NULL)
    saved_tunable = xstrdup (saved_tunable);
  char *tunable = xasprintf ("glibc.cpu.hwcaps=%s%s",
			     test_disable ? "-" : "", hwcap->name);
  TEST_COMPARE (setenv ("GLIBC_TUNABLES", tunable, 1), 0);

  pid_t pid;
  TEST_COMPARE (spawn_runner (&pid, argv), 0);

  if (saved_tunable != NULL)
    {
      TEST_COMPARE (setenv ("GLIBC_TUNABLES", saved_tunable, 1), 0);
      free (saved_tunable);
    }
  else
    TEST_COMPARE (unsetenv ("GLIBC_TUNABLES"), 0);

  int status;
  TEST_COMPARE (xwaitpid (pid, &status, 0), pid);
  TEST_VERIFY (WIFEXITED (status));
  TEST_VERIFY (!WIFSIGNALED (status));
  TEST_COMPARE (WEXITSTATUS (status), 0);

  free (argv);
  free (tunable);
}

static int
do_test (int argc, char *argv[])
{
  if (restart)
    handle_restart (argc - 1, &argv[1]);

  TEST_VERIFY_EXIT (argc >= 2);

  char **spargs = xmalloc ((argc + 2) * sizeof (*spargs));
  int i = 0;
  for (; i < argc - 1; ++i)
    spargs[i] = argv[i + 1];
  spargs[i++] = (char *) "--direct";
  spargs[i++] = (char *) "--restart";
  spargs[i] = NULL;

  const struct cpu_features *cpu_features = __get_cpu_features ();
  uint64_t hwprobe_extensions = get_hwprobe_extensions ();

  printf ("mvendorid: 0x%llx\n",
	  (unsigned long long) cpu_features->mvendorid);
  printf ("marchid: 0x%llx\n",
	  (unsigned long long) cpu_features->marchid);
  printf ("mimpid: 0x%llx\n",
	  (unsigned long long) cpu_features->mimpid);
  printf ("extensions: 0x%llx\n",
	  (unsigned long long) cpu_features->extensions);
  printf ("cpuperf: 0x%llx\n",
	  (unsigned long long) cpu_features->cpuperf);
  printf ("profile_level: %u\n", cpu_features->profile_level);
  printf ("zicclsm: %d\n",
	  riscv_hwprobe_have_zicclsm (cpu_features->cpuperf));
  printf ("zicboz_block_size: %u\n", cpu_features->zicboz_block_size);
  printf ("zicbop_block_size: %u\n", cpu_features->zicbop_block_size);
  printf ("fast_unaligned: %d\n",
	  riscv_hwprobe_fast_unaligned (cpu_features->cpuperf));
  printf ("hwprobe extensions: 0x%llx\n",
	  (unsigned long long) hwprobe_extensions);
  print_active_hwcaps ("hardware-supported glibc.cpu.hwcaps",
		       hwprobe_extensions);
  print_active_hwcaps ("active glibc.cpu.hwcaps", cpu_features->extensions);

  TEST_VERIFY ((cpu_features->extensions & ~hwprobe_extensions) == 0);

  for (size_t n = 0; n < array_length (riscv_hwcaps); ++n)
    {
      run_test (spargs, &riscv_hwcaps[n], hwprobe_extensions, false);
      if (hwcap_active (hwprobe_extensions, riscv_hwcaps[n].bit))
	run_test (spargs, &riscv_hwcaps[n], hwprobe_extensions, true);
    }

  free (spargs);
  return 0;
}

#define TEST_FUNCTION_ARGV do_test
#include <support/test-driver.c>
