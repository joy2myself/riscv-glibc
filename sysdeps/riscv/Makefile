ifeq ($(subdir),misc)
sysdep_headers += sys/asm.h
endif

ifeq ($(subdir),string)
gen-as-const-headers += locale-defines.sym
endif

ifeq ($(subdir),resolv)
CFLAGS-ns_rr_cursor_next.c += -fno-tree-vectorize
CFLAGS-res-noaaaa.c += -fno-tree-vectorize
CFLAGS-res_send.c += -fno-tree-vectorize
CFLAGS-dns-host.c += -fno-tree-vectorize
CFLAGS-dns-network.c += -fno-tree-vectorize
endif

ifeq ($(subdir),nss)
CFLAGS-files-hosts.c += -fno-tree-vectorize
CFLAGS-getaddrinfo.c += -fno-tree-vectorize
CFLAGS-getnameinfo.c += -fno-tree-vectorize
endif

ifeq ($(subdir),elf)
gen-as-const-headers += dl-link.sym
sysdep-dl-routines += \
  dl-get-cpu-features \
  # sysdep-dl-routines

ifneq ($(have-tunables),no)
tests += tst-hwcap-tunables
tst-hwcap-tunables-ARGS = -- $(host-test-program-cmd)
endif

$(objpfx)tst-glibc-hwcaps: $(objpfx)libmarkermod2-1.so \
  $(objpfx)libmarkermod3-1.so $(objpfx)libmarkermod4-1.so
$(objpfx)tst-glibc-hwcaps.out: \
  $(objpfx)libmarkermod2.so \
    $(objpfx)glibc-hwcaps/rva20u64/libmarkermod2.so \
  $(objpfx)libmarkermod3.so \
    $(objpfx)glibc-hwcaps/rva20u64/libmarkermod3.so \
    $(objpfx)glibc-hwcaps/rva22u64/libmarkermod3.so \
  $(objpfx)libmarkermod4.so \
    $(objpfx)glibc-hwcaps/rva20u64/libmarkermod4.so \
    $(objpfx)glibc-hwcaps/rva22u64/libmarkermod4.so \
    $(objpfx)glibc-hwcaps/rva23u64/libmarkermod4.so \

$(objpfx)glibc-hwcaps/rva20u64/libmarkermod2.so: $(objpfx)libmarkermod2-2.so
	$(make-target-directory)
	cp $< $@
$(objpfx)glibc-hwcaps/rva20u64/libmarkermod3.so: $(objpfx)libmarkermod3-2.so
	$(make-target-directory)
	cp $< $@
$(objpfx)glibc-hwcaps/rva22u64/libmarkermod3.so: $(objpfx)libmarkermod3-3.so
	$(make-target-directory)
	cp $< $@
$(objpfx)glibc-hwcaps/rva20u64/libmarkermod4.so: $(objpfx)libmarkermod4-2.so
	$(make-target-directory)
	cp $< $@
$(objpfx)glibc-hwcaps/rva22u64/libmarkermod4.so: $(objpfx)libmarkermod4-3.so
	$(make-target-directory)
	cp $< $@
$(objpfx)glibc-hwcaps/rva23u64/libmarkermod4.so: $(objpfx)libmarkermod4-4.so
	$(make-target-directory)
	cp $< $@

ifeq (no,$(build-hardcoded-path-in-tests))
# This is an ld.so.cache test, and RPATH/RUNPATH in the executable
# interferes with its test objectives.
tests-container += tst-glibc-hwcaps-cache
endif

endif

# RISC-V's assembler also needs to know about PIC as it changes the definition
# of some assembler macros.
ASFLAGS-.os += $(pic-ccflag)

ifeq (no,$(riscv-r-align))
ASFLAGS-.os += -Wa,-mno-relax
ASFLAGS-.o += -Wa,-mno-relax
sysdep-CFLAGS += -mno-relax
endif
