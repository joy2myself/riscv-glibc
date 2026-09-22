ifeq ($(subdir),mathvec)

libmvec-sysdep_routines = \
	$(addprefix $(libmvec-vprefix)v_d_, $(libmvec-double-funcs)) \
	$(addsuffix _data, $(addprefix v_d_, $(libmvec-double-funcs-with-data))) \
	$(addprefix v_f_, $(libmvec-float-funcs)) \
	$(addsuffix _data, $(addprefix v_f_, $(libmvec-float-funcs-with-data))) \

libmvec-source-original-routines = \
  v_f_acosf_original \
  v_f_asinf_original

libmvec-sysdep_routines += $(libmvec-source-original-routines)

$(foreach f,$(libmvec-source-original-routines), \
  $(eval CFLAGS-$(f).c += -ffp-contract=off))

define riscv64-vector-cflags-template
CFLAGS-$(1).c += -march=rv64gcv -Wno-maybe-uninitialized -Wno-undef
endef

define riscv32-vector-cflags-template
CFLAGS-$(1).c += -march=rv32gcv -Wno-maybe-uninitialized -Wno-undef
endef

ifeq ($(config-machine),riscv64)
$(foreach f,$(libmvec-sysdep_routines), $(eval $(call riscv64-vector-cflags-template,$(f))))
endif

ifeq ($(config-machine),riscv32)
$(foreach f,$(libmvec-sysdep_routines), $(eval $(call riscv32-vector-cflags-template,$(f))))
endif

endif

ifeq ($(build-mathvec),yes)

bench-libmvec := $(bench-libmvec-double) $(bench-libmvec-float)

ifeq (${STATIC-BENCHTESTS},yes)
libmvec-benchtests = $(common-objpfx)mathvec/libmvec.a $(common-objpfx)math/libm.a
else
libmvec-benchtests = $(libmvec) $(libm)
endif

$(addprefix $(objpfx)bench-,$(bench-libmvec)): $(libmvec-benchtests)
$(foreach b,$(bench-libmvec),$(eval LDLIBS-bench-$(b) += $(libmvec)))

bench-libmvec-deps = $(..)benchtests/bench-libmvec-skeleton.c \
		     $(..)sysdeps/riscv/rvd/bench-libmvec-arch.h \
		     $(..)sysdeps/riscv/rvd/scripts/bench_libmvec_rvv.py \
		     bench-timing.h Makefile

$(objpfx)bench-float-%.c: $(bench-libmvec-deps)
	$(PYTHON) $(..)sysdeps/riscv/rvd/scripts/bench_libmvec_rvv.py $(basename $(@F)) > $@-tmp
	mv -f $@-tmp $@

$(objpfx)bench-double-%.c: $(bench-libmvec-deps)
	$(PYTHON) $(..)sysdeps/riscv/rvd/scripts/bench_libmvec_rvv.py $(basename $(@F)) > $@-tmp
	mv -f $@-tmp $@

define riscv64-bench-vector-cflags-template
CFLAGS-bench-$(1).c += -march=rv64gcv -Wno-missing-prototypes -Wno-undef
endef

define riscv32-bench-vector-cflags-template
CFLAGS-bench-$(1).c += -march=rv32gcv -Wno-missing-prototypes -Wno-undef
endef

ifeq ($(config-machine),riscv64)
$(foreach b,$(bench-libmvec), $(eval $(call riscv64-bench-vector-cflags-template,$(b))))
endif

ifeq ($(config-machine),riscv32)
$(foreach b,$(bench-libmvec), $(eval $(call riscv32-bench-vector-cflags-template,$(b))))
endif

ifneq ($(multi-arch),no)
ifeq ($(subdir)$(build-mathvec),mathyes)

ifneq ($(have-test-clang),yes)

tests += \
  $(libmvec-abi-func-tests)     \

endif

endif

ifeq ($(subdir),benchtests)

ifeq ($(build-mathvec),yes)

bench-libmvec := $(bench-libmvec-float) $(bench-libmvec-double)

endif

endif

endif

endif
