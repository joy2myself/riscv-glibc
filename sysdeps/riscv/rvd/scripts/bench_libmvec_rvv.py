#!/usr/bin/python3
# Copyright (C) 2026 Free Software Foundation, Inc.
# This file is part of the GNU C Library.
#
# The GNU C Library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# The GNU C Library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with the GNU C Library; if not, see
# <https://www.gnu.org/licenses/>.

"""Generate RVV libmvec bench source that exercises all exported variants."""

import re
import sys
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Tuple, Union

TEMPLATE = r"""
#include <math.h>
#include <riscv_vector.h>

#pragma GCC diagnostic ignored "-Wmissing-prototypes"

struct args_1
{{
  {stype} arg0[{nelems}];
}};

struct variant
{{
  const char *name;
  size_t stride;
  double (*call) (const struct args_1 *, size_t idx);
  const struct args_1 *in;
}};

static const struct args_1 in0 = {{
{in_data}
}};

{protos}

{callers}

static const struct variant variants[] = {{
{variant_table}
}};

#define NUM_VARIANTS (sizeof (variants) / sizeof (variants[0]))
#define NUM_SAMPLES(i) ({nelems} / variants[i].stride)
#define VARIANT(i) (variants[i].name)
#define STRIDE (variants[v].stride)

static double volatile ret;

#define BENCH_FUNC(i, j) ({{ ret = variants[i].call (variants[i].in, (j) * variants[i].stride); }})
#define FUNCNAME "{func}"

#include <bench-libmvec-skeleton.c>
"""


def _parse_variants_macro(lines: List[str], macro: str) -> List[Tuple[int, int]]:
    variants: List[Tuple[int, int]] = []
    in_macro = False
    for line in lines:
        if line.startswith(f"#define {macro}"):
            in_macro = True
        if not in_macro:
            continue
        for lmul, simdlen in re.findall(r"X\s*\(\s*(\d+)\s*,\s*(\d+)\s*\)", line):
            variants.append((int(lmul), int(simdlen)))
        if in_macro and not line.rstrip().endswith("\\"):
            in_macro = False
            if variants:
                break
    if not variants:
        raise SystemExit(f"failed to parse {macro} from v_math_variants.h")
    return variants


def _load_variants_from_header() -> Dict[str, List[Tuple[int, int]]]:
    header = Path(__file__).resolve().parent.parent / "v_math_variants.h"
    try:
        lines = header.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise SystemExit(f"cannot read {header}: {exc}") from exc

    return {
        "double": _parse_variants_macro(lines, "RVV_LIBMVEC_VARIANTS_F64"),
        "float": _parse_variants_macro(lines, "RVV_LIBMVEC_VARIANTS_F32"),
    }


VARIANTS = _load_variants_from_header()

PRECISION_TABLE = {
    "double": {
        "stype": "double",
        "elem_width": 64,
        "variants": VARIANTS["double"],
    },
    "float": {
        "stype": "float",
        "elem_width": 32,
        "variants": VARIANTS["float"],
    },
}

Variant = Dict[str, Union[str, int]]


def _make_vector_variants(
    func: str, elem_width: int, variants: Iterable[Tuple[int, int]]
) -> Iterator[Variant]:
    """Yield vector variant metadata for the given function."""
    vtype_prefix = f"vfloat{elem_width}"
    setvl_prefix = f"e{elem_width}"
    for lmul, simdlen in variants:
        m_suffix = f"m{lmul}"
        vtype = f"{vtype_prefix}{m_suffix}_t"
        setvl_suffix = f"{setvl_prefix}{m_suffix}"
        asm_sym = f"_ZGVr{lmul}N{simdlen}v_{func}"
        yield {
            "name": asm_sym,
            "vtype": vtype,
            "setvl_suffix": setvl_suffix,
            "simdlen": simdlen,
            "asm_sym": asm_sym,
            "kind": "vector",
        }


def _parse_bench_name(name: str) -> Tuple[str, str]:
    """Parse a bench target name into (precision, function)."""
    # Accept names like "double-log" (bench list) or "bench-double-log"
    # (target name in benchtests).
    if name.startswith("bench-"):
        name = name[len("bench-"):]

    parts = name.split("-")
    if len(parts) < 2:
        raise SystemExit(f"unexpected name format: {name}")

    prec = parts[0]
    func = "_".join(parts[1:])
    return prec, func


def _load_input_data(func: str) -> Tuple[str, int]:
    """Load inputs for FUNC, returning the initializer and element count."""
    input_filename = f"{func}-inputs"
    with open(
        f"../benchtests/libmvec/{input_filename}",
        encoding="utf-8",
    ) as f:
        input_file = f.readlines()

    arity_line = [l for l in input_file if l.startswith("## args: ")]
    if not arity_line:
        raise SystemExit("missing arity line in input file")

    in_vals = (l.strip() for l in input_file if l and not l.startswith("#"))
    in_vals = (l.split(", ") for l in in_vals)
    in_vals = list(zip(*in_vals))
    in_data = ",\n".join(
        "{" + (", ".join(val for val in col) + "}") for col in in_vals
    )
    return in_data, len(in_vals[0])


def _vector_variant_snippets(
    idx: int,
    variant: Variant,
) -> Tuple[str, str, str]:
    """Build prototype, caller, and table entry for a vector variant."""
    vtype = variant["vtype"]
    setvl_suffix = variant["setvl_suffix"]
    simdlen = variant["simdlen"]
    asm_sym = variant["asm_sym"]
    name_sym = variant["name"]

    if "float64" in vtype:
        elem_width = "64"
        elem_type = "f64"
    elif "float32" in vtype:
        elem_width = "32"
        elem_type = "f32"
    else:
        raise SystemExit(f"unsupported vtype: {vtype}")

    m_suffix = vtype.replace(f"vfloat{elem_width}", "").replace("_t", "")

    proto = (
        f"extern {vtype} {name_sym}({vtype}) __asm__(\"{asm_sym}\");"
    )
    caller = (
        f"static double\n"
        f"call_{idx} (const struct args_1 *in, size_t idx)\n"
        f"{{\n"
        f"  size_t vl = __riscv_vsetvl_{setvl_suffix} ({simdlen});\n"
        f"  {vtype} x = "
        f"__riscv_vle{elem_width}_v_{elem_type}{m_suffix} "
        f"(&in->arg0[idx], vl);\n"
        f"  {vtype} r = {name_sym} (x);\n"
        f"  return __riscv_vfmv_f_s_{elem_type}{m_suffix}_{elem_type} "
        f"(r);\n"
        f"}}\n"
    )
    entry = f'  {{"{asm_sym}", {simdlen}, call_{idx}, &in0}}'
    return proto, caller, entry


def _build_variants(
    func: str,
    elem_width: int,
    variants: Iterable[Tuple[int, int]],
) -> Tuple[List[str], List[str], List[str]]:
    """Build prototypes, callers, and table entries for all variants."""
    protos: List[str] = []
    callers: List[str] = []
    variant_entries: List[str] = []

    variants_list = list(_make_vector_variants(func, elem_width, variants))
    prefix = "v_d_" if elem_width == 64 else "v_f_"
    original = Path(__file__).resolve().parent.parent / (
        prefix + func + "_original.c"
    )
    if original.is_file():
        variants_list.extend(_make_vector_variants(func + "_original", elem_width, variants))

    variants_list.append(
        {
            "name": "scalar",
            "asm_sym": func,
            "simdlen": 1,
            "kind": "scalar",
        }
    )

    for idx, variant in enumerate(variants_list):
        if variant["kind"] == "scalar":
            callers.append(
                f"static double\n"
                f"call_{idx} (const struct args_1 *in, size_t idx)\n"
                f"{{\n"
                f"  return {func} (in->arg0[idx]);\n"
                f"}}\n"
            )
            variant_entries.append(
                f'  {{"{variant["asm_sym"]}", 1, call_{idx}, &in0}}'
            )
            continue

        proto, caller, entry = _vector_variant_snippets(idx, variant)
        protos.append(proto)
        callers.append(caller)
        variant_entries.append(entry)

    return protos, callers, variant_entries

def main(name: str):
    prec, func = _parse_bench_name(name)

    config = PRECISION_TABLE.get(prec)
    if not config:
        raise SystemExit(f"unsupported precision: {prec}")

    stype = config["stype"]
    elem_width = config["elem_width"]
    in_data, nelems = _load_input_data(func)
    protos, callers, variant_entries = _build_variants(
        func,
        elem_width,
        config["variants"],
    )

    print(
        TEMPLATE.format(
            stype=stype,
            nelems=nelems,
            in_data=in_data,
            protos="\n".join(protos),
            callers="\n".join(callers),
            variant_table=",\n".join(variant_entries),
            func=func,
        )
    )


if __name__ == "__main__":
    main(sys.argv[1])
