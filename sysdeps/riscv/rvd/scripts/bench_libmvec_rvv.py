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

struct {args_type}
{{
{arg_fields}
}};

struct variant
{{
  const char *name;
  size_t stride;
  double (*call) (const struct {args_type} *, size_t idx);
  const struct {args_type} *in;
}};

static const struct {args_type} in0 = {{
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

FUNC_CONFIGS: Dict[str, Dict[str, object]] = {
    "sincos": {
        "abi": "ret2",
        "input_alias": "sin",
        "scalar_impl": "sincos",
    },
    "sincosf": {
        "abi": "ret2",
        "input_alias": "sinf",
        "scalar_impl": "sincosf",
    },
    "sincospi": {
        "abi": "ret2",
        "input_alias": "sinpi",
        "scalar_impl": "sincospi",
    },
    "sincospif": {
        "abi": "ret2",
        "input_alias": "sinpif",
        "scalar_impl": "sincospif",
    },
}


def _func_config(func: str) -> Dict[str, object]:
    config: Dict[str, object] = {
        "abi": "ret1",
        "input_alias": func,
        "scalar_impl": "direct",
    }
    config.update(FUNC_CONFIGS.get(func, {}))
    return config


def _make_vector_variants(
    func: str,
    elem_width: int,
    variants: Iterable[Tuple[int, int]],
    config: Dict[str, object],
) -> Iterator[Variant]:
    """Yield vector variant metadata for the given function."""
    vtype_prefix = f"vfloat{elem_width}"
    setvl_prefix = f"e{elem_width}"
    for lmul, simdlen in variants:
        m_suffix = f"m{lmul}"
        vtype = f"{vtype_prefix}{m_suffix}_t"
        setvl_suffix = f"{setvl_prefix}{m_suffix}"
        abi = str(config["abi"])
        if abi == "ret2":
            abi_suffix = "vl8l8" if elem_width == 64 else "vl4l4"
        elif abi == "arg2":
            abi_suffix = "vv"
        else:
            abi_suffix = "v"
        asm_sym = f"_ZGVr{lmul}N{simdlen}{abi_suffix}_{func}"
        yield {
            "name": asm_sym,
            "vtype": vtype,
            "setvl_suffix": setvl_suffix,
            "simdlen": simdlen,
            "asm_sym": asm_sym,
            "kind": "vector",
            "abi": abi,
            "lmul": lmul,
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


def _load_input_data(
    func: str, config: Dict[str, object]
) -> Tuple[str, int, int]:
    """Load inputs for FUNC, returning initializer, count, and arity."""
    input_filename = f"{config['input_alias']}-inputs"
    source_root = Path(__file__).resolve().parents[4]
    input_paths = (
        source_root / "benchtests" / "libmvec" / input_filename,
        source_root / "benchtests" / input_filename,
    )
    input_path = next((path for path in input_paths if path.is_file()), None)
    if input_path is None:
        raise SystemExit(f"cannot find benchmark input file for {func}")
    with input_path.open(encoding="utf-8") as f:
        input_file = f.readlines()

    arity_line = next(
        (l for l in input_file if l.startswith("## args: ")), None
    )
    if arity_line is None:
        raise SystemExit("missing arity line in input file")
    arity = len(arity_line[len("## args: "):].strip().split(":"))
    if arity not in (1, 2):
        raise SystemExit(f"unsupported input arity: {arity}")

    rows = [
        [value.strip() for value in line.strip().split(",")]
        for line in input_file
        if line.strip() and not line.startswith("#")
    ]
    if not rows or any(len(row) != arity for row in rows):
        raise SystemExit("input row count does not match declared arity")
    in_vals = list(zip(*rows))
    in_data = ",\n".join(
        "{" + (", ".join(val for val in col) + "}") for col in in_vals
    )
    return in_data, len(rows), arity


def _vector_variant_snippets(
    idx: int,
    variant: Variant,
    arity: int,
    args_type: str,
) -> Tuple[str, str, str]:
    """Build prototype, caller, and table entry for a vector variant."""
    vtype = variant["vtype"]
    setvl_suffix = variant["setvl_suffix"]
    simdlen = variant["simdlen"]
    asm_sym = variant["asm_sym"]
    name_sym = variant["name"]
    abi = str(variant["abi"])

    if "float64" in vtype:
        elem_width = "64"
        elem_type = "f64"
    elif "float32" in vtype:
        elem_width = "32"
        elem_type = "f32"
    else:
        raise SystemExit(f"unsupported vtype: {vtype}")

    m_suffix = vtype.replace(f"vfloat{elem_width}", "").replace("_t", "")

    if abi == "ret2":
        proto = (
            f"extern void {name_sym}({vtype}, {vtype} *, {vtype} *) "
            f"__asm__(\"{asm_sym}\");"
        )
    else:
        call_args = ", ".join([vtype] * arity)
        proto = (
            f"extern {vtype} {name_sym}({call_args}) "
            f"__asm__(\"{asm_sym}\");"
        )
    vector_args = [
        f"  {vtype} arg{i} = "
        f"__riscv_vle{elem_width}_v_{elem_type}{m_suffix} "
        f"(&in->arg{i}[idx], vl);"
        for i in range(arity)
    ]
    if abi == "ret2":
        caller = (
            f"static double\n"
            f"call_{idx} (const struct {args_type} *in, size_t idx)\n"
            f"{{\n"
            f"  size_t vl = __riscv_vsetvl_{setvl_suffix} ({simdlen});\n"
            f"{chr(10).join(vector_args)}\n"
            f"  {vtype} r0, r1;\n"
            f"  {name_sym} (arg0, &r0, &r1);\n"
            f"  return __riscv_vfmv_f_s_{elem_type}{m_suffix}_{elem_type} "
            f"(r0) + __riscv_vfmv_f_s_{elem_type}{m_suffix}_{elem_type} (r1);\n"
            f"}}\n"
        )
    else:
        call_expr = f"{name_sym} (" + ", ".join(
            f"arg{i}" for i in range(arity)
        ) + ")"
        caller = (
            f"static double\n"
            f"call_{idx} (const struct {args_type} *in, size_t idx)\n"
            f"{{\n"
            f"  size_t vl = __riscv_vsetvl_{setvl_suffix} ({simdlen});\n"
            f"{chr(10).join(vector_args)}\n"
            f"  {vtype} r = {call_expr};\n"
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
    config: Dict[str, object],
    arity: int,
    args_type: str,
) -> Tuple[List[str], List[str], List[str]]:
    """Build prototypes, callers, and table entries for all variants."""
    protos: List[str] = []
    callers: List[str] = []
    variant_entries: List[str] = []

    variants_list = list(
        _make_vector_variants(func, elem_width, variants, config)
    )
    prefix = "v_d_" if elem_width == 64 else "v_f_"
    original = Path(__file__).resolve().parent.parent / (
        prefix + func + "_original.c"
    )
    if original.is_file():
        variants_list.extend(
            _make_vector_variants(
                func + "_original", elem_width, variants, config
            )
        )

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
            scalar_impl = str(config["scalar_impl"])
            if scalar_impl == "sincos":
                caller = (
                    f"static double\n"
                    f"call_{idx} (const struct {args_type} *in, size_t idx)\n"
                    f"{{\n"
                    f"  double s, c;\n"
                    f"  sincos (in->arg0[idx], &s, &c);\n"
                    f"  return s + c;\n"
                    f"}}\n"
                )
            elif scalar_impl == "sincosf":
                caller = (
                    f"static double\n"
                    f"call_{idx} (const struct {args_type} *in, size_t idx)\n"
                    f"{{\n"
                    f"  float s, c;\n"
                    f"  sincosf (in->arg0[idx], &s, &c);\n"
                    f"  return (double) s + (double) c;\n"
                    f"}}\n"
                )
            elif scalar_impl == "sincospi":
                caller = (
                    f"static double\n"
                    f"call_{idx} (const struct {args_type} *in, size_t idx)\n"
                    f"{{\n"
                    f"  return sinpi (in->arg0[idx]) + cospi (in->arg0[idx]);\n"
                    f"}}\n"
                )
            elif scalar_impl == "sincospif":
                caller = (
                    f"static double\n"
                    f"call_{idx} (const struct {args_type} *in, size_t idx)\n"
                    f"{{\n"
                    f"  return (double) sinpif (in->arg0[idx]) + "
                    f"(double) cospif (in->arg0[idx]);\n"
                    f"}}\n"
                )
            else:
                scalar_func = func if scalar_impl == "direct" else scalar_impl
                caller = (
                    f"static double\n"
                    f"call_{idx} (const struct {args_type} *in, size_t idx)\n"
                    f"{{\n"
                    f"  return {scalar_func} (" + ", ".join(
                        f"in->arg{i}[idx]" for i in range(arity)
                    ) + ");\n"
                    f"}}\n"
                )
            callers.append(caller)
            variant_entries.append(
                f'  {{"{variant["asm_sym"]}", 1, call_{idx}, &in0}}'
            )
            continue

        proto, caller, entry = _vector_variant_snippets(
            idx, variant, arity, args_type
        )
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
    func_config = _func_config(func)
    in_data, nelems, arity = _load_input_data(func, func_config)
    args_type = f"args_{arity}"
    arg_fields = "\n".join(
        f"  {stype} arg{i}[{nelems}];" for i in range(arity)
    )
    protos, callers, variant_entries = _build_variants(
        func,
        elem_width,
        config["variants"],
        func_config,
        arity,
        args_type,
    )

    print(
        TEMPLATE.format(
            stype=stype,
            nelems=nelems,
            args_type=args_type,
            arg_fields=arg_fields,
            in_data=in_data,
            protos="\n".join(protos),
            callers="\n".join(callers),
            variant_table=",\n".join(variant_entries),
            func=func,
        )
    )


if __name__ == "__main__":
    main(sys.argv[1])
