#!/usr/bin/env python3

import argparse
import re
import sys
import numpy as np
import xarray as xr


def compare_variable(name, a, b, atol, rtol):
    result = {
        "name": name,
        "status": "",
        "max_abs": np.nan,
        "rms": np.nan,
        "max_rel": np.nan,
        "nfailed": 0,
        "ntotal": 0,
        "max_index": None,
        "value1": np.nan,
        "value2": np.nan,
    }

    # Check shape
    if a.shape != b.shape:
        result["status"] = f"SHAPE DIFFERENT: {a.shape} vs {b.shape}"
        return result

    # Skip non-numeric variables
    if not (
        np.issubdtype(a.dtype, np.number)
        and np.issubdtype(b.dtype, np.number)
    ):
        same = np.array_equal(a.values, b.values)
        result["status"] = "SAME" if same else "DIFFERENT"
        return result

    x = np.asarray(a.values)
    y = np.asarray(b.values)

    # Check NaN pattern
    nan_pattern_same = np.array_equal(np.isnan(x), np.isnan(y))

    # Points where both values are finite
    valid = np.isfinite(x) & np.isfinite(y)

    if not np.any(valid):
        result["status"] = "SAME" if nan_pattern_same else "NaN PATTERN DIFF"
        return result

    # Difference array, preserving original shape
    diff = np.full(x.shape, np.nan, dtype=np.float64)
    diff[valid] = x[valid].astype(np.float64) - y[valid].astype(np.float64)

    absdiff = np.abs(diff)

    # Maximum difference and its location
    flat_index = np.nanargmax(absdiff)
    max_index = np.unravel_index(flat_index, x.shape)

    result["max_index"] = max_index
    result["max_abs"] = absdiff[max_index]
    result["value1"] = x[max_index]
    result["value2"] = y[max_index]

    # Valid flattened values
    xv = x[valid].astype(np.float64)
    yv = y[valid].astype(np.float64)

    d = xv - yv
    ad = np.abs(d)

    result["rms"] = np.sqrt(np.mean(d * d))

    # Relative difference
    scale = np.maximum(np.abs(xv), np.abs(yv))
    tiny = np.finfo(np.float64).tiny

    rel = ad / np.maximum(scale, tiny)
    result["max_rel"] = np.max(rel)

    # Tolerance test
    failed = ad > (atol + rtol * np.abs(yv))

    result["nfailed"] = np.count_nonzero(failed)
    result["ntotal"] = len(xv)

    if result["nfailed"] == 0 and nan_pattern_same:
        result["status"] = "OK"
    else:
        result["status"] = "DIFF"

    return result


def io_name_map(yaml_path):
    """{name in the file: long name}, from "field io names" in the top-level "output" block
    of a JEDI YAML (empty when no YAML is given or it has no such block)."""
    if not yaml_path:
        return {}
    try:
        import yaml
    except ImportError:
        sys.exit("compare-nc.py: --yaml1/--yaml2 need PyYAML (python -m pip install pyyaml)")
    with open(yaml_path) as f:
        doc = yaml.safe_load(f) or {}
    output = doc.get("output") or {}
    names = output.get("field io names") or {}
    return {str(io_name): str(long_name) for long_name, io_name in names.items()}


def to_long_names(ds, mapping, label):
    """Rename the variables of ds that use an io name back to their long name."""
    rename = {}
    for io_name, long_name in mapping.items():
        if io_name == long_name or io_name not in ds.variables:
            continue
        if long_name in ds.variables:
            print(f"WARNING: {label} has both {io_name} and {long_name}; {io_name} is not renamed")
            continue
        rename[io_name] = long_name
    if rename:
        print(f"\n{label}: io names mapped to long names (from its YAML's field io names):")
        for io_name, long_name in sorted(rename.items()):
            print(f"  {io_name:20s} -> {long_name}")
        ds = ds.rename_vars(rename)
    return ds


def main():
    parser = argparse.ArgumentParser(
        description="""
Compare two NetCDF files variable-by-variable.

Results are sorted from the largest maximum absolute difference
to the smallest.
"""
    )

    parser.add_argument(
        "file1",
        help="First NetCDF file"
    )

    parser.add_argument(
        "file2",
        help="Second NetCDF file"
    )

    parser.add_argument(
        "--atol",
        type=float,
        default=0.0,
        help="Absolute tolerance (default: 0)"
    )

    parser.add_argument(
        "--rtol",
        type=float,
        default=0.0,
        help="Relative tolerance (default: 0)"
    )

    parser.add_argument(
        "--yaml1",
        help="JEDI YAML that wrote file 1: its output \"field io names\" map file names to long names"
    )

    parser.add_argument(
        "--yaml2",
        help="JEDI YAML that wrote file 2 (same use as --yaml1)"
    )

    args = parser.parse_args()

    print("\nOpening files...")
    print(f"FILE 1: {args.file1}")
    print(f"FILE 2: {args.file2}")
    print(f"ATOL  : {args.atol:g}")
    print(f"RTOL  : {args.rtol:g}")
    if args.yaml1 or args.yaml2:
        print(f"YAML 1: {args.yaml1 or '-'}")
        print(f"YAML 2: {args.yaml2 or '-'}")

    ds1 = xr.open_dataset(
        args.file1,
        decode_times=False
    )

    ds2 = xr.open_dataset(
        args.file2,
        decode_times=False
    )

    # Compare by long name: a file written with "field io names" stores io names (T, delp, ...)
    ds1 = to_long_names(ds1, io_name_map(args.yaml1), "file 1")
    ds2 = to_long_names(ds2, io_name_map(args.yaml2), "file 2")

    vars1 = set(ds1.variables)
    vars2 = set(ds2.variables)

    common_vars = sorted(vars1 & vars2)
    only1 = sorted(vars1 - vars2)
    only2 = sorted(vars2 - vars1)

    print("\n------------------------------------------------------------")
    print("VARIABLE SUMMARY")
    print("------------------------------------------------------------")

    print(f"Variables in file 1 : {len(vars1)}")
    print(f"Variables in file 2 : {len(vars2)}")
    print(f"Common variables    : {len(common_vars)}")

    if only1:
        print("\nOnly in file 1:")
        for name in only1:
            print(f"  {name}")

    if only2:
        print("\nOnly in file 2:")
        for name in only2:
            print(f"  {name}")

    print("\nComparing common variables...")

    results = []

    for name in common_vars:
        print(f"  comparing {name}")
        result = compare_variable(
            name,
            ds1[name],
            ds2[name],
            args.atol,
            args.rtol,
        )
        results.append(result)

    # Sort largest max absolute difference first.
    # Variables without a numerical max_abs go to the bottom.
    results.sort(
        key=lambda r:
        -np.inf if np.isnan(r["max_abs"]) else r["max_abs"],
        reverse=True,
    )

    print("\n")
    print("=" * 140)
    print("RESULTS SORTED BY MAXIMUM ABSOLUTE DIFFERENCE")
    print("=" * 140)

    header = (
        f"{'Variable':30s} "
        f"{'Status':8s} "
        f"{'MaxAbs':>14s} "
        f"{'RMS':>14s} "
        f"{'MaxRel':>14s} "
        f"{'Nfail':>16s}"
    )

    print(header)
    print("-" * 140)

    for r in results:

        if np.isnan(r["max_abs"]):
            print(
                f'{r["name"]:30s} '
                f'{r["status"]}'
            )

        else:
            print(
                f'{r["name"]:30s} '
                f'{r["status"]:8s} '
                f'{r["max_abs"]:14.6e} '
                f'{r["rms"]:14.6e} '
                f'{r["max_rel"]:14.6e} '
                f'{r["nfailed"]:8d}/{r["ntotal"]:<7d}'
            )

    print("\n")
    print("=" * 140)
    print("LOCATIONS OF MAXIMUM DIFFERENCE")
    print("=" * 140)

    for r in results:

        if np.isnan(r["max_abs"]):
            continue

        print(f"\nVariable: {r['name']}")
        print(f"  status     = {r['status']}")
        print(f"  max_abs    = {r['max_abs']:.12e}")
        print(f"  max_rel    = {r['max_rel']:.12e}")
        print(f"  rms        = {r['rms']:.12e}")
        print(f"  index      = {r['max_index']}")
        print(f"  file1 value= {r['value1']}")
        print(f"  file2 value= {r['value2']}")

    ds1.close()
    ds2.close()

    # Final verdict, read by compare_run.sh.  File 1 is the control: identical means every
    # control data variable is also in file 2 and is OK or SAME.  Extra variables in file 2
    # are listed but do not fail.  Axis variables (xaxis_N, yaxis_N, zaxis_N, Time) only
    # label the dimensions: the FMS and regional writers name them differently (e.g. yaxis_1
    # vs yaxis_2) and may leave their values unset, so they are listed but ignored.
    def is_axis(name):
        return re.fullmatch(r"[xyz]axis_\d+", name) is not None or name == "Time"

    differing = [r["name"] for r in results
                 if r["status"] not in ("OK", "SAME") and not is_axis(r["name"])]
    missing = [n for n in only1 if not is_axis(n)]    # control variables not in file 2
    extra = [n for n in only2 if not is_axis(n)]      # file 2 variables not in the control
    axis_notes = ([f"{r['name']} ({r['status']})" for r in results
                   if r["status"] not in ("OK", "SAME") and is_axis(r["name"])]
                  + [f"{n} (only in file 1)" for n in only1 if is_axis(n)]
                  + [f"{n} (only in file 2)" for n in only2 if is_axis(n)])

    print("\n" + "=" * 140)
    print("SUMMARY (file 1 = control; axis variables excluded from the verdict)")
    print("=" * 140)
    print(f"Control variables that differ    : {', '.join(differing) if differing else 'none'}")
    print(f"Control variables missing in 2   : {', '.join(missing) if missing else 'none'}")
    print(f"Extra variables in 2 (not judged): {', '.join(extra) if extra else 'none'}")
    print(f"Axis variables (not judged)      : {', '.join(axis_notes) if axis_notes else 'all match'}")
    if not differing and not missing:
        print("OVERALL: IDENTICAL")
    else:
        print(f"OVERALL: DIFFERENT ({len(differing)} control variables differ, "
              f"{len(missing)} control variables missing in file 2)")


if __name__ == "__main__":
    main()
