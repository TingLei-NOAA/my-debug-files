#!/usr/bin/env python3

import argparse
import re
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

    args = parser.parse_args()

    print("\nOpening files...")
    print(f"FILE 1: {args.file1}")
    print(f"FILE 2: {args.file2}")
    print(f"ATOL  : {args.atol:g}")
    print(f"RTOL  : {args.rtol:g}")

    ds1 = xr.open_dataset(
        args.file1,
        decode_times=False
    )

    ds2 = xr.open_dataset(
        args.file2,
        decode_times=False
    )

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

    # Final verdict, read by compare_run.sh.  Axis variables (xaxis_N, yaxis_N, zaxis_N, Time)
    # only label the dimensions: the FMS and regional writers name them differently
    # (e.g. yaxis_1 vs yaxis_2) and may leave their values unset, so they are reported but
    # do not decide the verdict.  Identical = every other variable is in both files and is
    # OK or SAME.
    def is_axis(name):
        return re.fullmatch(r"[xyz]axis_\d+", name) is not None or name == "Time"

    differing = [r["name"] for r in results
                 if r["status"] not in ("OK", "SAME") and not is_axis(r["name"])]
    missing1 = [n for n in only2 if not is_axis(n)]   # data variables absent from file 1
    missing2 = [n for n in only1 if not is_axis(n)]   # data variables absent from file 2
    axis_notes = ([f"{r['name']} ({r['status']})" for r in results
                   if r["status"] not in ("OK", "SAME") and is_axis(r["name"])]
                  + [f"{n} (only in file 1)" for n in only1 if is_axis(n)]
                  + [f"{n} (only in file 2)" for n in only2 if is_axis(n)])

    print("\n" + "=" * 140)
    print("SUMMARY (axis variables excluded from the verdict)")
    print("=" * 140)
    print(f"Data variables that differ   : {', '.join(differing) if differing else 'none'}")
    print(f"Data variables only in file 1: {', '.join(missing2) if missing2 else 'none'}")
    print(f"Data variables only in file 2: {', '.join(missing1) if missing1 else 'none'}")
    print(f"Axis variables, ignored      : {', '.join(axis_notes) if axis_notes else 'all match'}")
    if not differing and not missing1 and not missing2:
        print("OVERALL: IDENTICAL")
    else:
        print(f"OVERALL: DIFFERENT ({len(differing)} data variables differ, "
              f"{len(missing1) + len(missing2)} data variables in only one file)")


if __name__ == "__main__":
    main()
