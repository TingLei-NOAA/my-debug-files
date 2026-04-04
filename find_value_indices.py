"""
Find grid indices for maximum values and values close to a target in an FV3 NetCDF field.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import netCDF4 as nc
import numpy as np


def infer_level_slice(var, time_idx: int, k_idx: int) -> tuple:
    slicer = []
    for dim_name, dim_size in zip(var.dimensions, var.shape):
        dim = dim_name.lower()
        if dim in {"time", "t", "nt"}:
            idx = time_idx
        elif dim.startswith("zaxis") or dim in {"lev", "pfull", "phalf", "nz"}:
            idx = k_idx
        elif dim.startswith("yaxis") or dim in {"grid_yt", "grid_y"}:
            slicer.append(slice(None))
            continue
        elif dim.startswith("xaxis") or dim in {"grid_xt", "grid_x"}:
            slicer.append(slice(None))
            continue
        else:
            raise SystemExit(f"Unsupported dim {dim_name} on variable {var.name}")
        if idx < 0 or idx >= dim_size:
            raise SystemExit(f"Index {idx} out of bounds for dim {dim_name} size {dim_size}")
        slicer.append(idx)
    return tuple(slicer)


def main() -> None:
    parser = argparse.ArgumentParser(description="Find max and near-target indices in a 2D FV3 field slice.")
    parser.add_argument("--netcdf", type=Path, required=True, help="Input FV3 NetCDF file")
    parser.add_argument("--var", required=True, help="Variable name")
    parser.add_argument("--time", type=int, default=0, help="0-based time index")
    parser.add_argument("--k", type=int, required=True, help="0-based vertical index")
    parser.add_argument("--target", type=float, default=0.35, help="Target value to search around")
    parser.add_argument("--tol", type=float, default=0.01, help="Absolute tolerance around target")
    parser.add_argument("--top-n", type=int, default=20, help="Number of closest points to report")
    args = parser.parse_args()

    with nc.Dataset(args.netcdf) as ds:
        if args.var not in ds.variables:
            raise SystemExit(f"Variable not found: {args.var}")
        var = ds.variables[args.var]
        slicer = infer_level_slice(var, args.time, args.k)
        field = np.asarray(var[slicer])

    if field.ndim != 2:
        raise SystemExit(f"Expected 2D field after slicing, got shape {field.shape}")

    max_j, max_i = np.unravel_index(int(np.nanargmax(field)), field.shape)
    max_val = float(np.nanmax(field))
    print(f"max_value={max_val:.10g} at i={max_i + 1}, j={max_j + 1} (1-based)")
    print(f"max_value={max_val:.10g} at i={max_i}, j={max_j} (0-based)")

    mask = np.isfinite(field)
    all_j, all_i = np.where(mask)
    all_vals = field[mask]
    diffs = np.abs(all_vals - args.target)
    order = np.argsort(diffs)

    print(f"\nTop {min(args.top_n, order.size)} points closest to target={args.target} (tol={args.tol}):")
    print("# rank value abs_diff i_1based j_1based i_0based j_0based")
    shown = 0
    for rank_idx in order:
        val = float(all_vals[rank_idx])
        diff = float(diffs[rank_idx])
        i0 = int(all_i[rank_idx])
        j0 = int(all_j[rank_idx])
        if diff <= args.tol or shown < args.top_n:
            shown += 1
            print(f"{shown} {val:.10g} {diff:.10g} {i0 + 1} {j0 + 1} {i0} {j0}")
        if shown >= args.top_n and diff > args.tol:
            break


if __name__ == "__main__":
    main()
