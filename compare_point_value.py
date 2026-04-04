"""
Compare a direct FV3 NetCDF value at one point with the maximum value in an extracted profile file.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import netCDF4 as nc
import numpy as np


def read_profile_max(path: Path) -> float:
    values = []
    in_complete = False
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith("#"):
            if "complete profile" in line.lower():
                in_complete = True
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        if in_complete or parts[0].isdigit():
            values.append(float(parts[1]))
    if not values:
        raise SystemExit(f"No profile values found in {path}")
    return float(np.max(values))


def infer_slicer(var, time_idx: int, k_idx: int, j_idx: int, i_idx: int) -> tuple:
    slicer = []
    for dim_name, dim_size in zip(var.dimensions, var.shape):
        dim = dim_name.lower()
        if dim in {"time", "t", "nt"}:
            idx = time_idx
        elif dim.startswith("zaxis") or dim in {"lev", "pfull", "phalf", "nz"}:
            idx = k_idx
        elif dim.startswith("yaxis") or dim in {"grid_yt", "grid_y"}:
            idx = j_idx
        elif dim.startswith("xaxis") or dim in {"grid_xt", "grid_x"}:
            idx = i_idx
        else:
            raise SystemExit(f"Unsupported dim {dim_name} on variable {var.name}")
        if idx < 0 or idx >= dim_size:
            raise SystemExit(f"Index {idx} out of bounds for dim {dim_name} size {dim_size}")
        slicer.append(idx)
    return tuple(slicer)


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare direct FV3 point value with extracted profile max.")
    parser.add_argument("--netcdf", type=Path, required=True, help="Input FV3 NetCDF file")
    parser.add_argument("--var", required=True, help="Variable name")
    parser.add_argument("--time", type=int, default=0, help="0-based time index")
    parser.add_argument("--k", type=int, required=True, help="0-based vertical index")
    parser.add_argument("--j", type=int, required=True, help="0-based j index")
    parser.add_argument("--i", type=int, required=True, help="0-based i index")
    parser.add_argument("--profile", type=Path, required=True, help="Extracted profile text file")
    args = parser.parse_args()

    with nc.Dataset(args.netcdf) as ds:
        if args.var not in ds.variables:
            raise SystemExit(f"Variable not found: {args.var}")
        var = ds.variables[args.var]
        slicer = infer_slicer(var, args.time, args.k, args.j, args.i)
        direct_val = float(np.asarray(var[slicer]).squeeze())

    profile_max = read_profile_max(args.profile)

    print(f"netcdf_file={args.netcdf}")
    print(f"variable={args.var}")
    print(f"point=(time={args.time}, k={args.k}, j={args.j}, i={args.i})")
    print(f"direct_value={direct_val:.10g}")
    print(f"profile_file={args.profile}")
    print(f"profile_max={profile_max:.10g}")
    print(f"difference(profile_max-direct)={profile_max - direct_val:.10g}")


if __name__ == "__main__":
    main()
