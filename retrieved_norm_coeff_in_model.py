"""
Interpolate sparse profile values onto the full FV3 model grid.

Inputs:
  - profiles_out/profile_subdomain_####.txt files (complete profile section)
  - FV3 dyn NetCDF (defines xaxis_1, yaxis_1, zaxis_1)

Output:
  - NetCDF with interpolated profile values on (zaxis_1, yaxis_1, xaxis_1)
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import netCDF4 as nc

try:
    from scipy.interpolate import LinearNDInterpolator
except ImportError as e:  # pragma: no cover
    raise SystemExit("scipy is required for bilinear interpolation (LinearNDInterpolator)") from e


def read_profile_file(path: Path) -> tuple[int, int, int, np.ndarray]:
    idx = None
    i_idx = None
    j_idx = None
    data = []
    in_full = False
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith("#"):
            if "complete profile" in line:
                in_full = True
            continue
        if idx is None:
            parts = line.split()
            if len(parts) < 4:
                raise SystemExit(f"Bad header line in {path}")
            idx = int(parts[0])
            i_idx = int(parts[2])
            j_idx = int(parts[3])
            continue
        if in_full:
            parts = line.split()
            if len(parts) < 2:
                continue
            data.append(float(parts[1]))
    if idx is None or i_idx is None or j_idx is None:
        raise SystemExit(f"Missing header in {path}")
    if not data:
        raise SystemExit(f"No complete profile data found in {path}")
    return idx, i_idx, j_idx, np.asarray(data, dtype=float)


def main() -> None:
    parser = argparse.ArgumentParser(description="Interpolate profiles onto full FV3 grid.")
    parser.add_argument("--profiles-dir", type=Path, default=Path("profiles_out"), help="Directory with profile_subdomain_*.txt")
    parser.add_argument("--fv3-dyn", type=Path, required=True, help="FV3 dynamics NetCDF (grid definition)")
    parser.add_argument("--output", type=Path, default=Path("norm_coeff_full.nc"), help="Output NetCDF")
    args = parser.parse_args()

    files = sorted(args.profiles_dir.glob("profile_subdomain_*.txt"))
    if not files:
        raise SystemExit(f"No profile files found in {args.profiles_dir}")

    indices = []
    ij = []
    profiles = []
    for path in files:
        idx, i_idx, j_idx, prof = read_profile_file(path)
        indices.append(idx)
        ij.append((i_idx, j_idx))
        profiles.append(prof)

    profiles = np.asarray(profiles, dtype=float)
    with nc.Dataset(args.fv3_dyn) as ds:
        if "xaxis_1" not in ds.dimensions or "yaxis_1" not in ds.dimensions or "zaxis_1" not in ds.dimensions:
            raise SystemExit("Missing xaxis_1/yaxis_1/zaxis_1 dimensions in fv3 dyn")
        nx = len(ds.dimensions["xaxis_1"])
        ny = len(ds.dimensions["yaxis_1"])
        nz = len(ds.dimensions["zaxis_1"])
        if profiles.shape[1] != nz:
            raise SystemExit(f"Profile size {profiles.shape[1]} does not match zaxis_1 {nz}")
        if "xaxis_1" not in ds.variables or "yaxis_1" not in ds.variables:
            raise SystemExit("Missing xaxis_1/yaxis_1 coordinate variables in fv3 dyn")
        xcoord = np.array(ds.variables["xaxis_1"][:])
        ycoord = np.array(ds.variables["yaxis_1"][:])

    pts = []
    for i_idx, j_idx in ij:
        if i_idx < 1 or i_idx > nx or j_idx < 1 or j_idx > ny:
            raise SystemExit(f"Index out of range: i={i_idx} j={j_idx} for nx={nx} ny={ny}")
        pts.append((xcoord[i_idx - 1], ycoord[j_idx - 1]))
    pts = np.asarray(pts, dtype=float)

    grid_x, grid_y = np.meshgrid(xcoord, ycoord)
    full = np.full((nz, ny, nx), np.nan, dtype=float)

    for k in range(nz):
        values = profiles[:, k]
        interp = LinearNDInterpolator(pts, values)
        full[k, :, :] = interp(grid_x, grid_y)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with nc.Dataset(args.output, "w") as out:
        out.createDimension("zaxis_1", nz)
        out.createDimension("yaxis_1", ny)
        out.createDimension("xaxis_1", nx)
        out.createVariable("xaxis_1", xcoord.dtype, ("xaxis_1",))[:] = xcoord
        out.createVariable("yaxis_1", ycoord.dtype, ("yaxis_1",))[:] = ycoord
        out.createVariable("zaxis_1", "i4", ("zaxis_1",))[:] = np.arange(1, nz + 1)
        var = out.createVariable("norm_coeff", "f4", ("zaxis_1", "yaxis_1", "xaxis_1"), zlib=True, complevel=2)
        var[:] = full
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
