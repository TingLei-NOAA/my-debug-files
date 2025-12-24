"""
Plot model grid and filtering grid distributions with equal-i/j lines.

This script reads:
  - FV3 model grid spec (prefers grid_xt/grid_yt; falls back to grid_lont/grid_latt)
  - Filtering grid dumps: mgbf_filtering_grid_latlon_<rank>.txt (core nlon x nlat)
and produces a two-panel map to compare the grids.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import netCDF4 as nc


def parse_field(path: Path, nlon: int, nlat: int, input_order: str) -> tuple[np.ndarray, np.ndarray]:
    text = path.read_text()
    m = re.search(r"values:\s*\[(.*?)\]", text, re.DOTALL)
    if not m:
        raise ValueError(f"No values array in {path}")
    numbers = [float(x) for x in re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", m.group(1))]
    if len(numbers) < 2:
        raise ValueError(f"Not enough numbers in {path}")
    vals = np.asarray(numbers)
    if vals.size % 2:
        vals = vals[:-1]
    core = vals[: 2 * nlon * nlat]
    lon = core[0::2]
    lat = core[1::2]
    if input_order == "lonlat":
        lon = lon.reshape(nlat, nlon)
        lat = lat.reshape(nlat, nlon)
    else:  # latlon
        lon = lon.reshape(nlon, nlat).T
        lat = lat.reshape(nlon, nlat).T
    return lon, lat


def stitch(grids: list[tuple[Path, np.ndarray, np.ndarray, int]], tiles_x: int, tiles_y: int, nlon: int, nlat: int, rank_order: str):
    if len(grids) != tiles_x * tiles_y:
        raise ValueError(f"Expected {tiles_x*tiles_y} tiles, found {len(grids)}")
    full_lon = np.full((tiles_y * nlat, tiles_x * nlon), np.nan)
    full_lat = np.full_like(full_lon, np.nan)
    for path, lon_tile, lat_tile, rank in grids:
        row, col = rank_to_row_col(rank, tiles_x, tiles_y, rank_order)
        y0 = row * nlat
        x0 = col * nlon
        full_lon[y0:y0 + nlat, x0:x0 + nlon] = lon_tile
        full_lat[y0:y0 + nlat, x0:x0 + nlon] = lat_tile
    return full_lon, full_lat


def rank_to_row_col(rank: int, tiles_x: int, tiles_y: int, rank_order: str) -> tuple[int, int]:
    if rank_order == "col":
        col = rank // tiles_y
        row = rank % tiles_y
    elif rank_order == "col_flip_cols":
        col = tiles_x - 1 - (rank // tiles_y)
        row = rank % tiles_y
    elif rank_order == "col_flip_rows":
        col = rank // tiles_y
        row = tiles_y - 1 - (rank % tiles_y)
    elif rank_order == "col_rev":
        col = tiles_x - 1 - (rank // tiles_y)
        row = tiles_y - 1 - (rank % tiles_y)
    elif rank_order == "row":
        row = rank // tiles_x
        col = rank % tiles_x
    elif rank_order == "row_flip_cols":
        row = rank // tiles_x
        col = tiles_x - 1 - (rank % tiles_x)
    elif rank_order == "row_flip_rows":
        row = tiles_y - 1 - (rank // tiles_x)
        col = rank % tiles_x
    elif rank_order == "row_rev":
        row = tiles_y - 1 - (rank // tiles_x)
        col = tiles_x - 1 - (rank % tiles_x)
    else:
        raise ValueError(f"Unsupported rank_order {rank_order}")
    return row, col


def read_model_grid(fv3_spec: Path) -> tuple[np.ndarray, np.ndarray]:
    with nc.Dataset(fv3_spec) as ds:
        if "grid_xt" in ds.variables and "grid_yt" in ds.variables:
            lon_1d = np.array(ds.variables["grid_xt"][:])
            lat_1d = np.array(ds.variables["grid_yt"][:])
            lon, lat = np.meshgrid(lon_1d, lat_1d)
        else:
            lon = np.array(ds.variables["grid_lont"][:])
            lat = np.array(ds.variables["grid_latt"][:])
    mask = ~np.isfinite(lon) | ~np.isfinite(lat)
    lon[mask] = np.nan
    lat[mask] = np.nan
    return lon, lat


def normalize_lon_180(lon_deg: np.ndarray) -> np.ndarray:
    return (lon_deg + 180.0) % 360.0 - 180.0


def plot_lines(ax, lon: np.ndarray, lat: np.ndarray, i_stride: int, j_stride: int, lw: float, alpha: float):
    ny, nx = lon.shape
    for j in range(0, ny, j_stride):
        ax.plot(lon[j, :], lat[j, :], color="tab:gray", linewidth=lw, alpha=alpha)
    for i in range(0, nx, i_stride):
        ax.plot(lon[:, i], lat[:, i], color="tab:gray", linewidth=lw, alpha=alpha)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot model and filtering grid distributions with equal-i/j lines.")
    parser.add_argument("--pattern", default="mgbf_filtering_grid_latlon_*.txt", help="Glob for filtering tiles")
    parser.add_argument("--nlon", type=int, default=14, help="Core nlon per tile")
    parser.add_argument("--nlat", type=int, default=14, help="Core nlat per tile")
    parser.add_argument("--tiles-x", type=int, default=22, help="Tiles west->east")
    parser.add_argument("--tiles-y", type=int, default=22, help="Tiles south->north")
    parser.add_argument("--input-order", type=str, default="lonlat", choices=["lonlat", "latlon"], help="Storage order in filtering files")
    parser.add_argument(
        "--rank-order",
        type=str,
        default="row",
        choices=["col", "col_flip_cols", "col_flip_rows", "col_rev", "row", "row_flip_cols", "row_flip_rows", "row_rev"],
        help="Rank to tile mapping",
    )
    parser.add_argument("--fv3-spec", type=Path, default=Path("fv3_grid_spec"), help="FV3 grid spec NetCDF")
    parser.add_argument("--output", type=Path, default=Path("grid_distributions.png"), help="Output image")
    parser.add_argument("--model-i-stride", type=int, default=20, help="Stride for model equal-i lines")
    parser.add_argument("--model-j-stride", type=int, default=20, help="Stride for model equal-j lines")
    parser.add_argument("--filt-i-stride", type=int, default=10, help="Stride for filtering equal-i lines")
    parser.add_argument("--filt-j-stride", type=int, default=10, help="Stride for filtering equal-j lines")
    args = parser.parse_args()

    files = sorted(Path(".").glob(args.pattern))
    grids = []
    for path in files:
        m = re.search(r"_([0-9]+)(?:\.[^.]+)?$", path.name)
        if not m:
            print(f"[rank-parse-fail] {path.name}")
            continue
        rank = int(m.group(1))
        lon_tile, lat_tile = parse_field(path, args.nlon, args.nlat, args.input_order)
        grids.append((path, lon_tile, lat_tile, rank))

    lon_filt, lat_filt = stitch(grids, args.tiles_x, args.tiles_y, args.nlon, args.nlat, rank_order=args.rank_order)
    lon_filt = normalize_lon_180(lon_filt)
    lon_model, lat_model = read_model_grid(args.fv3_spec)
    lon_model = normalize_lon_180(lon_model)

    try:
        import matplotlib.pyplot as plt
    except ImportError as e:  # pragma: no cover
        raise SystemExit("matplotlib is required for plotting") from e

    fig, axes = plt.subplots(1, 2, figsize=(12, 6), dpi=160, sharex=True, sharey=True)
    ax0, ax1 = axes

    plot_lines(ax0, lon_model, lat_model, args.model_i_stride, args.model_j_stride, lw=0.4, alpha=0.7)
    ax0.set_title("Model grid (equal-i/j lines)")
    ax0.set_xlabel("Longitude")
    ax0.set_ylabel("Latitude")
    ax0.grid(True, linewidth=0.3, alpha=0.4)

    plot_lines(ax1, lon_filt, lat_filt, args.filt_i_stride, args.filt_j_stride, lw=0.5, alpha=0.7)
    ax1.set_title("Filtering grid (equal-i/j lines)")
    ax1.set_xlabel("Longitude")
    ax1.grid(True, linewidth=0.3, alpha=0.4)

    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output)
    plt.close(fig)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
