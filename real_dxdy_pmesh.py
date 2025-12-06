"""
Assemble MGBF filtering grid dumps, stitch tiles by MPI rank (column-major), and
plot dx/dy using pcolormesh to avoid contour triangulation artifacts.

Assumptions:
  - Filenames: mgbf_filtering_grid_latlon_<rank>.txt (rank 0-based)
  - Each file contains a 14x14 core (first 196 pairs are the core)
  - Tile layout: 22 x 22 tiles; ranks increase south->north (fast), west->east (slow)

Outputs (default dir: dr-figures):
  real_dx_km_pmesh.png
  real_dy_km_pmesh.png
  real_dx_over_dy_pmesh.png
  stitched lon/lat diagnostics (if --dump-stitched)
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np


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
        lon = lon.reshape(nlon, nlat).T
        lat = lat.reshape(nlon, nlat).T
    else:
        lon = lon.reshape(nlat, nlon)
        lat = lat.reshape(nlat, nlon)
    return lon, lat


def haversine(lon1: np.ndarray, lat1: np.ndarray, lon2: np.ndarray, lat2: np.ndarray) -> np.ndarray:
    lon1_rad = np.radians(lon1)
    lon2_rad = np.radians(lon2)
    lat1_rad = np.radians(lat1)
    lat2_rad = np.radians(lat2)
    dlon = lon2_rad - lon1_rad
    dlat = lat2_rad - lat1_rad
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1_rad) * np.cos(lat2_rad) * np.sin(dlon / 2) ** 2
    a = np.clip(a, 0.0, 1.0)
    c = 2 * np.arcsin(np.sqrt(a))
    return 6_371_000.0 * c


def compute_dx_dy(lon: np.ndarray, lat: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    dx = haversine(lon[:, :-1], lat[:, :-1], lon[:, 1:], lat[:, 1:])
    dy = haversine(lon[:-1, :], lat[:-1, :], lon[1:, :], lat[1:, :])
    dx_full = np.concatenate([dx, dx[:, -1:,]], axis=1)
    dy_full = np.concatenate([dy, dy[-1:, :]], axis=0)
    ratio = dx_full / dy_full
    return dx_full, dy_full, ratio


def _unwrap_lon_for_plot(lon: np.ndarray) -> np.ndarray:
    """Unwrap longitudes row-wise to avoid dateline jumps in pcolormesh."""
    lon_unw = np.degrees(np.unwrap(np.radians(lon), axis=1))
    # recentre to median of first column to keep within a reasonable band
    center = np.median(lon_unw[:, 0])
    lon_wrap = ((lon_unw - center + 180.0) % 360.0) - 180.0 + center
    return lon_wrap


def plot_pmesh(field: np.ndarray, lon: np.ndarray, lat: np.ndarray, title: str, outfile: Path, cmap: str, units: str):
    fig = plt.figure(figsize=(10, 8))
    ax = plt.axes(projection=ccrs.PlateCarree())
    lon_plot = _unwrap_lon_for_plot(lon)
    pm = ax.pcolormesh(lon_plot, lat, field, cmap=cmap, shading="auto", transform=ccrs.PlateCarree())
    # set extent to the data bounds to reduce spurious fill
    ax.set_extent([lon_plot.min(), lon_plot.max(), lat.min(), lat.max()], crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.gridlines(draw_labels=True, linestyle=":")
    cbar = fig.colorbar(pm, ax=ax, orientation="vertical", pad=0.02)
    cbar.set_label(units)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(outfile, dpi=200)
    plt.close(fig)
    print(f"Wrote {outfile}")


def stitch(grids: list[tuple[Path, np.ndarray, np.ndarray, int]], tiles_x: int, tiles_y: int, nlon: int, nlat: int, col_major: bool):
    if len(grids) != tiles_x * tiles_y:
        raise ValueError(f"Expected {tiles_x*tiles_y} tiles, found {len(grids)}")
    full_lon = np.full((tiles_y * nlat, tiles_x * nlon), np.nan)
    full_lat = np.full_like(full_lon, np.nan)
    for path, lon_tile, lat_tile, rank in grids:
        if col_major:
            col = rank // tiles_y
            row = rank % tiles_y
        else:
            row = rank // tiles_x
            col = rank % tiles_x
        y0 = row * nlat
        x0 = col * nlon
        full_lon[y0:y0 + nlat, x0:x0 + nlon] = lon_tile
        full_lat[y0:y0 + nlat, x0:x0 + nlon] = lat_tile
    return full_lon, full_lat


def main():
    parser = argparse.ArgumentParser(description="Stitch real MGBF grid and plot dx/dy with pcolormesh")
    parser.add_argument("--pattern", default="mgbf_filtering_grid_latlon_*.txt", help="Glob for input tiles")
    parser.add_argument("--nlon", type=int, default=14, help="Core nlon per tile")
    parser.add_argument("--nlat", type=int, default=14, help="Core nlat per tile")
    parser.add_argument("--tiles-x", type=int, default=22, help="Tiles west->east")
    parser.add_argument("--tiles-y", type=int, default=22, help="Tiles south->north")
    parser.add_argument("--input-order", type=str, default="lonlat", choices=["lonlat", "latlon"], help="Storage order in files")
    parser.add_argument("--rank-order", type=str, default="col", choices=["col", "row"], help="Map rank to tiles: col=column-major (south->north fast)")
    parser.add_argument("--output-dir", type=Path, default=Path("dr-figures"), help="Output directory")
    parser.add_argument("--dump-stitched", action="store_true", help="Dump stitched lon/lat arrays and plots")
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

    lon_full, lat_full = stitch(grids, args.tiles_x, args.tiles_y, args.nlon, args.nlat, col_major=(args.rank_order == "col"))
    dx, dy, ratio = compute_dx_dy(lon_full, lat_full)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    plot_pmesh(dx / 1000.0, lon_full, lat_full, "Real dx (km)", args.output_dir / "real_dx_km_pmesh.png", cmap="magma", units="km")
    plot_pmesh(dy / 1000.0, lon_full, lat_full, "Real dy (km)", args.output_dir / "real_dy_km_pmesh.png", cmap="magma", units="km")
    plot_pmesh(ratio, lon_full, lat_full, "Real dx/dy", args.output_dir / "real_dx_over_dy_pmesh.png", cmap="coolwarm", units="ratio")

    if args.dump_stitched:
        np.savetxt(args.output_dir / "stitched_lon.txt", lon_full, fmt="%.8f")
        np.savetxt(args.output_dir / "stitched_lat.txt", lat_full, fmt="%.8f")
        plot_pmesh(lon_full, lon_full, lat_full, "Stitched lon", args.output_dir / "stitched_lon.png", cmap="coolwarm", units="deg")
        plot_pmesh(lat_full, lon_full, lat_full, "Stitched lat", args.output_dir / "stitched_lat.png", cmap="coolwarm", units="deg")


if __name__ == "__main__":
    main()
