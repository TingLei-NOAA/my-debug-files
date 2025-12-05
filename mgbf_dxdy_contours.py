"""
Build dx/dy contour plots for the MGBF filtering grid.

The script reads the `FieldImpl[...] array=[values: ...]` text dumps that
match a glob such as `mgbf_filtering*latlon_*.txt`, remaps the scattered
points back onto the rotated structured grid described in the fv3jedi
configuration, computes great–circle distances in the x/y directions,
and writes three contour plots:
    1) dx over the full domain
    2) dy over the full domain
    3) dx/dy ratio over the full domain
"""

from __future__ import annotations

import argparse
import glob
import math
import re
from pathlib import Path
from typing import Iterable, Tuple

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np


# Grid description taken from the provided StructuredColumns space.
GRID_SPEC = {
    "nx": 308,
    "ny": 308,
    "x_start": -60.0310278618037,
    "x_end": 60.0310278618037,
    "y_start": -36.6574659140194,
    "y_end": 36.6574659140194,
    "halo": 1,
    "north_pole_lon": 67.4999923706,
    "north_pole_lat": 34.9999966097,
}


def parse_fieldimpl_latlon(path: Path) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract lon/lat pairs from a FieldImpl dump.
    Returns two 1D arrays (lon, lat). If parsing fails, both arrays are empty.
    """
    text = path.read_text()
    match = re.search(r"values:\s*\[(.*?)\]", text, re.DOTALL)
    if not match:
        return np.array([]), np.array([])

    numbers = [float(x) for x in re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", match.group(1))]
    if len(numbers) < 2:
        return np.array([]), np.array([])

    values = np.asarray(numbers)
    if values.size % 2:
        values = values[:-1]

    lon = values[0::2]
    lat = values[1::2]
    return lon, lat


def collect_latlon(pattern: str) -> Tuple[np.ndarray, np.ndarray, list[Path]]:
    """Read all files that match the pattern and concatenate lon/lat arrays."""
    lons: list[np.ndarray] = []
    lats: list[np.ndarray] = []
    files = sorted(Path(".").glob(pattern))
    for path in files:
        lon, lat = parse_fieldimpl_latlon(path)
        if lon.size and lat.size:
            lons.append(lon)
            lats.append(lat)
        else:
            print(f"Warning: no lon/lat values parsed from {path}")

    if not lons:
        return np.array([]), np.array([]), files

    return np.concatenate(lons), np.concatenate(lats), files


def build_rotated_axes(spec: dict) -> Tuple[np.ndarray, np.ndarray, float, float]:
    """
    Construct the rotated-grid x/y axes (including halo) and their spacings.
    Halo cells are created by extending the start/end by one native spacing.
    """
    dx = (spec["x_end"] - spec["x_start"]) / (spec["nx"] - 1)
    dy = (spec["y_end"] - spec["y_start"]) / (spec["ny"] - 1)
    nx = spec["nx"] + 2 * spec["halo"]
    ny = spec["ny"] + 2 * spec["halo"]

    x0 = spec["x_start"] - spec["halo"] * dx
    y0 = spec["y_start"] - spec["halo"] * dy
    x = x0 + np.arange(nx) * dx
    y = y0 + np.arange(ny) * dy
    return x, y, dx, dy


def map_points_to_grid(
    lon: np.ndarray, lat: np.ndarray, spec: dict, rotated_crs: ccrs.CRS
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Map scattered geodetic lon/lat pairs onto the expected rotated grid cells.
    Returns (lon_grid, lat_grid, expected_lon_grid, expected_lat_grid).
    """
    x_axis, y_axis, dx_rot, dy_rot = build_rotated_axes(spec)
    nx = x_axis.size
    ny = y_axis.size
    lon_grid = np.full((ny, nx), np.nan)
    lat_grid = np.full((ny, nx), np.nan)

    # Convert input geodetic coords to rotated coords, then snap to the nearest cell.
    rotated_points = rotated_crs.transform_points(ccrs.Geodetic(), lon, lat)
    rot_lon = rotated_points[:, 0]
    rot_lat = rotated_points[:, 1]

    x_idx = np.rint((rot_lon - x_axis[0]) / dx_rot).astype(int)
    y_idx = np.rint((rot_lat - y_axis[0]) / dy_rot).astype(int)

    rejected = 0
    collisions = 0
    for i, j, g_lon, g_lat in zip(x_idx, y_idx, lon, lat):
        if 0 <= i < nx and 0 <= j < ny:
            if math.isnan(lon_grid[j, i]):
                lon_grid[j, i] = g_lon
                lat_grid[j, i] = g_lat
            else:
                collisions += 1
        else:
            rejected += 1

    if rejected:
        print(f"Warning: {rejected} points fell outside the expected grid bounds")
    if collisions:
        print(f"Warning: {collisions} grid cells received more than one point")

    xx, yy = np.meshgrid(x_axis, y_axis)
    expected_geo = ccrs.Geodetic().transform_points(rotated_crs, xx, yy)
    expected_lon = expected_geo[..., 0]
    expected_lat = expected_geo[..., 1]
    return lon_grid, lat_grid, expected_lon, expected_lat


def haversine(lon1: np.ndarray, lat1: np.ndarray, lon2: np.ndarray, lat2: np.ndarray) -> np.ndarray:
    """Great-circle distance (meters) between paired lon/lat arrays."""
    lon1_rad = np.radians(lon1)
    lon2_rad = np.radians(lon2)
    lat1_rad = np.radians(lat1)
    lat2_rad = np.radians(lat2)
    dlon = lon2_rad - lon1_rad
    dlat = lat2_rad - lat1_rad
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1_rad) * np.cos(lat2_rad) * np.sin(dlon / 2) ** 2
    c = 2 * np.arcsin(np.sqrt(a))
    return 6_371_000.0 * c


def compute_dx_dy(lon_grid: np.ndarray, lat_grid: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute dx, dy (meters) and dx/dy ratio on the native grid."""
    dx = np.full_like(lon_grid, np.nan, dtype=float)
    dy = np.full_like(lon_grid, np.nan, dtype=float)

    dx[:, :-1] = haversine(lon_grid[:, :-1], lat_grid[:, :-1], lon_grid[:, 1:], lat_grid[:, 1:])
    dy[:-1, :] = haversine(lon_grid[:-1, :], lat_grid[:-1, :], lon_grid[1:, :], lat_grid[1:, :])
    ratio = dx / dy
    return dx, dy, ratio


def plot_field(
    field: np.ndarray,
    lon: np.ndarray,
    lat: np.ndarray,
    title: str,
    outfile: Path,
    cmap: str = "viridis",
    units: str | None = None,
) -> None:
    """Contour-plot a 2D field on a PlateCarree map."""
    fig = plt.figure(figsize=(10, 8))
    ax = plt.axes(projection=ccrs.PlateCarree())
    valid = np.isfinite(field) & np.isfinite(lon) & np.isfinite(lat)
    if not np.any(valid):
        print(f"Skipping plot for {outfile}: no valid data")
        return

    cf = ax.contourf(lon, lat, field, 40, transform=ccrs.PlateCarree(), cmap=cmap)
    ax.coastlines()
    ax.gridlines(draw_labels=True, linestyle=":")
    cbar = fig.colorbar(cf, ax=ax, orientation="vertical", pad=0.02)
    label = units if units else ""
    if label:
        cbar.set_label(label)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(outfile, dpi=200)
    plt.close(fig)
    print(f"Wrote {outfile}")


def print_alignment_stats(lon_grid: np.ndarray, lat_grid: np.ndarray, exp_lon: np.ndarray, exp_lat: np.ndarray) -> None:
    """Report how well the ingested points line up with the configured grid."""
    valid = np.isfinite(lon_grid) & np.isfinite(lat_grid)
    coverage = valid.sum()
    print(f"Mapped {coverage} points onto the grid")
    if not coverage:
        return

    # Normalize longitudes before differencing.
    lon_diff = ((lon_grid[valid] - exp_lon[valid] + 180.0) % 360.0) - 180.0
    lat_diff = lat_grid[valid] - exp_lat[valid]
    print(
        f"Alignment vs rotated grid (lon/lat degrees): "
        f"mean abs [{np.mean(np.abs(lon_diff)):.4f}, {np.mean(np.abs(lat_diff)):.4f}] "
        f"max abs [{np.max(np.abs(lon_diff)):.4f}, {np.max(np.abs(lat_diff)):.4f}]"
    )


def run(pattern: str, output_dir: Path) -> None:
    lon, lat, files = collect_latlon(pattern)
    if not lon.size:
        print(f"No lat/lon points found for pattern '{pattern}'. Files seen: {[str(f) for f in files]}")
        return

    print(f"Loaded {lon.size} points from {len(files)} files")
    rotated_crs = ccrs.RotatedPole(
        pole_longitude=GRID_SPEC["north_pole_lon"], pole_latitude=GRID_SPEC["north_pole_lat"]
    )
    lon_grid, lat_grid, exp_lon, exp_lat = map_points_to_grid(lon, lat, GRID_SPEC, rotated_crs)
    print_alignment_stats(lon_grid, lat_grid, exp_lon, exp_lat)

    dx, dy, ratio = compute_dx_dy(lon_grid, lat_grid)
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_field(dx / 1000.0, lon_grid, lat_grid, "dx (km)", output_dir / "mgbf_dx_km.png", cmap="magma", units="km")
    plot_field(dy / 1000.0, lon_grid, lat_grid, "dy (km)", output_dir / "mgbf_dy_km.png", cmap="magma", units="km")
    plot_field(ratio, lon_grid, lat_grid, "dx/dy", output_dir / "mgbf_dx_over_dy.png", cmap="coolwarm")


def main(argv: Iterable[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Plot dx, dy, and dx/dy for mgbf filtering rotated grid")
    parser.add_argument(
        "--pattern",
        default="mgbf_filtering*latlon_*.txt",
        help="Glob for FieldImpl lat/lon dump files (default: %(default)s)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("figures"),
        help="Directory to write plots (default: %(default)s)",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)
    run(args.pattern, args.output_dir)


if __name__ == "__main__":
    main()
