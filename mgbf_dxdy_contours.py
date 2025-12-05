"""
Build dx/dy contour plots for the MGBF filtering grid.

This version does **not** rely on an external grid description. It:
 1) reads all FieldImpl lat/lon dumps matching a glob (default:
    `mgbf_filtering*latlon_*.txt`),
 2) infers the structured grid by grouping points with nearly-equal
    latitudes into rows and sorting longitudes within each row,
 3) computes great-circle dx (east-west) and dy (north-south),
 4) writes three contour plots: dx, dy, and dx/dy.
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


def parse_fieldimpl_latlon(
    path: Path, core_nlon: int | None = None, core_nlat: int | None = None
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract lon/lat pairs from a FieldImpl dump.
    If core_nlon/core_nlat are provided, only the first (core_nlon * core_nlat) points
    are retained (assuming halo points follow). Returns 2D arrays (lon_grid, lat_grid).
    If parsing fails, both arrays are empty.
    """
    text = path.read_text()
    match = re.search(r"values:\s*\[(.*?)\]", text, re.DOTALL)
    if not match:
        return np.array([[]]), np.array([[]])

    numbers = [float(x) for x in re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", match.group(1))]
    if len(numbers) < 2:
        return np.array([[]]), np.array([[]])

    values = np.asarray(numbers)
    if values.size % 2:
        values = values[:-1]

    lon = values[0::2]
    lat = values[1::2]

    if core_nlon and core_nlat:
        core_expected = core_nlon * core_nlat
        if lon.size < core_expected or lat.size < core_expected:
            raise ValueError(
                f"{path}: not enough points for expected core grid: "
                f"found {lon.size} but need at least {core_expected} (including halo)"
            )
        lon = lon[:core_expected].reshape(core_nlat, core_nlon)
        lat = lat[:core_expected].reshape(core_nlat, core_nlon)
    else:
        # No shape expectations; return as 1D
        lon = lon
        lat = lat

    return lon, lat


def collect_grids(pattern: str, core_nlon: int | None, core_nlat: int | None):
    """Read all files that match the pattern and return per-file lon/lat grids."""
    grids = []
    files = sorted(Path(".").glob(pattern))
    for path in files:
        lon, lat = parse_fieldimpl_latlon(path, core_nlon=core_nlon, core_nlat=core_nlat)
        if lon.size and lat.size:
            grids.append((path, lon, lat))
        else:
            print(f"Warning: no lon/lat values parsed from {path}")
    return grids, files


def stitch_grids(
    grids: list[tuple[Path, np.ndarray, np.ndarray]],
    tile_nlon: int,
    tile_nlat: int,
    tiles_x: int,
    tiles_y: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Assemble individual tile grids into one stitched domain.
    tiles_x runs west->east (fast index), tiles_y runs south->north (rows).
    """
    expected_tiles = tiles_x * tiles_y
    if len(grids) != expected_tiles:
        raise ValueError(f"Expected {expected_tiles} tiles but found {len(grids)}")

    full_nlon = tiles_x * tile_nlon
    full_nlat = tiles_y * tile_nlat
    full_lon = np.full((full_nlat, full_nlon), np.nan)
    full_lat = np.full((full_nlat, full_nlon), np.nan)

    for idx, (_, lon_tile, lat_tile) in enumerate(grids):
        row = idx // tiles_x  # south to north
        col = idx % tiles_x   # west to east
        y0 = row * tile_nlat
        x0 = col * tile_nlon
        full_lon[y0:y0 + tile_nlat, x0:x0 + tile_nlon] = lon_tile
        full_lat[y0:y0 + tile_nlat, x0:x0 + tile_nlon] = lat_tile

    validate_latlon(full_lon.flatten(), full_lat.flatten())
    return full_lon, full_lat


def validate_latlon(lon: np.ndarray, lat: np.ndarray) -> None:
    """Validate lat/lon arrays and raise if abnormal values are present."""
    if lon.size != lat.size:
        raise ValueError(f"Mismatched lon/lat sizes: {lon.size} vs {lat.size}")
    if not np.all(np.isfinite(lon)) or not np.all(np.isfinite(lat)):
        bad = (~np.isfinite(lon)) | (~np.isfinite(lat))
        raise ValueError(f"Found non-finite lon/lat entries at indices: {np.where(bad)[0].tolist()}")
    if np.any((lat < -90) | (lat > 90)):
        bad = np.where((lat < -90) | (lat > 90))[0].tolist()
        raise ValueError(f"Found latitudes outside [-90, 90] at indices: {bad}")


def infer_structured_grid(
    lon: np.ndarray, lat: np.ndarray, tol: float = 1e-5
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build 2D lon/lat grids by grouping points with nearly-equal latitudes into rows
    and sorting longitudes within each row. All points are kept; rows are padded
    with NaN where column counts differ.
    """
    order = np.lexsort((lon, lat))  # primary sort by lat, then lon
    lon = lon[order]
    lat = lat[order]

    row_lons: list[list[float]] = []
    row_lats: list[list[float]] = []
    current_row_lon: list[float] = [lon[0]]
    current_row_lat: list[float] = [lat[0]]
    current_lat_ref = lat[0]

    for lo, la in zip(lon[1:], lat[1:]):
        if abs(la - current_lat_ref) <= tol:
            current_row_lon.append(lo)
            current_row_lat.append(la)
        else:
            row_lons.append(sorted(current_row_lon))
            mean_lat = sum(current_row_lat) / len(current_row_lat)
            row_lats.append([mean_lat] * len(current_row_lat))
            current_row_lon = [lo]
            current_row_lat = [la]
            current_lat_ref = la

    row_lons.append(sorted(current_row_lon))
    mean_lat = sum(current_row_lat) / len(current_row_lat)
    row_lats.append([mean_lat] * len(current_row_lat))

    if not row_lons:
        return np.array([[]]), np.array([[]])

    max_cols = max(len(r) for r in row_lons)
    padded_lon = np.full((len(row_lons), max_cols), np.nan)
    padded_lat = np.full((len(row_lats), max_cols), np.nan)
    for idx, (lons_row, lats_row) in enumerate(zip(row_lons, row_lats)):
        padded_lon[idx, : len(lons_row)] = lons_row
        padded_lat[idx, : len(lats_row)] = lats_row

    return padded_lon, padded_lat


def haversine(lon1: np.ndarray, lat1: np.ndarray, lon2: np.ndarray, lat2: np.ndarray) -> np.ndarray:
    """Great-circle distance (meters) between paired lon/lat arrays. Raises on invalid inputs."""
    if not (np.all(np.isfinite(lon1)) and np.all(np.isfinite(lat1)) and np.all(np.isfinite(lon2)) and np.all(np.isfinite(lat2))):
        raise ValueError("NaN/inf values passed to haversine; inputs must be fully finite.")

    lon1_rad = np.radians(lon1)
    lon2_rad = np.radians(lon2)
    lat1_rad = np.radians(lat1)
    lat2_rad = np.radians(lat2)
    dlon = lon2_rad - lon1_rad
    dlat = lat2_rad - lat1_rad
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1_rad) * np.cos(lat2_rad) * np.sin(dlon / 2) ** 2
    if np.any((a < -1e-12) | (a > 1 + 1e-12)):
        bad_idx = np.where((a < -1e-12) | (a > 1 + 1e-12))[0].tolist()
        raise ValueError(f"haversine encountered invalid 'a' values outside [0,1]; indices: {bad_idx}, min={a.min()}, max={a.max()}")
    a = np.clip(a, 0.0, 1.0)
    c = 2 * np.arcsin(np.sqrt(a))
    return 6_371_000.0 * c


def compute_dx_dy(lon_grid: np.ndarray, lat_grid: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute dx, dy (meters) and dx/dy ratio on the native grid."""
    dx = np.full_like(lon_grid, np.nan, dtype=float)
    dy = np.full_like(lon_grid, np.nan, dtype=float)

    # Forward differences for interior points
    dx[:, :-1] = haversine(lon_grid[:, :-1], lat_grid[:, :-1], lon_grid[:, 1:], lat_grid[:, 1:])
    dy[:-1, :] = haversine(lon_grid[:-1, :], lat_grid[:-1, :], lon_grid[1:, :], lat_grid[1:, :])
    # Boundary differences using inside neighbor
    dx[:, -1] = haversine(lon_grid[:, -1], lat_grid[:, -1], lon_grid[:, -2], lat_grid[:, -2])
    dy[-1, :] = haversine(lon_grid[-1, :], lat_grid[-1, :], lon_grid[-2, :], lat_grid[-2, :])

    # Validate all entries
    if np.any(~np.isfinite(dx)):
        bad = np.argwhere(~np.isfinite(dx)).tolist()
        raise ValueError(f"Non-finite dx encountered at indices (row,col): {bad[:10]}...")
    if np.any(~np.isfinite(dy)):
        bad = np.argwhere(~np.isfinite(dy)).tolist()
        raise ValueError(f"Non-finite dy encountered at indices (row,col): {bad[:10]}...")
    if np.any(dy == 0):
        zero_idx = np.argwhere(dy == 0).tolist()
        raise ValueError(f"Zero dy encountered at indices (row,col): {zero_idx[:10]}...")

    ratio = np.divide(dx, dy, out=np.full_like(dx, np.nan, dtype=float), where=np.isfinite(dx) & np.isfinite(dy) & (dy != 0))
    return dx, dy, ratio


def plot_field_grid(
    field: np.ndarray,
    lon: np.ndarray,
    lat: np.ndarray,
    title: str,
    outfile: Path,
    cmap: str = "viridis",
    units: str | None = None,
) -> None:
    """Contour plot on the stitched structured grid."""
    valid = np.isfinite(field) & np.isfinite(lon) & np.isfinite(lat)
    if not np.any(valid):
        print(f"Skipping plot for {outfile}: no valid data")
        return

    fig = plt.figure(figsize=(10, 8))
    ax = plt.axes(projection=ccrs.PlateCarree())
    cf = ax.contourf(lon, lat, field, 60, transform=ccrs.PlateCarree(), cmap=cmap)
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


def run(pattern: str, output_dir: Path, expected_nlon: int | None, expected_nlat: int | None, tiles_x: int, tiles_y: int) -> None:
    grids, files = collect_grids(pattern, core_nlon=expected_nlon, core_nlat=expected_nlat)
    if not grids:
        print(f"No lat/lon points found for pattern '{pattern}'. Files seen: {[str(f) for f in files]}")
        return

    for path, lon_grid, lat_grid in grids:
        validate_latlon(lon_grid.flatten(), lat_grid.flatten())
        if expected_nlat and expected_nlon:
            if lon_grid.shape != (expected_nlat, expected_nlon):
                print(
                    f"Warning: {path.name}: expected grid {expected_nlat}x{expected_nlon}, "
                    f"but got {lon_grid.shape[0]}x{lon_grid.shape[1]}"
                )
    if expected_nlon is None or expected_nlat is None:
        raise ValueError("expected_nlon and expected_nlat must be provided to stitch the domain.")

    lon_full, lat_full = stitch_grids(grids, expected_nlon, expected_nlat, tiles_x=tiles_x, tiles_y=tiles_y)
    dx_full, dy_full, ratio_full = compute_dx_dy(lon_full, lat_full)

    output_dir.mkdir(parents=True, exist_ok=True)
    plot_field_grid(dx_full / 1000.0, lon_full, lat_full, "dx (km) across domain", output_dir / "mgbf_dx_km.png", cmap="magma", units="km")
    plot_field_grid(dy_full / 1000.0, lon_full, lat_full, "dy (km) across domain", output_dir / "mgbf_dy_km.png", cmap="magma", units="km")
    plot_field_grid(ratio_full, lon_full, lat_full, "dx/dy across domain", output_dir / "mgbf_dx_over_dy.png", cmap="coolwarm")


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
    parser.add_argument(
        "--expected-nlon",
        type=int,
        default=14,
        help="Expected number of longitudes per subdomain grid (default: %(default)s)",
    )
    parser.add_argument(
        "--expected-nlat",
        type=int,
        default=14,
        help="Expected number of latitudes per subdomain grid (default: %(default)s)",
    )
    parser.add_argument(
        "--tiles-x",
        type=int,
        default=24,
        help="Number of subdomains west->east (default: %(default)s)",
    )
    parser.add_argument(
        "--tiles-y",
        type=int,
        default=24,
        help="Number of subdomains south->north (default: %(default)s)",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)
    run(args.pattern, args.output_dir, args.expected_nlon, args.expected_nlat, tiles_x=args.tiles_x, tiles_y=args.tiles_y)


if __name__ == "__main__":
    main()
