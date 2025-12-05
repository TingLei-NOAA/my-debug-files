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
    path: Path,
    core_nlon: int | None = None,
    core_nlat: int | None = None,
    input_order: str = "latlon",
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract lon/lat pairs from a FieldImpl dump.
    If core_nlon/core_nlat are provided, only the first (core_nlon * core_nlat) points
    are retained (assuming halo points follow). Returns 2D arrays (lon_grid, lat_grid).
    input_order controls how the flat core is reshaped:
      - "latlon": reshape to (nlat, nlon) directly (lat-major rows)
      - "lonlat": reshape to (nlon, nlat) then transpose to (nlat, nlon)
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
        lon = lon[:core_expected]
        lat = lat[:core_expected]
        if input_order == "latlon":
            lon = lon.reshape(core_nlat, core_nlon)
            lat = lat.reshape(core_nlat, core_nlon)
        elif input_order == "lonlat":
            lon = lon.reshape(core_nlon, core_nlat).T
            lat = lat.reshape(core_nlon, core_nlat).T
        else:
            raise ValueError(f"Unsupported input_order '{input_order}'. Use 'latlon' or 'lonlat'.")
    else:
        # No shape expectations; return as 1D
        lon = lon
        lat = lat

    return lon, lat


def collect_grids(
    pattern: str,
    core_nlon: int | None,
    core_nlat: int | None,
    input_order: str,
):
    """Read all files that match the pattern and return per-file lon/lat grids plus optional indices."""
    grids = []
    files = sorted(Path(".").glob(pattern))
    for path in files:
        lon, lat = parse_fieldimpl_latlon(path, core_nlon=core_nlon, core_nlat=core_nlat, input_order=input_order)
        if lon.size and lat.size:
            # Parse MPI rank from filename suffix (expects ..._<rank>.txt)
            m = re.search(r"_([0-9]+)\\.txt$", path.name)
            if not m:
                raise ValueError(f"Could not extract rank index from {path.name} (expected ..._<rank>.txt)")
            rank = int(m.group(1))
            grids.append((path, lon, lat, rank, None))
        else:
            print(f"Warning: no lon/lat values parsed from {path}")
    return grids, files


def stitch_grids(
    grids: list[tuple[Path, np.ndarray, np.ndarray, int | None, int | None]],
    tile_nlon: int,
    tile_nlat: int,
    tiles_x: int,
    tiles_y: int,
    rank_order: str = "row",
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

    resolved = []
    for path, lon_tile, lat_tile, rank, _ in grids:
        if rank is None:
            raise ValueError("Missing rank for tile placement.")
        if rank < 0 or rank >= tiles_x * tiles_y:
            raise ValueError(f"Rank {rank} out of range for tiles {tiles_x}x{tiles_y}")
        if rank_order.lower() == "row":
            row = rank // tiles_x
            col = rank % tiles_x
        else:  # col-major
            col = rank // tiles_y
            row = rank % tiles_y
        resolved.append((path, lon_tile, lat_tile, col, row))

    for _, lon_tile, lat_tile, col, row in resolved:
        y0 = row * tile_nlat
        x0 = col * tile_nlon
        full_lon[y0:y0 + tile_nlat, x0:x0 + tile_nlon] = lon_tile
        full_lat[y0:y0 + tile_nlat, x0:x0 + tile_nlon] = lat_tile

    validate_latlon(full_lon.flatten(), full_lat.flatten())
    return full_lon, full_lat


def infer_indices_from_centers(
    grids: list[tuple[Path, np.ndarray, np.ndarray, int | None, int | None]],
    tiles_x: int,
    tiles_y: int,
    tol_decimals: int = 6,
) -> list[tuple[Path, np.ndarray, np.ndarray, int, int]]:
    """Assign (tx, ty) based on tile centroid lon/lat ordering (west->east, south->north)."""
    lon_centers = []
    lat_centers = []
    for _, lon_tile, lat_tile, _, _ in grids:
        lon_centers.append(np.nanmean(lon_tile))
        lat_centers.append(np.nanmean(lat_tile))
    lon_centers = np.array(lon_centers)
    lat_centers = np.array(lat_centers)

    # Unwrap longitudes to avoid dateline jumps, then round to clusters.
    lon_unwrapped = np.degrees(np.unwrap(np.radians(lon_centers)))
    uniq_lon = np.unique(np.round(lon_unwrapped, tol_decimals))
    uniq_lat = np.unique(np.round(lat_centers, tol_decimals))
    if len(uniq_lon) != tiles_x or len(uniq_lat) != tiles_y:
        raise ValueError(
            f"Auto layout: unique lon centers={len(uniq_lon)} (need {tiles_x}), "
            f"unique lat centers={len(uniq_lat)} (need {tiles_y})."
        )
    uniq_lon_sorted = np.sort(uniq_lon)
    uniq_lat_sorted = np.sort(uniq_lat)

    remapped = []
    for (path, lon_tile, lat_tile, _, _) , lon_c, lat_c in zip(grids, lon_unwrapped, lat_centers):
        lon_key = np.round(lon_c, tol_decimals)
        lat_key = np.round(lat_c, tol_decimals)
        tx = int(np.where(uniq_lon_sorted == lon_key)[0][0])
        ty = int(np.where(uniq_lat_sorted == lat_key)[0][0])
        remapped.append((path, lon_tile, lat_tile, tx, ty))
    # Note: ty grows northward because lat_key increases northward in sorting
    return remapped


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


def seam_diagnostics(lon_full: np.ndarray, lat_full: np.ndarray, tile_nlon: int, tile_nlat: int, tiles_x: int, tiles_y: int) -> dict:
    """Compute seam gaps between tiles to catch ordering/placement issues."""
    # Vertical seams (between columns)
    gaps_vert = []
    for c in range(1, tiles_x):
        left_col = c * tile_nlon - 1
        right_col = c * tile_nlon
        gaps_vert.append(haversine(lon_full[:, left_col], lat_full[:, left_col], lon_full[:, right_col], lat_full[:, right_col]))
    # Horizontal seams (between rows)
    gaps_horz = []
    for r in range(1, tiles_y):
        bottom_row = r * tile_nlat - 1
        top_row = r * tile_nlat
        gaps_horz.append(haversine(lon_full[bottom_row, :], lat_full[bottom_row, :], lon_full[top_row, :], lat_full[top_row, :]))

    def summarize(gaps):
        if not gaps:
            return {"mean": np.nan, "max": np.nan}
        merged = np.concatenate(gaps)
        return {"mean": float(np.nanmean(merged)), "max": float(np.nanmax(merged))}

    vert = summarize(gaps_vert)
    horz = summarize(gaps_horz)
    return {"vertical": vert, "horizontal": horz}


def monotonic_checks(lon_full: np.ndarray, lat_full: np.ndarray, tol_deg: float) -> None:
    """
    Ensure longitude increases west->east and latitude increases south->north.
    Longitudes are unwrapped along x to handle dateline crossings.
    """
    # Unwrap longitudes along x-direction to avoid artificial jumps at -180/180
    lon_unwrapped = np.degrees(np.unwrap(np.radians(lon_full), axis=1))
    dlon = np.diff(lon_unwrapped, axis=1)
    dlat = np.diff(lat_full, axis=0)
    min_dlon = float(np.nanmin(dlon))
    max_dlon = float(np.nanmax(dlon))
    min_dlat = float(np.nanmin(dlat))
    max_dlat = float(np.nanmax(dlat))
    print(f"Monotonic check: lon step min/max = {min_dlon:.6f}/{max_dlon:.6f} deg; lat step min/max = {min_dlat:.6f}/{max_dlat:.6f} deg")
    if min_dlon < -tol_deg:
        raise ValueError(f"Longitude decreases eastward (min step {min_dlon:.6f} deg < -{tol_deg} deg). Check tile ordering/layout.")
    if min_dlat < -tol_deg:
        raise ValueError(f"Latitude decreases northward (min step {min_dlat:.6f} deg < -{tol_deg} deg). Check tile ordering/layout.")


def plot_field_grid(
    field: np.ndarray,
    lon: np.ndarray,
    lat: np.ndarray,
    title: str,
    outfile: Path,
    cmap: str = "viridis",
    units: str | None = None,
    cbar_label: str | None = None,
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
    label = cbar_label if cbar_label is not None else (units if units else "")
    cbar.set_label(label)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(outfile, dpi=200)
    plt.close(fig)
    print(f"Wrote {outfile}")


def run(
    pattern: str,
    output_dir: Path,
    expected_nlon: int | None,
    expected_nlat: int | None,
    tiles_x: int,
    tiles_y: int,
    input_order: str,
    seam_threshold_km: float,
    monotonic_tol_deg: float,
    auto_layout: bool,
    auto_layout_decimals: int,
    rank_order: str,
    dump_stitched: bool,
    plot_subdomains: bool,
) -> None:
    grids, files = collect_grids(
        pattern,
        core_nlon=expected_nlon,
        core_nlat=expected_nlat,
        input_order=input_order,
    )
    if not grids:
        print(f"No lat/lon points found for pattern '{pattern}'. Files seen: {[str(f) for f in files]}")
        return

    if auto_layout:
        grids = infer_indices_from_centers(grids, tiles_x=tiles_x, tiles_y=tiles_y, tol_decimals=auto_layout_decimals)

    for path, lon_grid, lat_grid, tx, ty in grids:
        validate_latlon(lon_grid.flatten(), lat_grid.flatten())
        if expected_nlat and expected_nlon:
            if lon_grid.shape != (expected_nlat, expected_nlon):
                print(
                    f"Warning: {path.name}: expected grid {expected_nlat}x{expected_nlon}, "
                    f"but got {lon_grid.shape[0]}x{lon_grid.shape[1]}"
                )
        dx_tile, dy_tile, ratio_tile = compute_dx_dy(lon_grid, lat_grid)
        max_dx_km = np.nanmax(dx_tile) / 1000.0
        if max_dx_km > 1000.0:
            raise ValueError(f"{path.name}: dx exceeds 1000 km (max={max_dx_km:.2f} km)")

        if plot_subdomains:
            stem = path.stem
            suffix = ""
            if tx is not None and ty is not None:
                suffix = f"_tx{tx}_ty{ty}"
            plot_field_grid(
                dx_tile / 1000.0,
                lon_grid,
                lat_grid,
                f"{stem}{suffix} dx (km)",
                output_dir / f"{stem}{suffix}_dx_km.png",
                cmap="magma",
                units="km",
            )
            plot_field_grid(
                dy_tile / 1000.0,
                lon_grid,
                lat_grid,
                f"{stem}{suffix} dy (km)",
                output_dir / f"{stem}{suffix}_dy_km.png",
                cmap="magma",
                units="km",
            )
            plot_field_grid(
                ratio_tile,
                lon_grid,
                lat_grid,
                f"{stem}{suffix} dx/dy",
                output_dir / f"{stem}{suffix}_dx_over_dy.png",
                cmap="coolwarm",
            )
    if expected_nlon is None or expected_nlat is None:
        raise ValueError("expected_nlon and expected_nlat must be provided to stitch the domain.")

    lon_full, lat_full = stitch_grids(
        grids, expected_nlon, expected_nlat, tiles_x=tiles_x, tiles_y=tiles_y, rank_order=rank_order
    )
    seam_stats = seam_diagnostics(lon_full, lat_full, expected_nlon, expected_nlat, tiles_x, tiles_y)
    for lbl, stats in seam_stats.items():
        if stats["mean"] == stats["mean"]:  # not NaN
            print(f"{lbl.capitalize()} seam gaps (meters): mean={stats['mean']:.2f}, max={stats['max']:.2f}")
    if (
        seam_stats["vertical"]["max"] / 1000.0 > seam_threshold_km
        or seam_stats["horizontal"]["max"] / 1000.0 > seam_threshold_km
    ):
        raise ValueError(
            f"Seam gaps exceed {seam_threshold_km} km (vertical max={seam_stats['vertical']['max']/1000.0:.2f} km, "
            f"horizontal max={seam_stats['horizontal']['max']/1000.0:.2f} km). Check tile ordering/indices."
        )

    monotonic_checks(lon_full, lat_full, tol_deg=monotonic_tol_deg)

    if dump_stitched:
        output_dir.mkdir(parents=True, exist_ok=True)
        np.savetxt(output_dir / "stitched_lon.txt", lon_full, fmt="%.8f")
        np.savetxt(output_dir / "stitched_lat.txt", lat_full, fmt="%.8f")
        print(f"Wrote stitched lon/lat to {output_dir / 'stitched_lon.txt'} and {output_dir / 'stitched_lat.txt'}")
        # Quick plots of stitched lon/lat surfaces
        plot_field_grid(
            lon_full,
            lon_full,
            lat_full,
            "Stitched Lon",
            output_dir / "stitched_lon.png",
            cmap="coolwarm",
            units="deg",
            cbar_label="Longitude (deg)",
        )
        plot_field_grid(
            lat_full,
            lon_full,
            lat_full,
            "Stitched Lat",
            output_dir / "stitched_lat.png",
            cmap="coolwarm",
            units="deg",
            cbar_label="Latitude (deg)",
        )

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
        default=Path("dr-figures"),
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
        default=22,
        help="Number of subdomains west->east (default: %(default)s)",
    )
    parser.add_argument(
        "--tiles-y",
        type=int,
        default=22,
        help="Number of subdomains south->north (default: %(default)s)",
    )
    parser.add_argument(
        "--input-order",
        type=str,
        default="lonlat",
        choices=["latlon", "lonlat"],
        help="Order of the 2D core in the file: 'latlon' (rows = lat, cols = lon) or 'lonlat' (rows = lon, cols = lat).",
    )
    parser.add_argument(
        "--rank-order",
        type=str,
        default="row",
        choices=["row", "col"],
        help="If a single rank is parsed, map it row-major ('row': west->east then south->north) or column-major ('col').",
    )
    parser.add_argument(
        "--auto-layout",
        action="store_true",
        help="(Deprecated) Infer tile positions from centroid lon/lat; rank-based placement is recommended.",
    )
    parser.add_argument(
        "--auto-layout-decimals",
        type=int,
        default=3,
        help="(Deprecated) Rounding decimals for centroid clustering when using --auto-layout",
    )
    parser.add_argument(
        "--seam-threshold-km",
        type=float,
        default=1000.0,
        help="Raise an error if any seam gap exceeds this threshold (km). Default: %(default)s",
    )
    parser.add_argument(
        "--monotonic-tol-deg",
        type=float,
        default=1e-4,
        help="Tolerance (degrees) for detecting decreases in lon/lat across stitched grid (default: %(default)s)",
    )
    parser.add_argument(
        "--dump-stitched",
        action="store_true",
        help="Write stitched lon/lat grids to text files for inspection.",
    )
    parser.add_argument(
        "--plot-subdomains",
        action="store_true",
        help="Also plot dx/dy/dx_over_dy for each subdomain separately.",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)
    run(
        args.pattern,
        args.output_dir,
        args.expected_nlon,
        args.expected_nlat,
        tiles_x=args.tiles_x,
        tiles_y=args.tiles_y,
        input_order=args.input_order,
        seam_threshold_km=args.seam_threshold_km,
        monotonic_tol_deg=args.monotonic_tol_deg,
        auto_layout=args.auto_layout,
        auto_layout_decimals=args.auto_layout_decimals,
        rank_order=args.rank_order,
        dump_stitched=args.dump_stitched,
        plot_subdomains=args.plot_subdomains,
    )


if __name__ == "__main__":
    main()
