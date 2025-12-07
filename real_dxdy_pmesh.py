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
        # lon varies fastest (west->east), then lat: reshape to (nlon, nlat) then transpose -> (nlat, nlon)
        lon = lon.reshape(nlon, nlat).T
        lat = lat.reshape(nlon, nlat).T
    else:  # latlon
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


def seam_diagnostics(lon_full: np.ndarray, lat_full: np.ndarray, nlon: int, nlat: int, tiles_x: int, tiles_y: int) -> dict:
    gaps_vert = []
    for c in range(1, tiles_x):
        lc = c * nlon - 1
        rc = c * nlon
        gaps_vert.append(haversine(lon_full[:, lc], lat_full[:, lc], lon_full[:, rc], lat_full[:, rc]))
    gaps_horz = []
    for r in range(1, tiles_y):
        br = r * nlat - 1
        tr = r * nlat
        gaps_horz.append(haversine(lon_full[br, :], lat_full[br, :], lon_full[tr, :], lat_full[tr, :]))

    def summarize(gaps):
        if not gaps:
            return {"mean": np.nan, "max": np.nan}
        merged = np.concatenate(gaps)
        return {"mean": float(np.nanmean(merged)), "max": float(np.nanmax(merged))}

    return {"vertical": summarize(gaps_vert), "horizontal": summarize(gaps_horz)}


def monotonic_checks(lon_full: np.ndarray, lat_full: np.ndarray, tol_deg: float = 1e-3) -> None:
    lon_unw = np.degrees(np.unwrap(np.radians(lon_full), axis=1))
    dlon = np.diff(lon_unw, axis=1)
    dlat = np.diff(lat_full, axis=0)
    min_dlon = float(np.nanmin(dlon))
    min_dlat = float(np.nanmin(dlat))
    print(f"Monotonic check: dlon min={min_dlon:.6f} deg, dlat min={min_dlat:.6f} deg")
    if min_dlon < -tol_deg:
        print(f"Warning: longitude decreases eastward (min step {min_dlon:.6f} deg < -{tol_deg} deg)")
    if min_dlat < -tol_deg:
        print(f"Warning: latitude decreases northward (min step {min_dlat:.6f} deg < -{tol_deg} deg)")


def enforce_orientation(lon_full: np.ndarray, lat_full: np.ndarray) -> tuple[np.ndarray, np.ndarray, dict]:
    """Flip rows/cols if needed so lon increases west->east and lat increases south->north."""
    info = {"flip_cols": False, "flip_rows": False}
    lon_unw = np.degrees(np.unwrap(np.radians(lon_full), axis=1))
    if np.nanmean(np.diff(lon_unw, axis=1)) < 0:
        lon_full = np.fliplr(lon_full)
        lat_full = np.fliplr(lat_full)
        info["flip_cols"] = True
    if np.nanmean(np.diff(lat_full, axis=0)) < 0:
        lon_full = np.flipud(lon_full)
        lat_full = np.flipud(lat_full)
        info["flip_rows"] = True
    return lon_full, lat_full, info


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
    parser.add_argument("--swap-axes", action="store_true", help="Transpose stitched lon/lat before computing dx/dy (for layout debugging)")
    parser.add_argument("--seam-threshold-km", type=float, default=5000.0, help="Warn if seam gaps exceed this (km)")
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
    if args.swap_axes:
        lon_full = lon_full.T
        lat_full = lat_full.T
    # Enforce orientation so lon increases eastward, lat northward
    lon_full, lat_full, flip_info = enforce_orientation(lon_full, lat_full)
    if flip_info["flip_cols"] or flip_info["flip_rows"]:
        print(f"Applied orientation fix: flip_cols={flip_info['flip_cols']}, flip_rows={flip_info['flip_rows']}")
    dx, dy, ratio = compute_dx_dy(lon_full, lat_full)

    seam_stats = seam_diagnostics(lon_full, lat_full, args.nlon, args.nlat, args.tiles_x, args.tiles_y)
    for lbl, stats in seam_stats.items():
        if stats["mean"] == stats["mean"]:
            print(f"{lbl.capitalize()} seam gaps (meters): mean={stats['mean']:.2f}, max={stats['max']:.2f}")
    if seam_stats["vertical"]["max"]/1000.0 > args.seam_threshold_km or seam_stats["horizontal"]["max"]/1000.0 > args.seam_threshold_km:
        print(f"Warning: seam gaps exceed {args.seam_threshold_km} km; layout may be mis-ordered.")
    monotonic_checks(lon_full, lat_full)

    # Diagnostics
    print(f"dx km mean/min/max: {np.nanmean(dx)/1000.0:.3f} / {np.nanmin(dx)/1000.0:.3f} / {np.nanmax(dx)/1000.0:.3f}")
    print(f"dy km mean/min/max: {np.nanmean(dy)/1000.0:.3f} / {np.nanmin(dy)/1000.0:.3f} / {np.nanmax(dy)/1000.0:.3f}")
    # Variation with latitude (per row)
    dx_row_std = np.nanstd(dx/1000.0, axis=1)
    dy_row_std = np.nanstd(dy/1000.0, axis=1)
    print(f"dx row-std km: mean={np.nanmean(dx_row_std):.3f}, min={np.nanmin(dx_row_std):.3f}, max={np.nanmax(dx_row_std):.3f}")
    print(f"dy row-std km: mean={np.nanmean(dy_row_std):.3f}, min={np.nanmin(dy_row_std):.3f}, max={np.nanmax(dy_row_std):.3f}")

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
