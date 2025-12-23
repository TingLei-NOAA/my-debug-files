"""
Match each filtering subdomain (tile) center to the closest FV3 model grid point.

Inputs:
  - Filtering grid dumps: mgbf_filtering_grid_latlon_<rank>.txt (first 14x14 points are the core)
  - FV3 model grid: fv3_grid_spec (prefers grid_xt/grid_yt when present; falls back to grid_latt/grid_lont)

Outputs:
  - Text file with one line per filtering subdomain (tile):
      label rank row col center_filt_i center_filt_j center_filt_lon center_filt_lat
      model_i model_j model_lon model_lat dist_km
    where center_filt_i/center_filt_j are the stitched filtering-grid indices (0-based, x then y),
    and model_i/model_j are the FV3 T-grid indices (grid_xt, grid_yt).
  - Optional plot of selected model grid points vs filtering centers.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import netCDF4 as nc

try:
    from scipy.spatial import cKDTree
except ImportError as e:  # pragma: no cover
    raise SystemExit("scipy is required for nearest-neighbor search (cKDTree)") from e


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


def read_model_grid(fv3_spec: Path) -> tuple[np.ndarray, np.ndarray]:
    with nc.Dataset(fv3_spec) as ds:
        if "grid_xt" in ds.variables and "grid_yt" in ds.variables:
            lon_1d = np.array(ds.variables["grid_xt"][:])
            lat_1d = np.array(ds.variables["grid_yt"][:])
            lon, lat = np.meshgrid(lon_1d, lat_1d)
        else:
            lon = np.array(ds.variables["grid_lont"][:])
            lat = np.array(ds.variables["grid_latt"][:])
    # Mask missing values
    mask = ~np.isfinite(lon) | ~np.isfinite(lat)
    lon[mask] = np.nan
    lat[mask] = np.nan
    return lon, lat


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


def lonlat_to_xyz(lon_deg: np.ndarray, lat_deg: np.ndarray) -> np.ndarray:
    lon_rad = np.radians(lon_deg)
    lat_rad = np.radians(lat_deg)
    cos_lat = np.cos(lat_rad)
    x = cos_lat * np.cos(lon_rad)
    y = cos_lat * np.sin(lon_rad)
    z = np.sin(lat_rad)
    return np.stack([x, y, z], axis=-1)


def normalize_lon_180(lon_deg: np.ndarray) -> np.ndarray:
    return (lon_deg + 180.0) % 360.0 - 180.0


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


def main():
    parser = argparse.ArgumentParser(description="Map filtering grid points to nearest FV3 model grid points.")
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
    parser.add_argument("--output", type=Path, default=Path("filter_to_model_map.txt"), help="Output text file")
    parser.add_argument(
        "--output-simple",
        type=Path,
        default=Path("filter_to_model_map_indices.txt"),
        help="Simplified output: index + filtering/model i/j",
    )
    parser.add_argument("--plot", type=Path, default=None, help="Output plot for stitched filtering grid")
    parser.add_argument(
        "--plot-selected",
        type=Path,
        default=Path("selected_model_grids_4panel.png"),
        help="4-panel plot of filtering centers and their selected model grid points",
    )
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
    ny_filt, nx_filt = lon_filt.shape
    print(f"Filtering grid stitched: shape (ny, nx)=({ny_filt},{nx_filt})")

    if args.plot is not None:
        try:
            import matplotlib.pyplot as plt
        except ImportError as e:  # pragma: no cover
            raise SystemExit("matplotlib is required for plotting the filtering grid") from e

        fig, ax = plt.subplots(figsize=(8, 6), dpi=150)
        ax.scatter(lon_filt.ravel(), lat_filt.ravel(), s=1, c="tab:blue", alpha=0.6, rasterized=True)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_title("Stitched filtering grid (lat/lon)")
        ax.grid(True, linewidth=0.3, alpha=0.4)
        fig.tight_layout()
        fig.savefig(args.plot)
        plt.close(fig)
        print(f"Wrote {args.plot}")

    lon_model, lat_model = read_model_grid(args.fv3_spec)
    lon_model = normalize_lon_180(lon_model)
    ny_m, nx_m = lon_model.shape
    print(f"Model grid (T) shape (ny, nx)=({ny_m},{nx_m})")

    # Build KDTree on model grid in 3D cartesian
    xyz_model = lonlat_to_xyz(lon_model.ravel(), lat_model.ravel())
    valid_mask = np.isfinite(xyz_model).all(axis=1)
    if not valid_mask.all():
        xyz_model = xyz_model[valid_mask]
        lon_model_flat = lon_model.ravel()[valid_mask]
        lat_model_flat = lat_model.ravel()[valid_mask]
        valid_indices = np.nonzero(valid_mask)[0]
    else:
        lon_model_flat = lon_model.ravel()
        lat_model_flat = lat_model.ravel()
        valid_indices = np.arange(xyz_model.shape[0])
    tree = cKDTree(xyz_model)

    center_lon = []
    center_lat = []
    center_filt_i = []
    center_filt_j = []
    center_row = []
    center_col = []
    center_rank = []
    center_label = []
    center_label_with_index = []
    center_y = (args.nlat - 1) // 2
    center_x = (args.nlon - 1) // 2

    for path, lon_tile, lat_tile, rank in grids:
        row, col = rank_to_row_col(rank, args.tiles_x, args.tiles_y, args.rank_order)
        y = row * args.nlat + center_y
        x = col * args.nlon + center_x
        center_lon.append(lon_filt[y, x])
        center_lat.append(lat_filt[y, x])
        center_filt_i.append(x)
        center_filt_j.append(y)
        center_row.append(row)
        center_col.append(col)
        center_rank.append(rank)
        label = f"subdomain_r{row:02d}_c{col:02d}"
        center_label.append(label)
        center_label_with_index.append(f"{label}_{rank}")

    center_lon = np.asarray(center_lon)
    center_lat = np.asarray(center_lat)
    center_filt_i = np.asarray(center_filt_i)
    center_filt_j = np.asarray(center_filt_j)
    center_row = np.asarray(center_row)
    center_col = np.asarray(center_col)
    center_rank = np.asarray(center_rank)
    center_label = np.asarray(center_label)
    center_label_with_index = np.asarray(center_label_with_index)

    order = np.argsort(center_rank)
    center_lon = center_lon[order]
    center_lat = center_lat[order]
    center_filt_i = center_filt_i[order]
    center_filt_j = center_filt_j[order]
    center_row = center_row[order]
    center_col = center_col[order]
    center_rank = center_rank[order]
    center_label = center_label[order]
    center_label_with_index = center_label_with_index[order]

    xyz_centers = lonlat_to_xyz(center_lon, center_lat)
    dist_cart, idx = tree.query(xyz_centers, k=1)
    model_flat_idx = valid_indices[idx]
    model_j, model_i = np.unravel_index(model_flat_idx, (ny_m, nx_m))

    d_gc = haversine(center_lon, center_lat, lon_model_flat[idx], lat_model_flat[idx]) / 1000.0

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as f:
        f.write("# NOTE: row/col/center_filt_i/center_filt_j/model_i/model_j are 1-based (added +1).\n")
        f.write(
            "# label_with_index rank row col center_filt_i center_filt_j center_filt_lon center_filt_lat "
            "model_i model_j model_lon model_lat dist_km\n"
        )
        for label, rank, row, col, fi, fj, flon, flat, mi, mj, mlon, mlat, dkm in zip(
            center_label_with_index, center_rank, center_row, center_col, center_filt_i, center_filt_j,
            center_lon, center_lat, model_i, model_j, lon_model_flat[idx], lat_model_flat[idx], d_gc
        ):
            f.write(
                f"{label} {rank} {row + 1} {col + 1} {fi + 1} {fj + 1} {flon:.6f} {flat:.6f} "
                f"{mi + 1} {mj + 1} {mlon:.6f} {mlat:.6f} {dkm:.3f}\n"
            )
    print(f"Wrote {args.output}")

    args.output_simple.parent.mkdir(parents=True, exist_ok=True)
    with args.output_simple.open("w") as f:
        f.write("# NOTE: index/center_filt_i/center_filt_j/model_i/model_j are 1-based (added +1).\n")
        f.write("# index center_filt_i center_filt_j model_i model_j\n")
        for rank, fi, fj, mi, mj in zip(center_rank, center_filt_i, center_filt_j, model_i, model_j):
            f.write(f"{rank} {fi + 1} {fj + 1} {mi + 1} {mj + 1}\n")
    print(f"Wrote {args.output_simple}")
    for label, rank, fi, fj, mi, mj in zip(center_label_with_index, center_rank, center_filt_i, center_filt_j, model_i, model_j):
        print(f"{label} rank={rank} center_filt=({fi + 1},{fj + 1}) model=({mi + 1},{mj + 1})")

    if args.plot_selected is not None:
        try:
            import matplotlib.pyplot as plt
        except ImportError as e:  # pragma: no cover
            raise SystemExit("matplotlib is required for plotting selected model grids") from e

        fig, axes = plt.subplots(2, 2, figsize=(10, 8), dpi=150, sharex=True, sharey=True)
        axes = axes.ravel()
        model_i_arr = model_i.astype(float)
        model_j_arr = model_j.astype(float)
        i_split = np.nanmedian(model_i_arr)
        j_split = np.nanmedian(model_j_arr)
        subsets = [
            ("NW", (model_i_arr < i_split) & (model_j_arr >= j_split)),
            ("NE", (model_i_arr >= i_split) & (model_j_arr >= j_split)),
            ("SW", (model_i_arr < i_split) & (model_j_arr < j_split)),
            ("SE", (model_i_arr >= i_split) & (model_j_arr < j_split)),
        ]
        for k, (label, mask) in enumerate(subsets):
            ax = axes[k]
            ax.scatter(model_i_arr[mask] + 1, model_j_arr[mask] + 1, s=10, c="tab:orange", label="Selected model grids")
            ax.set_title(f"{label} quadrant (i split {i_split + 1:.1f}, j split {j_split + 1:.1f})")
            ax.grid(True, linewidth=0.3, alpha=0.4)

        for ax in axes[2:]:
            ax.set_xlabel("Model i")
        for ax in axes[0::2]:
            ax.set_ylabel("Model j")
        axes[0].legend(loc="best", frameon=False)
        fig.suptitle("Filtering centers vs selected model grids (4-panel)")
        fig.tight_layout()
        fig.savefig(args.plot_selected)
        plt.close(fig)
        print(f"Wrote {args.plot_selected}")


if __name__ == "__main__":
    main()
