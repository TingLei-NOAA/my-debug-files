import argparse
import os
from typing import List, Optional, Tuple

import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np


def build_xy_from_dxdy(dx: np.ndarray, dy: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    ny, nx = dx.shape
    x = np.zeros((ny, nx), dtype=float)
    y = np.zeros((ny, nx), dtype=float)
    x[:, 1:] = np.cumsum(dx[:, :-1], axis=1)
    y[1:, :] = np.cumsum(dy[:-1, :], axis=0)
    return x, y


def load_xy(grid_spec: str) -> Tuple[np.ndarray, np.ndarray]:
    with nc.Dataset(grid_spec) as ds:
        dx = ds.variables["dx"][:]
        dy = ds.variables["dy"][:]
    return build_xy_from_dxdy(dx, dy)


def read_field(path: str, varname: str, level: int, y_slice: slice, x_slice: slice) -> Tuple[np.ndarray, Optional[str]]:
    with nc.Dataset(path) as ds:
        var = ds.variables[varname]
        units = getattr(var, "units", None)
        if var.ndim == 4:
            data = var[0, level, y_slice, x_slice]
        elif var.ndim == 3:
            data = var[level, y_slice, x_slice]
        elif var.ndim == 2:
            if level not in (0, None):
                raise ValueError(f"Variable {varname} is 2D; use level=0.")
            data = var[y_slice, x_slice]
        else:
            raise ValueError(f"Unsupported rank {var.ndim} for {varname} in {path}")
    return np.array(data), units


def extract_window(X: int, Y: int, half: int, x_grid: np.ndarray, y_grid: np.ndarray) -> Tuple[slice, slice, np.ndarray, np.ndarray]:
    ny, nx = x_grid.shape
    x0 = max(0, X - half)
    x1 = min(nx, X + half + 1)
    y0 = max(0, Y - half)
    y1 = min(ny, Y + half + 1)
    xs = slice(x0, x1)
    ys = slice(y0, y1)
    return xs, ys, x_grid[ys, xs], y_grid[ys, xs]


def centers_to_edges(C: np.ndarray) -> np.ndarray:
    ny, nx = C.shape
    left = C[:, [0]] - (C[:, [1]] - C[:, [0]])
    right = C[:, [-1]] + (C[:, [-1]] - C[:, [-2]])
    Cx = np.concatenate([left, C, right], axis=1)
    top = Cx[[0], :] - (Cx[[1], :] - Cx[[0], :])
    bottom = Cx[[-1], :] + (Cx[[-1], :] - Cx[[-2], :])
    Cy = np.concatenate([top, Cx, bottom], axis=0)
    edges = 0.25 * (Cy[:-1, :-1] + Cy[:-1, 1:] + Cy[1:, :-1] + Cy[1:, 1:])
    return edges


def main():
    parser = argparse.ArgumentParser(description="Horizontal contours around multiple centers for two files.")
    parser.add_argument(
        "--file1",
        default="Data/mgbf_NA/20240527.010000.p484-inho-ref_L30_proc5_2G_35kmv6-mgbf-D1-p64_dirac_SABER_lam.fv_core.res.nc",
        help="First NetCDF file (default: run1.nc).",
    )
    parser.add_argument(
        "--file2",
        default="Data/mgbf_NA/20240527.010000.p484-inho_L30_proc5_2G_35kmv6-mgbf-D1-p64_dirac_SABER_lam.fv_core.res.nc",
        help="Second NetCDF file (default: run2.nc).",
    )
    parser.add_argument("--label1", default=None, help="Title label for first file.")
    parser.add_argument("--label2", default=None, help="Title label for second file.")
    parser.add_argument("--grid-spec", default="fv3_grid_dxdy.nc", help="NetCDF with dx/dy on T grid.")
    parser.add_argument("--variable", default="air_temperature", help="Variable name to plot.")
    parser.add_argument("--level", type=int, default=29, help="Vertical level index (0-based).")
    parser.add_argument("--x-indices", nargs="+", type=int, default=[149,1899], help="List of X indices (0-based).")
    parser.add_argument("--y-indices", nargs="+", type=int, default=[149,1799], help="List of Y indices (0-based).")
    parser.add_argument("--half-width", "--lgh", dest="half_width", type=int, default=50, help="Half-width (grid points) of the square window.")
    parser.add_argument("--titles", type=str, nargs="*", default=None, help="Titles for subplots; length should match centers*files. Default: 'xxxx'.")
    parser.add_argument("--output", default="dev_multi.png", help="Output figure filename.")
    args = parser.parse_args()

    # Check files
    for path in (args.file1, args.file2):
        if not os.path.exists(path):
            raise FileNotFoundError(f"Input file not found: {path}")

    x_grid, y_grid = load_xy(args.grid_spec)
    ny, nx = x_grid.shape

    if args.x_indices is None or args.y_indices is None:
        x_list = [nx // 2]
        y_list = [ny // 2]
    else:
        x_list = args.x_indices
        y_list = args.y_indices
    if len(x_list) != len(y_list):
        raise ValueError("x-indices and y-indices must have the same length.")
    centers = list(zip(x_list, y_list))

    files = [args.file1, args.file2]
    labels = [args.label1 or os.path.basename(args.file1), args.label2 or os.path.basename(args.file2)]

    nrows = len(centers)
    ncols = len(files)

    titles = args.titles
    if not titles or len(titles) != nrows * ncols:
        titles = ["xxxx"] * (nrows * ncols)

    # Precompute all windows and normalized fields
    entries = []
    global_min = np.inf
    global_max = -np.inf
    idx_title = 0
    curve_labels = ["label1", "label2", "label3", "label4"]
    for r, (Xc, Yc) in enumerate(centers):
        if not (0 <= Xc < nx and 0 <= Yc < ny):
            raise ValueError(f"Center (X,Y)=({Xc},{Yc}) outside grid ({nx},{ny}).")
        xs, ys, Xmesh, Ymesh = extract_window(Xc, Yc, args.half_width, x_grid, y_grid)
        center_xy = (x_grid[Yc, Xc], y_grid[Yc, Xc])
        for c, (fpath, flabel) in enumerate(zip(files, labels)):
            field, units = read_field(fpath, args.variable, args.level, ys, xs)
            max_val = np.nanmax(field)
            if not np.isfinite(max_val) or max_val == 0:
                field_norm = np.zeros_like(field)
            else:
                field_norm = field / max_val
            global_min = min(global_min, np.nanmin(field_norm))
            global_max = max(global_max, np.nanmax(field_norm))
            # profile along x through the maximum point in this window
            flat_idx = np.nanargmax(field_norm)
            jj, ii = np.unravel_index(flat_idx, field_norm.shape)
            x_profile = Xmesh[jj, :] / 1000.0
            y_profile = field_norm[jj, :]
            entries.append(
                {
                    "r": r,
                    "c": c,
                    "field": field_norm,
                    "units": units,
                    "xe": centers_to_edges(Xmesh) / 1000.0,
                    "ye": centers_to_edges(Ymesh) / 1000.0,
                    "center_xy": center_xy,
                    "title": titles[idx_title],
                    "src": fpath,
                    "orig_min": float(np.nanmin(field)),
                    "orig_max": float(np.nanmax(field)),
                    "prof_x": x_profile,
                    "prof_y": y_profile,
                }
            )
            idx_title += 1

    if not np.isfinite(global_min):
        global_min = 0.0
    if not np.isfinite(global_max) or global_max == 0:
        global_max = 1.0

    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 4.5 * nrows), constrained_layout=True, squeeze=False)

    for entry in entries:
        ax = axes[entry["r"], entry["c"]]
        pcm = ax.pcolormesh(entry["xe"], entry["ye"], entry["field"], shading="auto", cmap="turbo", vmin=global_min, vmax=global_max)
        ax.plot(entry["center_xy"][0] / 1000.0, entry["center_xy"][1] / 1000.0, "ro", ms=4, label="(X,Y)")
        ax.set_aspect("equal")
        ax.set_xlabel("X (km)")
        ax.set_ylabel("Y (km)")
        ax.set_title(entry["title"])
        ax.legend(loc="best", fontsize=8)
        cbar = fig.colorbar(pcm, ax=ax)
        cbar.set_label(entry["units"] or "")
        print(
            f"[{entry['src']}] center (X,Y)=({centers[entry['r']][0]},{centers[entry['r']][1]}) "
            f"window orig min/max=({entry['orig_min']:.3g},{entry['orig_max']:.3g}) "
            f"norm min/max=({np.nanmin(entry['field']):.3g},{np.nanmax(entry['field']):.3g})"
        )

    fig.suptitle(f"{args.variable} level={args.level}, window +/-{args.half_width} pts; centers={centers}")
    fig.savefig(args.output, dpi=150)
    print(f"Wrote {args.output}")

    # Second figure: profiles (x-direction through window maximum) per file
    fig2, axes2 = plt.subplots(1, ncols, figsize=(6 * ncols, 4), constrained_layout=True, squeeze=False)
    for c in range(ncols):
        ax = axes2[0, c]
        for entry in entries:
            if entry["c"] != c:
                continue
            lbl_idx = entry["r"]
            lbl = curve_labels[lbl_idx] if lbl_idx < len(curve_labels) else f"center{lbl_idx}"
            ax.plot(entry["prof_x"], entry["prof_y"], label=lbl)
        ax.set_xlabel("X (km)")
        ax.set_ylabel(f"{args.variable} (normalized)")
        ax.set_title(titles[c] if c < len(titles) else labels[c])
        ax.grid(True, linestyle=":")
        ax.legend()
    fig2.suptitle(f"{args.variable} level={args.level} profiles (x through max); centers={centers}")
    profile_out = os.path.splitext(args.output)[0] + "_profiles.png"
    fig2.savefig(profile_out, dpi=150)
    print(f"Wrote {profile_out}")


if __name__ == "__main__":
    main()
