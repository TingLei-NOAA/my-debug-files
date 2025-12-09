import argparse
import os
from typing import Tuple, Optional

import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np


def build_xy_from_dxdy(dx: np.ndarray, dy: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build physical x/y (meters) from dx/dy (meters) on a T-cell grid.
    Assumes dx, dy have shape (ny, nx) matching FV3 grid_yt/grid_xt.
    """
    ny, nx = dx.shape
    x = np.zeros((ny, nx), dtype=float)
    y = np.zeros((ny, nx), dtype=float)

    # Cumulative distance eastward and northward
    x[:, 1:] = np.cumsum(dx[:, :-1], axis=1)
    y[1:, :] = np.cumsum(dy[:-1, :], axis=0)
    return x, y


def load_xy(grid_spec: str) -> Tuple[np.ndarray, np.ndarray]:
    with nc.Dataset(grid_spec) as ds:
        dx = ds.variables["dx"][:]
        dy = ds.variables["dy"][:]
    return build_xy_from_dxdy(dx, dy)


def read_field(
    path: str,
    varname: str,
    level: int,
    y_slice: slice,
    x_slice: slice,
) -> Tuple[np.ndarray, Optional[str]]:
    """
    Read a 2D slice (y, x) from a NetCDF variable.
    Supports variables shaped (time, lev, y, x), (lev, y, x), or (y, x).
    """
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


def extract_window(
    X: int,
    Y: int,
    half: int,
    x_grid: np.ndarray,
    y_grid: np.ndarray,
) -> Tuple[slice, slice, np.ndarray, np.ndarray]:
    ny, nx = x_grid.shape
    x0 = max(0, X - half)
    x1 = min(nx, X + half + 1)
    y0 = max(0, Y - half)
    y1 = min(ny, Y + half + 1)
    xs = slice(x0, x1)
    ys = slice(y0, y1)
    return xs, ys, x_grid[ys, xs], y_grid[ys, xs]


def plot_two_panels(
    fig,
    axes,
    fields,
    meshes,
    titles,
    center_xy,
    units: Optional[str],
    cmap="turbo",
):
    def centers_to_edges(C: np.ndarray) -> np.ndarray:
        """
        Build an edge grid (ny+1, nx+1) from cell-center coordinates (ny, nx)
        by simple linear extrapolation on the borders and averaging corners.
        Works for curvilinear grids to satisfy pcolormesh edge input.
        """
        ny, nx = C.shape
        # Pad in x
        left = C[:, [0]] - (C[:, [1]] - C[:, [0]])
        right = C[:, [-1]] + (C[:, [-1]] - C[:, [-2]])
        Cx = np.concatenate([left, C, right], axis=1)
        # Pad in y
        top = Cx[[0], :] - (Cx[[1], :] - Cx[[0], :])
        bottom = Cx[[-1], :] + (Cx[[-1], :] - Cx[[-2], :])
        Cy = np.concatenate([top, Cx, bottom], axis=0)
        # Average four surrounding points to get edges
        edges = 0.25 * (
            Cy[:-1, :-1] + Cy[:-1, 1:] + Cy[1:, :-1] + Cy[1:, 1:]
        )
        return edges

    for ax, field, (Xmesh, Ymesh), title in zip(axes, fields, meshes, titles):
        xe = centers_to_edges(Xmesh) / 1000.0
        ye = centers_to_edges(Ymesh) / 1000.0
        pcm = ax.pcolormesh(xe, ye, field, shading="auto", cmap=cmap)
        ax.plot(center_xy[0] / 1000.0, center_xy[1] / 1000.0, "ro", ms=4, label="(X,Y)")
        ax.set_aspect("equal")
        ax.set_xlabel("X (km)")
        ax.set_ylabel("Y (km)")
        ax.set_title(title)
        ax.legend(loc="best", fontsize=8)
        cbar = fig.colorbar(pcm, ax=ax)
        cbar.set_label(units or "")


def main():
    parser = argparse.ArgumentParser(description="Horizontal contours around (X,Y) using physical FV3 x/y.")
    parser.add_argument(
        "--file1",
        default="Data/mgbf_NA/20240527.010000.p484-inho-ref_L30_proc5_2G_35kmv6-mgbf-D1-p64_dirac_SABER_lam.fv_core.res.nc",
        help="First NetCDF file (default: run1.nc).",
    )
    parser.add_argument(
        "--file2",
        default="Data/mgbf_NA/20240527.010000.p484-inho-test_L30_proc5_2G_35kmv6-mgbf-D1-p64_dirac_SABER_lam.fv_core.res.nc",
        help="Second NetCDF file (default: run2.nc).",
    )
    parser.add_argument("--label1", default=None, help="Title label for first plot.")
    parser.add_argument("--label2", default=None, help="Title label for second plot.")
    parser.add_argument("--grid-spec", default="fv3_grid_dxdy.nc", help="NetCDF with dx/dy on T grid.")
    parser.add_argument("--variable", default="air_temperature", help="Variable name to plot.")
    parser.add_argument("--level", type=int, default=0, help="Vertical level index (0-based).")
    parser.add_argument("--x-index", type=int, default=None, help="Center X grid index (0-based, grid_xt).")
    parser.add_argument("--y-index", type=int, default=None, help="Center Y grid index (0-based, grid_yt).")
    parser.add_argument(
        "--half-width",
        "--lgh",
        dest="half_width",
        type=int,
        default=10,
        help="Half-width (grid points) of the square window around (X,Y).",
    )
    parser.add_argument("--output", default="dev.png", help="Output figure filename.")
    args = parser.parse_args()

    # Sanity-check files exist
    for path in (args.file1, args.file2):
        if not os.path.exists(path):
            raise FileNotFoundError(f"Input file not found: {path}")

    x_grid, y_grid = load_xy(args.grid_spec)
    ny, nx = x_grid.shape

    if args.x_index is None or args.y_index is None:
        # Auto-center on the global max of file1 at the requested level
        with nc.Dataset(args.file1) as ds:
            var = ds.variables[args.variable]
            if var.ndim == 4:
                base = var[0, args.level, :, :]
            elif var.ndim == 3:
                base = var[args.level, :, :]
            elif var.ndim == 2:
                base = var[:, :]
            else:
                raise ValueError(f"Unsupported rank {var.ndim} for {args.variable} in {args.file1}")
        max_idx = np.nanargmax(base)
        Yc, Xc = np.unravel_index(max_idx, base.shape)
        print(
            f"Auto-centered on global max of {args.variable} (level {args.level}) "
            f"at idx (X,Y)=({Xc},{Yc}), value={base[Yc, Xc]:.3g}"
        )
    else:
        Xc = args.x_index
        Yc = args.y_index

    if not (0 <= Xc < nx and 0 <= Yc < ny):
        raise ValueError(f"(X,Y)=({Xc},{Yc}) outside grid ({nx}, {ny}).")

    xs, ys, Xmesh, Ymesh = extract_window(Xc, Yc, args.half_width, x_grid, y_grid)

    field1, units1 = read_field(args.file1, args.variable, args.level, ys, xs)
    field2, units2 = read_field(args.file2, args.variable, args.level, ys, xs)
    units = units1 or units2

    center_xy = (x_grid[Yc, Xc], y_grid[Yc, Xc])
    print(f"Center indices (X,Y) = ({Xc},{Yc}); physical (km) = ({center_xy[0]/1000:.3f}, {center_xy[1]/1000:.3f})")
    titles = [
        args.label1 or os.path.basename(args.file1),
        args.label2 or os.path.basename(args.file2),
    ]

    # Debug info about window and data ranges
    win_max1 = np.nanmax(field1)
    win_max2 = np.nanmax(field2)
    loc1 = np.unravel_index(np.nanargmax(field1), field1.shape)
    loc2 = np.unravel_index(np.nanargmax(field2), field2.shape)
    abs_loc1 = (xs.start + loc1[1], ys.start + loc1[0])
    abs_loc2 = (xs.start + loc2[1], ys.start + loc2[0])
    print(
        f"Window x: [{xs.start}:{xs.stop}), y: [{ys.start}:{ys.stop}), shape={field1.shape}\n"
        f"{args.variable} min/max file1=({np.nanmin(field1):.3g},{np.nanmax(field1):.3g}) "
        f"file2=({np.nanmin(field2):.3g},{np.nanmax(field2):.3g})\n"
        f"Window max file1={win_max1:.3g} at idx (X,Y)=({abs_loc1[0]},{abs_loc1[1]}); "
        f"file2={win_max2:.3g} at idx (X,Y)=({abs_loc2[0]},{abs_loc2[1]})"
    )

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    plot_two_panels(
        fig,
        axes,
        fields=[field1, field2],
        meshes=[(Xmesh, Ymesh), (Xmesh, Ymesh)],
        titles=titles,
        center_xy=center_xy,
        units=units,
    )
    fig.suptitle(f"{args.variable} level={args.level}, window +/-{args.half_width} pts around (X={Xc}, Y={Yc})")
    fig.savefig(args.output, dpi=150)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
