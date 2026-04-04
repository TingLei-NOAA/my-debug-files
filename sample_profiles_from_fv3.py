"""
Sample vertical profiles from FV3 NetCDF at selected i/j points.

Inputs:
  - filter_to_model_map_indices.txt (columns: index center_filt_i center_filt_j model_i model_j)
  - FV3 NetCDF file
  - variable list (file or comma-separated)
  - k index list (file or comma-separated)

Output:
  - One text file per index, containing (k, value) rows.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import netCDF4 as nc
import yaml

TRACER_VAR_LIST = [
    "water_vapor_mixing_ratio_wrt_moist_air",
    "cloud_liquid_ice",
    "cloud_liquid_water",
    "rain_water",
    "snow_water",
    "graupel",
    "ozone_mass_mixing_ratio",
    "rain_number_concentration",
    "cloud_ice_number_concentration",
    "cloud_droplet_number_concentration",
    "xvar_1",
    "xvar_2",
    "xvar_3",
    "xvar_4",
    "xvar_5",
]

DYN_VAR_LIST = [
    "upward_air_velocity",
    "eastward_wind",
    "northward_wind",
    "air_temperature",
    "air_pressure_thickness",
]


def read_list(path_or_csv: str, cast=int) -> list:
    p = Path(path_or_csv)
    if p.exists():
        items = []
        for line in p.read_text().splitlines():
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            for token in line.replace(",", " ").split():
                items.append(cast(token))
        return items
    items = [x for x in path_or_csv.replace(",", " ").split() if x]
    return [cast(x) for x in items]


def read_vars(path_or_csv: str) -> list[str]:
    p = Path(path_or_csv)
    if p.exists():
        vars_out = []
        for line in p.read_text().splitlines():
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            for token in line.replace(",", " ").split():
                vars_out.append(token)
        return vars_out
    return [x for x in path_or_csv.replace(",", " ").split() if x]


def read_dirac_yaml(path: Path) -> tuple[list[str], list[int]]:
    data = yaml.safe_load(path.read_text())
    if not isinstance(data, dict) or "dirac" not in data:
        raise SystemExit(f"Missing 'dirac' section in {path}")
    dirac = data["dirac"]
    if "ifdir" not in dirac or "ildir" not in dirac:
        raise SystemExit(f"Missing ifdir/ildir in {path}")
    vars_list = list(dirac["ifdir"])
    k_list = list(dirac["ildir"])
    if not vars_list or not k_list:
        raise SystemExit(f"Empty ifdir/ildir in {path}")
    if len(vars_list) != len(k_list):
        raise SystemExit(f"ifdir/ildir length mismatch in {path}")
    return vars_list, [int(x) for x in k_list]


def ensure_list(name: str, value):
    if isinstance(value, (int, float, str)):
        return [value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, list):
        return value
    if isinstance(value, tuple):
        return list(value)
    raise SystemExit(f"{name} must be a list; got {type(value).__name__}")


def parse_map_file(path: Path) -> list[tuple[int, int, int]]:
    rows = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 5:
            continue
        idx = int(parts[0])
        model_i = int(parts[3])
        model_j = int(parts[4])
        rows.append((idx, model_i, model_j))
    return rows


def parse_ij_pairs(path_or_csv: str) -> list[tuple[int, int]]:
    p = Path(path_or_csv)
    pairs: list[tuple[int, int]] = []
    tokens: list[str] = []
    if p.exists():
        for line in p.read_text().splitlines():
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            tokens.extend(line.replace(",", " ").split())
    else:
        tokens = [x for x in path_or_csv.replace(",", " ").replace(";", " ").split() if x]
    if len(tokens) % 2 != 0:
        raise SystemExit("Group list must have an even number of entries (i j pairs).")
    it = iter(tokens)
    for i_str, j_str in zip(it, it):
        pairs.append((int(i_str), int(j_str)))
    return pairs


def infer_index_order(var, i_idx: int, j_idx: int, k_idx: int, time_idx: int | None):
    dims = tuple(var.dimensions)
    shape = tuple(var.shape)
    slicer = []

    for dim_name, dim_size in zip(dims, shape):
        dim = dim_name.lower()
        if dim in {"time", "t", "nt"}:
            idx = time_idx if time_idx is not None else 0
        elif dim.startswith("zaxis") or dim in {"lev", "pfull", "phalf", "nz"}:
            idx = k_idx
        elif dim.startswith("yaxis") or dim in {"grid_yt", "grid_y"}:
            idx = j_idx
        elif dim.startswith("xaxis") or dim in {"grid_xt", "grid_x"}:
            idx = i_idx
        elif dim in {"lon", "longitude", "lat", "latitude"}:
            raise SystemExit(f"Unsupported 1D lat/lon dimension on variable {var.name()} with dims {dims}")
        else:
            raise SystemExit(f"Cannot infer dimension role for variable {var.name()} dim '{dim_name}' with dims {dims}")

        if idx < 0 or idx >= dim_size:
            raise SystemExit(
                f"Index out of bounds for variable {var.name()}: dim {dim_name} size={dim_size}, requested index={idx}"
            )
        slicer.append(idx)

    return tuple(slicer)


def infer_dim_roles(var) -> dict[str, int]:
    roles: dict[str, int] = {}
    dims = tuple(var.dimensions)
    for pos, dim_name in enumerate(dims):
        dim = dim_name.lower()
        if dim in {"time", "t", "nt"}:
            roles["t"] = pos
        elif dim.startswith("zaxis") or dim in {"lev", "pfull", "phalf", "nz"}:
            roles["k"] = pos
        elif dim.startswith("yaxis") or dim in {"grid_yt", "grid_y"}:
            roles["j"] = pos
        elif dim.startswith("xaxis") or dim in {"grid_xt", "grid_x"}:
            roles["i"] = pos
    return roles


def plot_profiles_for_group(
    ds_dyn: nc.Dataset,
    ds_tracer: nc.Dataset,
    points: list[tuple[int, int]],
    vars_list: list[str],
    k_list_in: list[int],
    k_list0: list[int],
    time_index: int,
    out_dir: Path,
) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError as e:  # pragma: no cover
        raise SystemExit("matplotlib is required for plotting profiles") from e

    out_dir.mkdir(parents=True, exist_ok=True)
    for (i_s, j_s) in points:
        fig, ax = plt.subplots(figsize=(6, 6), dpi=150)
        values = []
        for var_name, k_s in zip(vars_list, k_list0):
            if var_name in TRACER_VAR_LIST:
                ds = ds_tracer
            elif var_name in DYN_VAR_LIST:
                ds = ds_dyn
            else:
                raise SystemExit(f"Variable not in tracer/dyn lists: {var_name}")
            if var_name not in ds.variables:
                raise SystemExit(f"Variable not found in file: {var_name}")
            var = ds.variables[var_name]
            slicer = infer_index_order(var, i_s, j_s, k_s, time_index)
            val = var[slicer]
            try:
                val = float(val)
            except TypeError:
                val = float(np.asarray(val).squeeze())
            values.append(val)
        ax.plot(values, k_list_in, marker="o", linewidth=1.0)
        for x, y, name in zip(values, k_list_in, vars_list):
            ax.text(x, y, name, fontsize=7, ha="left", va="center")
        ax.set_xlabel("Value")
        ax.set_ylabel("k")
        ax.set_title(f"Profile at i={i_s + 1}, j={j_s + 1}")
        ax.grid(True, linewidth=0.3, alpha=0.4)
        fig.tight_layout()
        out_path = out_dir / f"profile_i{i_s + 1:04d}_j{j_s + 1:04d}.png"
        fig.savefig(out_path)
        plt.close(fig)
        print(f"Wrote {out_path}")


def select_dataset(var_name: str, ds_dyn: nc.Dataset, ds_tracer: nc.Dataset) -> nc.Dataset:
    if var_name in TRACER_VAR_LIST:
        return ds_tracer
    if var_name in DYN_VAR_LIST:
        return ds_dyn
    raise SystemExit(f"Variable not in tracer/dyn lists: {var_name}")


def sample_scalar(var, i_s: int, j_s: int, k_s: int, time_index: int) -> float:
    slicer = infer_index_order(var, i_s, j_s, k_s, time_index)
    val = var[slicer]
    try:
        return float(val)
    except TypeError:
        return float(np.asarray(val).squeeze())


def extract_level_field(var, k_s: int, time_index: int) -> np.ndarray:
    roles = infer_dim_roles(var)
    dims = tuple(var.dimensions)
    shape = tuple(var.shape)
    slicer = []
    for pos, (dim_name, dim_size) in enumerate(zip(dims, shape)):
        if roles.get("t") == pos:
            idx = time_index
            if idx < 0 or idx >= dim_size:
                raise SystemExit(f"time index out of bounds for {var.name}: {idx} not in [0,{dim_size})")
            slicer.append(idx)
        elif roles.get("k") == pos:
            if k_s < 0 or k_s >= dim_size:
                raise SystemExit(f"k index out of bounds for {var.name}: {k_s} not in [0,{dim_size})")
            slicer.append(k_s)
        elif roles.get("j") == pos or roles.get("i") == pos:
            slicer.append(slice(None))
        else:
            raise SystemExit(f"Cannot build horizontal field slice for {var.name} with dims {dims}")
    field = np.asarray(var[tuple(slicer)])
    if field.ndim != 2:
        raise SystemExit(f"Expected 2D horizontal field for {var.name}, got shape {field.shape}")
    return field


def plot_horizontal_diagnostic(
    ds_dyn: nc.Dataset,
    ds_tracer: nc.Dataset,
    points0: list[tuple[int, int, int]],
    var_name: str,
    k_in: int,
    k0: int,
    time_index: int,
    out_path: Path,
) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError as e:  # pragma: no cover
        raise SystemExit("matplotlib is required for diagnostic plotting") from e

    ds = select_dataset(var_name, ds_dyn, ds_tracer)
    if var_name not in ds.variables:
        raise SystemExit(f"Variable not found in file: {var_name}")
    var = ds.variables[var_name]

    xs = []
    ys = []
    vals = []
    for _, i_s, j_s in points0:
        xs.append(i_s + 1)
        ys.append(j_s + 1)
        vals.append(sample_scalar(var, i_s, j_s, k0, time_index))

    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    vals = np.asarray(vals, dtype=float)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 6), dpi=150)
    if xs.size >= 3:
        cf = ax.tricontourf(xs, ys, vals, levels=20, cmap="magma")
        fig.colorbar(cf, ax=ax, pad=0.02, label=f"{var_name} at k={k_in}")
    sc = ax.scatter(xs, ys, c=vals, s=10, cmap="magma", edgecolors="k", linewidths=0.1)
    if xs.size < 3:
        fig.colorbar(sc, ax=ax, pad=0.02, label=f"{var_name} at k={k_in}")
    ax.set_xlabel("Model i")
    ax.set_ylabel("Model j")
    ax.set_title(f"Sampled horizontal values: {var_name} at k={k_in}")
    ax.grid(True, linewidth=0.3, alpha=0.4)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    print(f"Wrote {out_path}")


def plot_full_field_with_samples(
    ds_dyn: nc.Dataset,
    ds_tracer: nc.Dataset,
    points0: list[tuple[int, int, int]],
    var_name: str,
    k_in: int,
    k0: int,
    time_index: int,
    out_path: Path,
) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError as e:  # pragma: no cover
        raise SystemExit("matplotlib is required for diagnostic plotting") from e

    ds = select_dataset(var_name, ds_dyn, ds_tracer)
    if var_name not in ds.variables:
        raise SystemExit(f"Variable not found in file: {var_name}")
    var = ds.variables[var_name]
    field = extract_level_field(var, k0, time_index)

    xs = []
    ys = []
    vals = []
    labels = []
    for _, i_s, j_s in points0:
        labels.append(_)
        xs.append(i_s)
        ys.append(j_s)
        vals.append(sample_scalar(var, i_s, j_s, k0, time_index))
    xs = np.asarray(xs, dtype=int)
    ys = np.asarray(ys, dtype=int)
    vals = np.asarray(vals, dtype=float)

    full_max_idx = np.unravel_index(int(np.nanargmax(field)), field.shape)
    full_max_val = float(np.nanmax(field))
    sampled_max_val = float(np.nanmax(vals))
    sampled_argmax = int(np.nanargmax(vals))
    sampled_max_i = int(xs[sampled_argmax]) + 1
    sampled_max_j = int(ys[sampled_argmax]) + 1
    print(
        f"[diag] {var_name} k={k_in}: full max={full_max_val:.6g} at (i={full_max_idx[1] + 1}, j={full_max_idx[0] + 1}); "
        f"sampled max={sampled_max_val:.6g} at (i={sampled_max_i}, j={sampled_max_j})"
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(10, 6), dpi=150)
    im = ax.imshow(field, origin="lower", cmap="magma", aspect="auto")
    fig.colorbar(im, ax=ax, pad=0.02, label=f"{var_name} at k={k_in}")
    sc = ax.scatter(xs, ys, c=vals, cmap="magma", s=16, edgecolors="cyan", linewidths=0.3)
    ax.scatter([full_max_idx[1]], [full_max_idx[0]], marker="x", c="lime", s=60, linewidths=1.0)
    for label, x, y in zip(labels, xs, ys):
        ax.text(x + 6, y + 6, str(label), color="white", fontsize=5, alpha=0.8)
    ax.set_xlabel("Model i (0-based)")
    ax.set_ylabel("Model j (0-based)")
    ax.set_title(f"Full field with sampled points: {var_name} at k={k_in}")
    ax.grid(True, linewidth=0.2, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    print(f"Wrote {out_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Sample vertical profiles from FV3 NetCDF at selected points.")
    parser.add_argument("--map", type=Path, default=Path("filter_to_model_map_indices.txt"), help="Map file with model i/j")
    parser.add_argument("--fv3-dyn", type=Path, default="Data/mgbf_NA/20240527.010000.p484-L30_proc5_2G_35kmv6-mgbf-D1-test-18var_dirac_SABER_lam.fv_core.res.nc",
                          help="FV3 dynamics NetCDF file")
    parser.add_argument("--fv3-tracer", type=Path, default="Data/mgbf_NA/20240527.010000.p484-L30_proc5_2G_35kmv6-mgbf-D1-test-18var_dirac_SABER_lam.fv_tracer.res.nc", 
                          help="FV3 tracer NetCDF file")
    parser.add_argument("--vars", help="Variable list (file or comma-separated)")
    parser.add_argument("--k-list", help="k index list (file or comma-separated)")
    parser.add_argument("--dirac-yaml", type=Path,default="test-18var.yaml", help="YAML file with dirac.ifdir and dirac.ildir")
    parser.add_argument("--k-base", type=int, default=1, choices=[0, 1], help="k index base in input list")
    parser.add_argument("--ij-base", type=int, default=1, choices=[0, 1], help="i/j base in map file")
    parser.add_argument("--time-index", type=int, default=0, help="Time index for variables with time dimension")
    parser.add_argument("--out-dir", type=Path, default=Path("profiles_out"), help="Output directory")
    parser.add_argument("--plot-group", help="Group of i/j points (file or comma/space list)")
    parser.add_argument("--plot_sub_domain_index",default="12,30,40", help="Subdomain indices to plot (file or comma/space list)")
    parser.add_argument("--plot-group-ij-base", type=int, default=1, choices=[0, 1], help="i/j base in group list")
    parser.add_argument("--plot-out-dir", type=Path, default=Path("profiles_plots"), help="Output directory for plots")
    parser.add_argument(
        "--plot-first-level-map",
        type=Path,
        help="Write a contour/scatter diagnostic plot for the first (var, k) pair over all sampled points",
    )
    parser.add_argument(
        "--plot-first-level-full-field",
        type=Path,
        help="Write a full horizontal field plot for the first (var, k) pair with sampled points overlaid",
    )
    args = parser.parse_args()

    if args.dirac_yaml:
        vars_list, k_list = read_dirac_yaml(args.dirac_yaml)
    else:
        if not args.vars or not args.k_list:
            raise SystemExit("Provide --vars and --k-list, or use --dirac-yaml.")
        k_list = read_list(args.k_list, cast=int)
        vars_list = read_vars(args.vars)
    vars_list = ensure_list("vars_list", vars_list)
    k_list = ensure_list("k_list", k_list)
    if not k_list:
        raise SystemExit("Empty k list")
    if not vars_list:
        raise SystemExit("Empty variable list")
    k_list = [int(x) for x in k_list]
    vars_list = [str(x) for x in vars_list]

    points = parse_map_file(args.map)
    if not points:
        raise SystemExit(f"No points found in {args.map}")

    # Convert to 0-based
    k0 = [k - args.k_base for k in k_list]
    points0 = [(idx, i - args.ij_base, j - args.ij_base) for idx, i, j in points]

    args.out_dir.mkdir(parents=True, exist_ok=True)

    with nc.Dataset(args.fv3_dyn) as ds_dyn, nc.Dataset(args.fv3_tracer) as ds_tracer:
        if args.plot_first_level_map:
            plot_horizontal_diagnostic(
                ds_dyn=ds_dyn,
                ds_tracer=ds_tracer,
                points0=points0,
                var_name=vars_list[0],
                k_in=k_list[0],
                k0=k0[0],
                time_index=args.time_index,
                out_path=args.plot_first_level_map,
            )
        if args.plot_first_level_full_field:
            plot_full_field_with_samples(
                ds_dyn=ds_dyn,
                ds_tracer=ds_tracer,
                points0=points0,
                var_name=vars_list[0],
                k_in=k_list[0],
                k0=k0[0],
                time_index=args.time_index,
                out_path=args.plot_first_level_full_field,
            )
        if args.plot_sub_domain_index:
            idx_list = read_list(args.plot_sub_domain_index, cast=int)
            idx_map = {idx: (i, j) for idx, i, j in points0}
            group_pts0 = []
            for idx in idx_list:
                if idx not in idx_map:
                    raise SystemExit(f"Subdomain index not found in map: {idx}")
                group_pts0.append(idx_map[idx])
            plot_profiles_for_group(
                ds_dyn=ds_dyn,
                ds_tracer=ds_tracer,
                points=group_pts0,
                vars_list=vars_list,
                k_list_in=k_list,
                k_list0=k0,
                time_index=args.time_index,
                out_dir=args.plot_out_dir,
            )
        elif args.plot_group:
            group_pts = parse_ij_pairs(args.plot_group)
            group_pts0 = [(i - args.plot_group_ij_base, j - args.plot_group_ij_base) for i, j in group_pts]
            plot_profiles_for_group(
                ds_dyn=ds_dyn,
                ds_tracer=ds_tracer,
                points=group_pts0,
                vars_list=vars_list,
                k_list_in=k_list,
                k_list0=k0,
                time_index=args.time_index,
                out_dir=args.plot_out_dir,
            )
        for idx, i_s, j_s in points0:
            out_path = args.out_dir / f"profile_subdomain_{idx:04d}.txt"
            with out_path.open("w") as f:
                f.write("# index  number of sampling vertical levels i j\n")
                f.write(f"{idx} {len(k_list)} {i_s + 1} {j_s + 1}\n")
                f.write("# k value\n")
                sampled = []
                for var_name, k_in, k_s in zip(vars_list, k_list, k0):
                    ds = select_dataset(var_name, ds_dyn, ds_tracer)
                    if var_name not in ds.variables:
                        raise SystemExit(f"Variable not found in file: {var_name}")
                    var = ds.variables[var_name]
                    val = sample_scalar(var, i_s, j_s, k_s, args.time_index)
                    sampled.append((k_in, val))
                    f.write(f"{k_in} {val:.6g}\n")
                f.write("# now the complete profile with fillled values\n")
                k_min = k_list[0]
                k_max = k_list[-1]
                k_vals = [k for k, _ in sampled]
                v_vals = [v for _, v in sampled]
                lo = 0
                for k_full in range(k_min, k_max + 1):
                    if k_full == k_vals[lo]:
                        v_full = v_vals[lo]
                    else:
                        while lo < len(k_vals) - 1 and k_vals[lo + 1] < k_full:
                            lo += 1
                        hi = min(lo + 1, len(k_vals) - 1)
                        if k_vals[hi] == k_vals[lo]:
                            v_full = v_vals[lo]
                        else:
                            v0 = v_vals[lo]
                            v1 = v_vals[hi]
                            k0_val = k_vals[lo]
                            k1_val = k_vals[hi]
                            v_full = v0 + (v1 - v0) * (k_full - k0_val) / (k1_val - k0_val)
                    f.write(f"{k_full} {v_full:.6g}\n")
            print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
