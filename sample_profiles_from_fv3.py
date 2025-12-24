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


def guess_dim_role(dim_name: str) -> str:
    name = dim_name.lower()
    if any(k in name for k in ["lev", "pfull", "phalf", "nz", "z"]):
        return "k"
    if any(k in name for k in ["lat", "yt", "y", "nj"]):
        return "j"
    if any(k in name for k in ["lon", "xt", "x", "ni"]):
        return "i"
    return ""


def infer_index_order(var, i_idx: int, j_idx: int, k_idx: int, time_idx: int | None):
    dims = var.dimensions
    roles = [guess_dim_role(d) for d in dims]
    order = []
    for d, r in zip(dims, roles):
        if r == "k":
            order.append(("k", d))
        elif r == "j":
            order.append(("j", d))
        elif r == "i":
            order.append(("i", d))
        elif d.lower() in ["time", "t", "nt"]:
            order.append(("t", d))
        else:
            order.append(("", d))

    # Build slicer with defaults, assume typical (time, k, j, i) or (k, j, i)
    slicer = []
    for role, d in order:
        if role == "t":
            slicer.append(time_idx if time_idx is not None else 0)
        elif role == "k":
            slicer.append(k_idx)
        elif role == "j":
            slicer.append(j_idx)
        elif role == "i":
            slicer.append(i_idx)
        else:
            slicer.append(0)
    return tuple(slicer)


def plot_profiles_for_group(
    ds: nc.Dataset,
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
            if var_name not in ds.variables:
                raise SystemExit(f"Variable not found: {var_name}")
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


def main() -> None:
    parser = argparse.ArgumentParser(description="Sample vertical profiles from FV3 NetCDF at selected points.")
    parser.add_argument("--map", type=Path, default=Path("filter_to_model_map_indices.txt"), help="Map file with model i/j")
    parser.add_argument("--fv3", type=Path, required=True, help="FV3 NetCDF file to sample")
    parser.add_argument("--vars", required=True, help="Variable list (file or comma-separated)")
    parser.add_argument("--k-list", required=True, help="k index list (file or comma-separated)")
    parser.add_argument("--k-base", type=int, default=1, choices=[0, 1], help="k index base in input list")
    parser.add_argument("--ij-base", type=int, default=1, choices=[0, 1], help="i/j base in map file")
    parser.add_argument("--time-index", type=int, default=0, help="Time index for variables with time dimension")
    parser.add_argument("--out-dir", type=Path, default=Path("profiles_out"), help="Output directory")
    parser.add_argument("--plot-group", help="Group of i/j points (file or comma/space list)")
    parser.add_argument("--plot-group-ij-base", type=int, default=1, choices=[0, 1], help="i/j base in group list")
    parser.add_argument("--plot-out-dir", type=Path, default=Path("profiles_plots"), help="Output directory for plots")
    args = parser.parse_args()

    k_list = read_list(args.k_list, cast=int)
    vars_list = read_vars(args.vars)
    if not k_list:
        raise SystemExit("Empty k list")
    if not vars_list:
        raise SystemExit("Empty variable list")

    points = parse_map_file(args.map)
    if not points:
        raise SystemExit(f"No points found in {args.map}")

    # Convert to 0-based
    k0 = [k - args.k_base for k in k_list]
    points0 = [(idx, i - args.ij_base, j - args.ij_base) for idx, i, j in points]

    args.out_dir.mkdir(parents=True, exist_ok=True)

    with nc.Dataset(args.fv3) as ds:
        if args.plot_group:
            group_pts = parse_ij_pairs(args.plot_group)
            group_pts0 = [(i - args.plot_group_ij_base, j - args.plot_group_ij_base) for i, j in group_pts]
            plot_profiles_for_group(
                ds=ds,
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
                f.write("# index i j\n")
                f.write(f"{idx} {i_s + 1} {j_s + 1}\n")
                f.write("# k value\n")
                for var_name, k_in, k_s in zip(vars_list, k_list, k0):
                    if var_name not in ds.variables:
                        raise SystemExit(f"Variable not found: {var_name}")
                    var = ds.variables[var_name]
                    slicer = infer_index_order(var, i_s, j_s, k_s, args.time_index)
                    val = var[slicer]
                    try:
                        val = float(val)
                    except TypeError:
                        val = float(np.asarray(val).squeeze())
                    f.write(f"{k_in} {val:.6g}\n")
            print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
