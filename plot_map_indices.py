"""
Plot model grid index points from a map file.

The file format is expected to have 3 header lines, then data lines with
at least 5 columns where columns 4 and 5 are model_i and model_j (1-based).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def read_model_indices(path: Path, skip_lines: int) -> tuple[np.ndarray, np.ndarray]:
    i_list = []
    j_list = []
    for line_no, line in enumerate(path.read_text().splitlines(), start=1):
        if line_no <= skip_lines:
            continue
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 5:
            continue
        i_list.append(int(parts[3]))
        j_list.append(int(parts[4]))
    if not i_list:
        raise SystemExit(f"No model indices found in {path}")
    return np.asarray(i_list), np.asarray(j_list)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot model index points from a map file.")
    parser.add_argument("--map", type=Path, default=Path("test-filter_to_model_map_indices.txt"), help="Map file to read")
    parser.add_argument("--skip-lines", type=int, default=3, help="Number of header lines to skip")
    parser.add_argument("--nx", type=int, default=3950, help="Grid size in i direction")
    parser.add_argument("--ny", type=int, default=2700, help="Grid size in j direction")
    parser.add_argument("--output", type=Path, default=Path("model_index_points.png"), help="Output image")
    parser.add_argument("--point-size", type=float, default=3.0, help="Marker size")
    args = parser.parse_args()

    i_idx, j_idx = read_model_indices(args.map, args.skip_lines)

    try:
        import matplotlib.pyplot as plt
    except ImportError as e:  # pragma: no cover
        raise SystemExit("matplotlib is required for plotting") from e

    fig, ax = plt.subplots(figsize=(9, 6), dpi=150)
    ax.scatter(i_idx, j_idx, s=args.point_size, c="tab:blue", alpha=0.7, rasterized=True)
    ax.set_xlim(1, args.nx)
    ax.set_ylim(1, args.ny)
    ax.set_xlabel("Model i")
    ax.set_ylabel("Model j")
    ax.set_title("Selected model grid points")
    ax.grid(True, linewidth=0.3, alpha=0.4)
    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output)
    plt.close(fig)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
