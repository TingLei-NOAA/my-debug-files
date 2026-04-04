"""
Find sampled points near a strip in model-index space and report their extracted profile maxima.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def read_map(path: Path) -> list[tuple[int, int, int]]:
    rows: list[tuple[int, int, int]] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 5:
            continue
        idx = int(parts[0])
        model_i = int(parts[3])
        model_j = int(parts[4])
        rows.append((idx, model_i, model_j))
    return rows


def read_profile_max(path: Path) -> float:
    values = []
    in_complete = False
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith("#"):
            if "complete profile" in line.lower():
                in_complete = True
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        if in_complete:
            values.append(float(parts[1]))
    if not values:
        raise SystemExit(f"No complete-profile values found in {path}")
    return max(values)


def main() -> None:
    parser = argparse.ArgumentParser(description="Find sampled profile maxima near a strip of model indices.")
    parser.add_argument("--map", type=Path, required=True, help="filter_to_model_map_indices.txt")
    parser.add_argument("--profiles-dir", type=Path, required=True, help="Directory containing profile_subdomain_*.txt")
    parser.add_argument("--i-min", type=int, required=True, help="Minimum 1-based model i")
    parser.add_argument("--i-max", type=int, required=True, help="Maximum 1-based model i")
    parser.add_argument("--j-min", type=int, default=1, help="Minimum 1-based model j")
    parser.add_argument("--j-max", type=int, default=10**9, help="Maximum 1-based model j")
    parser.add_argument("--top-n", type=int, default=20, help="Number of top sampled points to report")
    args = parser.parse_args()

    rows = read_map(args.map)
    matches = []
    for idx, model_i, model_j in rows:
        if not (args.i_min <= model_i <= args.i_max and args.j_min <= model_j <= args.j_max):
            continue
        profile = args.profiles_dir / f"profile_subdomain_{idx:04d}.txt"
        if not profile.exists():
            continue
        prof_max = read_profile_max(profile)
        matches.append((prof_max, idx, model_i, model_j, profile))

    if not matches:
        print("No sampled points found in the requested strip.")
        return

    matches.sort(reverse=True, key=lambda x: x[0])
    print(f"Found {len(matches)} sampled points in strip: i=[{args.i_min},{args.i_max}], j=[{args.j_min},{args.j_max}]")
    print("# rank sampled_profile_max subdomain_index model_i model_j profile_file")
    for rank, (prof_max, idx, model_i, model_j, profile) in enumerate(matches[: args.top_n], start=1):
        print(f"{rank} {prof_max:.10g} {idx} {model_i} {model_j} {profile}")


if __name__ == "__main__":
    main()
