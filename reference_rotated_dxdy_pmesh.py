"""
Generate reference dx/dy for the rotated_lonlat StructuredColumns grid and plot with pcolormesh
to avoid contour triangulation masking.

Grid:
  xspace: linear, N=308, start=-60.0310278618037, end=60.0310278618037
  yspace: linear, N=308, start=-36.6574659140194, end=36.6574659140194
  projection: rotated_lonlat, north_pole=[67.4999923706, 34.9999966097]

Outputs (default dir: dr-figures):
  reference_dx_km_pmesh.png
  reference_dy_km_pmesh.png
  reference_dx_over_dy_pmesh.png
  reference_lon_pmesh.png
  reference_lat_pmesh.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np


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
    lon_unw = np.degrees(np.unwrap(np.radians(lon), axis=1))
    center = np.median(lon_unw[:, 0])
    lon_wrap = ((lon_unw - center + 180.0) % 360.0) - 180.0 + center
    return lon_wrap


def plot_pmesh(field: np.ndarray, lon: np.ndarray, lat: np.ndarray, title: str, outfile: Path, cmap: str, units: str):
    fig = plt.figure(figsize=(10, 8))
    ax = plt.axes(projection=ccrs.PlateCarree())
    lon_plot = _unwrap_lon_for_plot(lon)
    pm = ax.pcolormesh(lon_plot, lat, field, cmap=cmap, shading="auto", transform=ccrs.PlateCarree())
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


def main():
    parser = argparse.ArgumentParser(description="Reference dx/dy with pcolormesh plotting")
    parser.add_argument("--output-dir", type=Path, default=Path("dr-figures"), help="Output directory")
    args = parser.parse_args()

    nx = 308
    ny = 308
    x_start = -60.0310278618037
    x_end = 60.0310278618037
    y_start = -36.6574659140194
    y_end = 36.6574659140194
    np_lon = 67.4999923706
    np_lat = 34.9999966097

    x = np.linspace(x_start, x_end, nx)
    y = np.linspace(y_start, y_end, ny)
    xx, yy = np.meshgrid(x, y)

    rotated_crs = ccrs.RotatedPole(pole_longitude=np_lon, pole_latitude=np_lat)
    geodetic = ccrs.Geodetic().transform_points(rotated_crs, xx, yy)
    lon = geodetic[..., 0]
    lat = geodetic[..., 1]

    dx, dy, ratio = compute_dx_dy(lon, lat)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    plot_pmesh(dx / 1000.0, lon, lat, "Reference dx (km)", args.output_dir / "reference_dx_km_pmesh.png", cmap="magma", units="km")
    plot_pmesh(dy / 1000.0, lon, lat, "Reference dy (km)", args.output_dir / "reference_dy_km_pmesh.png", cmap="magma", units="km")
    plot_pmesh(ratio, lon, lat, "Reference dx/dy", args.output_dir / "reference_dx_over_dy_pmesh.png", cmap="coolwarm", units="ratio")
    plot_pmesh(lon, lon, lat, "Reference lon", args.output_dir / "reference_lon_pmesh.png", cmap="coolwarm", units="deg")
    plot_pmesh(lat, lon, lat, "Reference lat", args.output_dir / "reference_lat_pmesh.png", cmap="coolwarm", units="deg")


if __name__ == "__main__":
    main()
