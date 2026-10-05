#!/usr/bin/env python3
"""
Line plots of the MPAS-JEDI Dirac test output (B applied to a unit impulse)
through the impulse, with an optional reference Gaussian.

Companion to plot_dirac.py (which draws contours): it imports that script's
mesh, field and Gaussian helpers, so keep the two files in the same directory.

The impulse location is given with --center LAT LON LEVEL, the same values as
dirLats / dirLons / dirLevs in the Dirac YAML (LEVEL is 1-based, as dirLevs).

One figure per variable, three panels:
  1  x line   field at the impulse level along the impulse latitude, vs signed
              east-west distance from the impulse (km)
  2  y line   same along the impulse longitude, vs north-south distance (km)
  3  profile  field in the impulse column vs model level (or height)
Panels 1 and 2 show the actual mesh-point values (each cell the line crosses,
plotted once at its distance from the impulse) and share the same y scale.

With --gauss LX LY LZ the reference Gaussian of plot_dirac.py,
    G = exp(-0.5*[(dx/LX)^2 + (dy/LY)^2 + (dk/LZ)^2]),
is drawn as a dashed curve on each panel with peak value 1 (the response of a
normalized B to a unit impulse); dotted lines mark +-one length scale, where
G = exp(-0.5) = 0.61.

Examples
--------
# Dirac at dirLats=35.0165, dirLons=269.9974, dirLevs=10, dirVars=eastward_wind:
python plot_dirac_lines.py mpas.Dirac...nc -v uReconstructZonal \\
       --center 35.0165443420410160 269.9974 10 --gauss 150 150 3 \\
       --grid-file conus.invariant.nc

# several variables, +-300 km, lower 30 levels on the profile:
python plot_dirac_lines.py mpas.Dirac...nc -v uReconstructZonal theta qv \\
       --center 35.0165443420410160 269.9974 10 --gauss 150 150 3 \\
       --halfwidth 300 --zrange 1 30 --grid-file conus.invariant.nc
"""

import argparse
import os

import numpy as np
from netCDF4 import Dataset

from plot_dirac import (KM_PER_DEG, Center, Field, Mesh, gaussian, report_center,
                        vertical_axis, wrap_lon)
import matplotlib.pyplot as plt   # after plot_dirac, which selects the Agg backend


# ---------------------------------------------------------------- arguments

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('file', help='MPAS-JEDI Dirac output netCDF file')
    p.add_argument('-v', '--vars', nargs='+', required=True, help='variables to plot')
    p.add_argument('--center', nargs=3, type=float, required=True,
                   metavar=('LAT', 'LON', 'LEVEL'),
                   help='impulse location: dirLats, dirLons, dirLevs (1-based) from the YAML')
    p.add_argument('--halfwidth', type=float, default=500.0,
                   help='x/y line half-length in km around the impulse (default 500)')
    p.add_argument('--grid-file', help='MPAS invariant/static file with lat/lon/zgrid')
    p.add_argument('--time', type=int, default=0, help='Time index (default 0)')
    p.add_argument('--npts', type=int, default=401,
                   help='sampling points used to find the cells a line crosses (default 401)')

    g = p.add_argument_group('vertical profile')
    g.add_argument('--zcoord', choices=['level', 'height'], default='level',
                   help='vertical axis; "height" needs zgrid (default level)')
    g.add_argument('--zrange', nargs=2, type=float, metavar=('ZMIN', 'ZMAX'),
                   help='vertical axis limits (levels, or km for height)')

    g = p.add_argument_group('Gaussian reference')
    g.add_argument('--gauss', nargs=3, type=float, metavar=('LX_KM', 'LY_KM', 'LZ_LEV'),
                   help='east-west, north-south (km) and vertical (model levels) length '
                        'scales of a peak-1 reference Gaussian')

    g = p.add_argument_group('output')
    g.add_argument('--outdir', default='.')
    g.add_argument('--prefix', default='', help='prefix for output file names')
    g.add_argument('--fmt', default='png')
    g.add_argument('--dpi', type=int, default=150)
    return p.parse_args()


# ---------------------------------------------------------------- data along the lines

def line_data(mesh, field, center, axis, npts):
    """Values at the impulse level of the mesh points the x or y line passes
    through, each once, with their signed distance (km) from the impulse along
    the line; sorted by distance."""
    lats, lons, _ = center.section(axis, npts)
    cells = np.unique(mesh.nearest(field.location, lats, lons))
    lat, lon = mesh.latlon(field.location)
    if axis == 'x':
        dist = KM_PER_DEG * np.cos(np.radians(center.lat)) * wrap_lon(lon[cells] - center.lon)
    else:
        dist = KM_PER_DEG * (lat[cells] - center.lat)
    order = np.argsort(dist)
    return dist[order], field.values[cells[order], center.level - 1]


def unit_gaussian(args, center, lat, lon, level):
    """The --gauss reference (unit amplitude) centred on the impulse."""
    lx, ly, lz = args.gauss
    return gaussian(lat, lon, level, center.lat, center.lon, center.level, lx, ly, lz)


# ---------------------------------------------------------------- panels

def draw_line(ax, args, mesh, field, center, axis):
    dist, values = line_data(mesh, field, center, axis, args.npts)
    ax.plot(dist, values, 'o-', markersize=3, label=field.name)
    if args.gauss:
        lats, lons, d = center.section(axis, 1001)
        length = args.gauss[0] if axis == 'x' else args.gauss[1]
        ax.plot(d, unit_gaussian(args, center, lats, lons, center.level), 'k--',
                label='Gaussian, L = %g km' % length)
        for s in (-length, length):
            ax.axvline(s, color='gray', linestyle=':', linewidth=0.8)
    ax.axhline(0.0, color='gray', linewidth=0.6)
    ax.axvline(0.0, color='k', linestyle=':', linewidth=0.8)
    ax.set_xlim(-center.hw, center.hw)
    if axis == 'x':
        ax.set_xlabel('x: distance from impulse (km), west -> east')
        ax.set_title('x line at lat %.2fN, level %d' % (center.lat, center.level), fontsize=10)
    else:
        ax.set_xlabel('y: distance from impulse (km), south -> north')
        ax.set_title('y line at lon %.2fE, level %d' % (center.lon, center.level), fontsize=10)
    ax.set_ylabel(field.label())
    ax.grid(True, linestyle=':')
    ax.legend(fontsize=8)


def draw_profile(ax, args, mesh, field, center, i_center):
    cell = (int(mesh.nearest('nCells', center.lat, center.lon)[0])
            if args.zcoord == 'height' else i_center)
    zz, zlabel = vertical_axis(args, mesh, [cell], field.nlev)
    z = zz[0]
    model_levels = np.arange(1, field.nlev + 1)
    ax.plot(field.values[i_center, :], z, 'o-', markersize=3, label=field.name)
    if args.gauss:
        lz = args.gauss[2]
        k = np.linspace(1, field.nlev, 1000)
        ax.plot(unit_gaussian(args, center, center.lat, center.lon, k),
                np.interp(k, model_levels, z), 'k--', label='Gaussian, L = %g levels' % lz)
        for kk in (center.level - lz, center.level + lz):
            if 1 <= kk <= field.nlev:
                ax.axhline(np.interp(kk, model_levels, z), color='gray', linestyle=':',
                           linewidth=0.8)
    ax.axvline(0.0, color='gray', linewidth=0.6)
    ax.axhline(z[center.level - 1], color='k', linestyle=':', linewidth=0.8)
    if args.zrange:
        ax.set_ylim(args.zrange)
    ax.set_xlabel(field.label())
    ax.set_ylabel(zlabel)
    ax.set_title('profile at (%.2fN, %.2fE)' % (center.lat, center.lon), fontsize=10)
    ax.grid(True, linestyle=':')
    ax.legend(fontsize=8)


def plot_variable(args, mesh, field, center, i_center):
    """One figure: x line, y line, vertical profile."""
    fig, (ax_x, ax_y, ax_z) = plt.subplots(
        1, 3, figsize=(17, 5.5), gridspec_kw={'width_ratios': [1.3, 1.3, 1.0], 'wspace': 0.3})
    draw_line(ax_x, args, mesh, field, center, 'x')
    draw_line(ax_y, args, mesh, field, center, 'y')
    draw_profile(ax_z, args, mesh, field, center, i_center)

    # Same y scale on the x and y lines, so their widths compare directly.
    ylo = min(ax_x.get_ylim()[0], ax_y.get_ylim()[0])
    yhi = max(ax_x.get_ylim()[1], ax_y.get_ylim()[1])
    ax_x.set_ylim(ylo, yhi)
    ax_y.set_ylim(ylo, yhi)

    gtxt = ('\nGaussian (peak 1) L = %g km x %g km x %g levels'
            % tuple(args.gauss) if args.gauss else '')
    fig.suptitle('%s, impulse at (%.2fN, %.2fE) level %d%s'
                 % (field.long_name, center.lat, center.lon, center.level, gtxt))
    fname = os.path.join(args.outdir, '%sdirac_lines_%s.%s' % (args.prefix, field.name, args.fmt))
    fig.savefig(fname, dpi=args.dpi, bbox_inches='tight')
    plt.close(fig)
    print('wrote', fname)


# ---------------------------------------------------------------- main

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    data_nc = Dataset(args.file, 'r')
    grid_nc = Dataset(args.grid_file, 'r') if args.grid_file else None
    mesh = Mesh(data_nc, grid_nc)
    center = Center(*args.center, halfwidth_km=args.halfwidth)

    for name in args.vars:
        field = Field(data_nc, name, args.time)
        if not 1 <= center.level <= field.nlev:
            raise SystemExit('ERROR: --center level %d outside 1..%d for %s'
                             % (center.level, field.nlev, name))
        i_center = int(mesh.nearest(field.location, center.lat, center.lon)[0])
        report_center(mesh, field, center, i_center)
        plot_variable(args, mesh, field, center, i_center)

    data_nc.close()
    if grid_nc is not None:
        grid_nc.close()


if __name__ == '__main__':
    main()
