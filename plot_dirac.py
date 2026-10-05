#!/usr/bin/env python3
"""
Plot MPAS-JEDI output fields (e.g. the Dirac test, B*delta) on an MPAS mesh,
global or regional.

Follows the conventions of plot_inc.py (horizontal maps via tricontourf on
latCell/lonCell) and plot_inc_CellorEdge_level.py (level-vs-position sections),
but with locations, domain, levels and vertical coordinate as parameters.

Three plot types, selected with --plot (one or more):
  horiz    horizontal map at one or more model levels
  vsec     vertical cross section along a great-circle line
  profile  vertical profile in the single column nearest to a point

Location arguments accept either explicit coordinates or the keyword "max",
meaning the cell (and level) where |field| is largest -- i.e. the Dirac
response peak -- for the variable given by --center-var (default: first -v).

Grid fields (latCell, lonCell, latEdge, lonEdge, zgrid) are read from the data
file; if mpas-jedi did not write them there, pass an MPAS init/restart file
with --grid-file.

Optionally (--gauss LX LY LZ) a reference Gaussian
    G = exp(-0.5*[(dx/LX)^2 + (dy/LY)^2 + (dk/LZ)^2])
is overlaid, where dx/dy are east-west/north-south distances in km from the
Gaussian centre and dk the distance in model levels. Maps and sections show
dashed contours of G at the fractions --gauss-fracs (G = 0.61 at one length
scale); profiles show A*G with A the field value at the Gaussian centre.

Examples
--------
# horizontal maps of qv and theta at levels 10 and 30, auto-zoomed 6 deg
# around the Dirac peak:
python plot_dirac.py mpas.Dirac...nc -v qv theta --plot horiz \\
       --levels 10 30 --box-around max --halfwidth 6

# map at the level of the peak, explicit domain:
python plot_dirac.py mpas.Dirac...nc -v qv --plot horiz --levels max \\
       --domain -110 -80 25 50

# zonal and meridional cross sections through the peak, height axis:
python plot_dirac.py mpas.Dirac...nc -v qv theta --plot vsec \\
       --through max --orient zonal meridional --halfwidth 6 --zcoord height \\
       --grid-file x1.157859.init.nc

# cross section between two points, and a column profile at a point:
python plot_dirac.py mpas.Dirac...nc -v theta --plot vsec profile \\
       --start 35 -100 --end 40 -90 --point "37.5,-95"

# compare the Dirac response with a Gaussian of 150 km x 150 km x 3 levels:
python plot_dirac.py mpas.Dirac...nc -v qv --plot horiz vsec profile \\
       --box-around max --through max --orient zonal meridional \\
       --gauss 150 150 3
"""

import argparse
import os
import sys

import numpy as np
from netCDF4 import Dataset
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    HAVE_CARTOPY = True
except ImportError:
    HAVE_CARTOPY = False

EARTH_RADIUS_KM = 6371.0

# MPAS horizontal dimension -> (lat name, lon name)
LOCATION_COORDS = {
    'nCells': ('latCell', 'lonCell'),
    'nEdges': ('latEdge', 'lonEdge'),
    'nVertices': ('latVertex', 'lonVertex'),
}


# ---------------------------------------------------------------- arguments

def parse_point(text):
    """'max' or 'LAT,LON' / 'LAT LON' -> 'max' or (lat, lon)."""
    if text == 'max':
        return 'max'
    lat, lon = (float(x) for x in text.replace(',', ' ').split())
    return (lat, lon)


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('file', help='MPAS-JEDI output netCDF file')
    p.add_argument('-v', '--vars', nargs='+', required=True, help='variables to plot')
    p.add_argument('--plot', nargs='+', default=['horiz'],
                   choices=['horiz', 'vsec', 'profile'])
    p.add_argument('--grid-file', help='MPAS file with lat/lon/zgrid if missing from FILE')
    p.add_argument('--time', type=int, default=0, help='Time index (default 0)')
    p.add_argument('--center-var', help='variable used to locate "max" (default: first -v)')

    g = p.add_argument_group('horizontal maps')
    g.add_argument('--levels', nargs='+', default=['max'],
                   help='1-based model levels, or "max" (default) for the peak level')
    g.add_argument('--domain', nargs=4, type=float,
                   metavar=('LONMIN', 'LONMAX', 'LATMIN', 'LATMAX'),
                   help='map extent in degrees (default: whole mesh)')
    g.add_argument('--box-around', type=parse_point, metavar='max|"LAT,LON"',
                   help='center the map on a point, with --halfwidth')
    g.add_argument('--halfwidth', type=float, default=5.0,
                   help='half-width in degrees for --box-around and --through (default 5)')

    g = p.add_argument_group('vertical cross sections')
    g.add_argument('--start', nargs=2, type=float, metavar=('LAT', 'LON'))
    g.add_argument('--end', nargs=2, type=float, metavar=('LAT', 'LON'))
    g.add_argument('--through', type=parse_point, metavar='max|"LAT,LON"',
                   help='section(s) through a point instead of --start/--end')
    g.add_argument('--orient', nargs='+', default=['zonal'], choices=['zonal', 'meridional'],
                   help='orientation(s) of sections built with --through')
    g.add_argument('--npts', type=int, default=200, help='points along the section')

    g = p.add_argument_group('profiles')
    g.add_argument('--point', type=parse_point, metavar='max|"LAT,LON"',
                   help='column location for --plot profile (default: max)')

    g = p.add_argument_group('Gaussian overlay')
    g.add_argument('--gauss', nargs=3, type=float, metavar=('LX_KM', 'LY_KM', 'LZ_LEV'),
                   help='overlay a Gaussian with these east-west, north-south (km) '
                        'and vertical (model levels) length scales')
    g.add_argument('--gauss-center', type=parse_point, default='max', metavar='max|"LAT,LON"',
                   help='Gaussian centre (default: max)')
    g.add_argument('--gauss-level', type=int,
                   help='1-based centre level (default: level of the peak)')
    g.add_argument('--gauss-amp', type=float,
                   help='profile amplitude (default: field value at the centre)')
    g.add_argument('--gauss-fracs', nargs='+', type=float, default=[0.25, 0.5, 0.75],
                   help='contour values of the unit Gaussian (default 0.25 0.5 0.75)')

    g = p.add_argument_group('vertical axis, colors, output')
    g.add_argument('--zcoord', choices=['level', 'height'], default='level',
                   help='vertical axis; "height" needs zgrid (default level)')
    g.add_argument('--zrange', nargs=2, type=float, metavar=('ZMIN', 'ZMAX'),
                   help='vertical axis limits (levels, or km for height)')
    g.add_argument('--cmap', help='colormap (default RdBu_r if signed, viridis otherwise)')
    g.add_argument('--vmax', type=float, help='color limit; symmetric for signed fields')
    g.add_argument('--ncontours', type=int, default=21)
    g.add_argument('--no-coast', action='store_true',
                   help='skip coastlines/borders (offline nodes without Natural Earth data)')
    g.add_argument('--outdir', default='.')
    g.add_argument('--prefix', default='', help='prefix for output file names')
    g.add_argument('--fmt', default='png')
    g.add_argument('--dpi', type=int, default=150)
    return p.parse_args()


# ---------------------------------------------------------------- mesh

class Mesh:
    """Coordinates of an MPAS mesh, looked up first in the data file, then the grid file."""

    def __init__(self, data_nc, grid_nc):
        self._sources = [nc for nc in (data_nc, grid_nc) if nc is not None]
        self._latlon = {}
        self._trees = {}
        self.zgrid = self._read('zgrid', required=False)   # (nCells, nVertLevelsP1), m

    def _read(self, name, required=True):
        for nc in self._sources:
            if name in nc.variables:
                return np.ma.filled(nc.variables[name][:].astype(np.float64), np.nan)
        if required:
            sys.exit('ERROR: "%s" not found; pass an MPAS init/restart file with --grid-file' % name)
        return None

    def latlon(self, location):
        """Degrees, lon wrapped to [-180, 180)."""
        if location not in self._latlon:
            latname, lonname = LOCATION_COORDS[location]
            lat = np.degrees(self._read(latname))
            lon = (np.degrees(self._read(lonname)) + 180.0) % 360.0 - 180.0
            self._latlon[location] = (lat, lon)
        return self._latlon[location]

    def nearest(self, location, lats, lons):
        """Index of the nearest mesh point for each (lat, lon)."""
        lat, lon = self.latlon(location)
        if location not in self._trees:
            self._trees[location] = _make_tree(_to_xyz(lat, lon))
        return self._trees[location](_to_xyz(np.atleast_1d(lats), np.atleast_1d(lons)))

    def layer_heights_km(self, cells):
        """Mid-layer heights (km) for the given cells, shape (len(cells), nVertLevels)."""
        z = self.zgrid[cells, :]
        return 0.5 * (z[:, :-1] + z[:, 1:]) / 1000.0


def _to_xyz(lat_deg, lon_deg):
    lat, lon = np.radians(lat_deg), np.radians(lon_deg)
    return np.column_stack((np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)))


def _make_tree(xyz):
    """Return a nearest-neighbour query function (scipy if available, else brute force)."""
    try:
        from scipy.spatial import cKDTree
        tree = cKDTree(xyz)
        return lambda q: tree.query(q)[1]
    except ImportError:
        return lambda q: np.array([np.argmax(xyz @ v) for v in q])


# ---------------------------------------------------------------- field

class Field:
    """One variable at one time: values (nLoc,) or (nLoc, nLev), plus metadata."""

    def __init__(self, nc, name, time):
        if name not in nc.variables:
            sys.exit('ERROR: variable "%s" not in file' % name)
        var = nc.variables[name]
        dims = var.dimensions
        self.name = name
        self.units = getattr(var, 'units', '')
        self.long_name = getattr(var, 'long_name', name)
        self.location = next((d for d in dims if d in LOCATION_COORDS), None)
        if self.location is None:
            sys.exit('ERROR: "%s" has no nCells/nEdges/nVertices dimension' % name)
        values = var[time] if dims[0] == 'Time' else var[:]
        self.values = np.ma.filled(values.astype(np.float64), np.nan)
        self.is3d = self.values.ndim == 2
        self.nlev = self.values.shape[1] if self.is3d else 1

    def at_level(self, level):
        """level is 1-based, as in plot_inc.py."""
        return self.values[:, level - 1] if self.is3d else self.values

    def peak(self):
        """(location index, 1-based level) of max |value|."""
        flat = np.nanargmax(np.abs(self.values))
        if self.is3d:
            i, k = np.unravel_index(flat, self.values.shape)
            return int(i), int(k) + 1
        return int(flat), 1

    def label(self):
        return '%s (%s)' % (self.name, self.units) if self.units else self.name


def color_setup(values, args):
    """Contour levels, colormap and norm: symmetric diverging map for signed fields."""
    finite = values[np.isfinite(values)]
    vmin, vmax = (float(finite.min()), float(finite.max())) if finite.size else (0.0, 1.0)
    signed = vmin < 0.0 < vmax
    if signed:
        lim = args.vmax if args.vmax else max(abs(vmin), abs(vmax))
        lo, hi = -lim, lim
        cmap = args.cmap or 'RdBu_r'
    else:
        lo, hi = (vmin, args.vmax) if args.vmax else (vmin, vmax)
        cmap = args.cmap or 'viridis'
    if hi <= lo:
        hi = lo + 1e-12
    return np.linspace(lo, hi, args.ncontours), cmap, mcolors.Normalize(lo, hi)


def stats_text(values):
    return 'min=%.3e  max=%.3e' % (np.nanmin(values), np.nanmax(values))


# ---------------------------------------------------------------- geometry helpers

def great_circle(lat1, lon1, lat2, lon2, npts):
    """Points along the great circle, and along-track distance (km)."""
    a, b = _to_xyz(np.array([lat1]), np.array([lon1]))[0], _to_xyz(np.array([lat2]), np.array([lon2]))[0]
    omega = np.arccos(np.clip(a @ b, -1.0, 1.0))
    t = np.linspace(0.0, 1.0, npts)
    if omega < 1e-12:
        xyz = np.tile(a, (npts, 1))
    else:
        xyz = (np.sin((1 - t) * omega)[:, None] * a + np.sin(t * omega)[:, None] * b) / np.sin(omega)
    lats = np.degrees(np.arcsin(np.clip(xyz[:, 2], -1, 1)))
    lons = np.degrees(np.arctan2(xyz[:, 1], xyz[:, 0]))
    return lats, lons, t * omega * EARTH_RADIUS_KM


def resolve_point(point, mesh, center):
    """'max' -> the peak location of the center field; otherwise (lat, lon) unchanged."""
    if point != 'max':
        return point
    i, _ = center.peak()
    lat, lon = mesh.latlon(center.location)
    return float(lat[i]), float(lon[i])


def sections_from_args(args, mesh, center):
    """List of (tag, start, end) for the requested cross sections."""
    if args.start and args.end:
        return [('', tuple(args.start), tuple(args.end))]
    lat0, lon0 = resolve_point(args.through or 'max', mesh, center)
    w = args.halfwidth
    sections = []
    for orient in args.orient:
        if orient == 'zonal':
            sections.append(('zonal', (lat0, lon0 - w), (lat0, lon0 + w)))
        else:
            sections.append(('meridional', (lat0 - w, lon0), (lat0 + w, lon0)))
    return sections


# ---------------------------------------------------------------- Gaussian reference

def gaussian(lat, lon, level, lat0, lon0, level0, lx_km, ly_km, lz_lev):
    """Unit-amplitude Gaussian exp(-0.5*[(dx/Lx)^2 + (dy/Ly)^2 + (dk/Lz)^2]).

    dx, dy: east-west and north-south distances (km) from (lat0, lon0) on the
    local tangent plane; dk = level - level0 in model levels. The value at one
    length scale from the centre is exp(-0.5) = 0.61. lat, lon, level broadcast.
    """
    dlon = (np.asarray(lon) - lon0 + 180.0) % 360.0 - 180.0
    dx = EARTH_RADIUS_KM * np.cos(np.radians(lat0)) * np.radians(dlon)
    dy = EARTH_RADIUS_KM * np.radians(np.asarray(lat) - lat0)
    dk = np.asarray(level, dtype=float) - level0
    return np.exp(-0.5 * ((dx / lx_km) ** 2 + (dy / ly_km) ** 2 + (dk / lz_lev) ** 2))


class GaussOverlay:
    """The --gauss reference: centre, length scales, and how it is drawn."""

    def __init__(self, args, mesh, center):
        self.lx, self.ly, self.lz = args.gauss
        self.lat0, self.lon0 = resolve_point(args.gauss_center, mesh, center)
        self.level0 = args.gauss_level or center.peak()[1]
        self.fracs = sorted(args.gauss_fracs)
        self.amp = args.gauss_amp
        self.mesh = mesh
        print('Gaussian: centre (%.3fN, %.3fE) level %d, L = %g km x %g km x %g levels'
              % (self.lat0, self.lon0, self.level0, self.lx, self.ly, self.lz))

    def __call__(self, lat, lon, level):
        return gaussian(lat, lon, level, self.lat0, self.lon0, self.level0,
                        self.lx, self.ly, self.lz)

    def amplitude(self, field):
        """--gauss-amp if given, else the field value at the Gaussian centre."""
        if self.amp is not None:
            return self.amp
        i = int(self.mesh.nearest(field.location, self.lat0, self.lon0)[0])
        return field.values[i, self.level0 - 1] if field.is3d else field.values[i]

    def contour(self, ax, x, y, g, triangulated=False, **kw):
        """Dashed contours of the unit Gaussian g at --gauss-fracs, labelled by fraction."""
        draw = ax.tricontour if triangulated else ax.contour
        cs = draw(x, y, g, levels=self.fracs, colors='k', linestyles='--',
                  linewidths=1.0, **kw)
        ax.clabel(cs, fmt='%.2g', fontsize=7)


def outname(args, *parts):
    stem = '_'.join(str(p) for p in parts if p != '')
    return os.path.join(args.outdir, '%s%s.%s' % (args.prefix, stem, args.fmt))


def vertical_axis(args, mesh, cells, nlev):
    """Vertical coordinate for the columns `cells`: (Z array (len(cells), nlev), label)."""
    if args.zcoord == 'height':
        if mesh.zgrid is None:
            print('WARNING: zgrid not found (use --grid-file); falling back to model level')
        else:
            return mesh.layer_heights_km(cells), 'Height (km)'
    return np.tile(np.arange(1, nlev + 1, dtype=float), (len(cells), 1)), 'Model level'


# ---------------------------------------------------------------- plots

def map_axes(fig, args, extent):
    if HAVE_CARTOPY:
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
        ax.set_extent(extent, crs=ccrs.PlateCarree())
        if not args.no_coast:
            ax.coastlines(resolution='50m', linewidth=0.6)
            ax.add_feature(cfeature.BORDERS, linewidth=0.4)
            ax.add_feature(cfeature.STATES, linewidth=0.3, edgecolor='gray')
        gl = ax.gridlines(draw_labels=True, linestyle='--', linewidth=0.4)
        gl.top_labels = gl.right_labels = False
        return ax, {'transform': ccrs.PlateCarree()}
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_aspect('equal')
    return ax, {}


def map_extent(args, mesh, field, center):
    if args.domain:
        return list(args.domain)
    if args.box_around:
        lat0, lon0 = resolve_point(args.box_around, mesh, center)
        w = args.halfwidth
        return [lon0 - w, lon0 + w, lat0 - w, lat0 + w]
    lat, lon = mesh.latlon(field.location)
    return [lon.min(), lon.max(), lat.min(), lat.max()]


def plot_horizontal(args, mesh, field, center, sections, gauss):
    lat, lon = mesh.latlon(field.location)
    extent = map_extent(args, mesh, field, center)
    # Triangulate only the points in (a margin around) the domain: faster, and
    # avoids triangles spanning the concave edges of a regional mesh.
    margin = 0.1 * max(extent[1] - extent[0], extent[3] - extent[2]) + 0.5
    inside = ((lon >= extent[0] - margin) & (lon <= extent[1] + margin) &
              (lat >= extent[2] - margin) & (lat <= extent[3] + margin))

    levels = []
    for lv in args.levels:
        levels.append(field.peak()[1] if lv == 'max' else int(lv))
    if not field.is3d:
        levels = [1]

    for level in sorted(set(levels)):
        if not 1 <= level <= field.nlev:
            print('WARNING: %s level %d out of range 1..%d, skipped' % (field.name, level, field.nlev))
            continue
        data = field.at_level(level)
        ok = inside & np.isfinite(data)
        clev, cmap, norm = color_setup(data[ok], args)

        fig = plt.figure(figsize=(9, 8))
        ax, tr = map_axes(fig, args, extent)
        cs = ax.tricontourf(lon[ok], lat[ok], data[ok], levels=clev, cmap=cmap,
                            norm=norm, extend='both', **tr)
        if gauss:
            g = gauss(lat[ok], lon[ok], level if field.is3d else gauss.level0)
            gauss.contour(ax, lon[ok], lat[ok], g, triangulated=True, **tr)
        for tag, (la1, lo1), (la2, lo2) in sections:
            slat, slon, _ = great_circle(la1, lo1, la2, lo2, 50)
            ax.plot(slon, slat, 'k--', linewidth=1.0, **tr)
        if args.box_around == 'max' or args.through == 'max':
            plat, plon = resolve_point('max', mesh, center)
            ax.plot(plon, plat, 'kx', markersize=8, **tr)
        fig.colorbar(cs, ax=ax, orientation='horizontal', pad=0.06, shrink=0.85,
                     label=field.label())
        lev_txt = ', level %d' % level if field.is3d else ''
        ax.set_title('%s%s\n%s' % (field.long_name, lev_txt, stats_text(data[ok])))
        fname = outname(args, 'horiz', field.name, 'L%02d' % level if field.is3d else '')
        fig.savefig(fname, dpi=args.dpi, bbox_inches='tight')
        plt.close(fig)
        print('wrote', fname)


def plot_vsection(args, mesh, field, tag, start, end, gauss):
    slat, slon, dist = great_circle(start[0], start[1], end[0], end[1], args.npts)
    idx = mesh.nearest(field.location, slat, slon)
    data = field.values[idx, :]                                   # (npts, nlev)
    cells = mesh.nearest('nCells', slat, slon) if args.zcoord == 'height' else idx
    zz, zlabel = vertical_axis(args, mesh, cells, field.nlev)
    xx = np.tile(dist[:, None], (1, field.nlev))
    clev, cmap, norm = color_setup(data, args)

    fig, ax = plt.subplots(figsize=(10, 6))
    cs = ax.contourf(xx, zz, data, levels=clev, cmap=cmap, norm=norm, extend='both')
    if gauss:
        levels = np.arange(1, field.nlev + 1)
        gauss.contour(ax, xx, zz, gauss(slat[:, None], slon[:, None], levels[None, :]))
    if args.zcoord == 'height' and mesh.zgrid is not None:
        terrain = mesh.zgrid[cells, 0] / 1000.0
        ax.fill_between(dist, 0.0, terrain, color='saddlebrown', zorder=3)
    if args.zrange:
        ax.set_ylim(args.zrange)
    ax.set_ylabel(zlabel)
    ax.set_xlabel('Distance (km) from (%.2fN, %.2fE) to (%.2fN, %.2fE)'
                  % (start[0], start[1], end[0], end[1]))
    fig.colorbar(cs, ax=ax, orientation='horizontal', pad=0.12, label=field.label())
    ax.set_title('%s %s cross section\n%s' % (field.long_name, tag, stats_text(data)))
    fname = outname(args, 'vsec', field.name, tag)
    fig.savefig(fname, dpi=args.dpi, bbox_inches='tight')
    plt.close(fig)
    print('wrote', fname)


def plot_profile(args, mesh, field, point, gauss):
    i = int(mesh.nearest(field.location, point[0], point[1])[0])
    cell = int(mesh.nearest('nCells', point[0], point[1])[0]) if args.zcoord == 'height' else i
    zz, zlabel = vertical_axis(args, mesh, [cell], field.nlev)
    lat, lon = mesh.latlon(field.location)

    fig, ax = plt.subplots(figsize=(5, 7))
    ax.plot(field.values[i, :], zz[0], 'o-', markersize=3, label=field.name)
    if gauss:
        g = gauss(lat[i], lon[i], np.arange(1, field.nlev + 1))
        ax.plot(gauss.amplitude(field) * g, zz[0], 'k--', label='Gaussian')
        ax.legend(fontsize=8)
    ax.axvline(0.0, color='gray', linewidth=0.6)
    if args.zrange:
        ax.set_ylim(args.zrange)
    ax.set_xlabel(field.label())
    ax.set_ylabel(zlabel)
    ax.grid(True, linestyle=':')
    ax.set_title('%s\n%s index %d at (%.2fN, %.2fE)'
                 % (field.long_name, field.location, i, lat[i], lon[i]))
    fname = outname(args, 'profile', field.name, '%.2fN_%.2fE' % (lat[i], lon[i]))
    fig.savefig(fname, dpi=args.dpi, bbox_inches='tight')
    plt.close(fig)
    print('wrote', fname)


# ---------------------------------------------------------------- main

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    if not HAVE_CARTOPY and 'horiz' in args.plot:
        print('WARNING: cartopy not available; maps drawn without projection or coastlines')

    data_nc = Dataset(args.file, 'r')
    grid_nc = Dataset(args.grid_file, 'r') if args.grid_file else None
    mesh = Mesh(data_nc, grid_nc)

    center = Field(data_nc, args.center_var or args.vars[0], args.time)
    ci, ck = center.peak()
    clat, clon = mesh.latlon(center.location)
    print('peak |%s| = %.4e at %s %d (%.3fN, %.3fE), level %d'
          % (center.name, abs(center.values[ci] if not center.is3d else center.values[ci, ck - 1]),
             center.location, ci, clat[ci], clon[ci], ck))

    sections = sections_from_args(args, mesh, center) if 'vsec' in args.plot else []
    gauss = GaussOverlay(args, mesh, center) if args.gauss else None

    for name in args.vars:
        field = center if name == center.name else Field(data_nc, name, args.time)
        if 'horiz' in args.plot:
            plot_horizontal(args, mesh, field, center, sections, gauss)
        if not field.is3d and ('vsec' in args.plot or 'profile' in args.plot):
            print('NOTE: %s is 2D; vertical plots skipped' % name)
            continue
        for tag, start, end in sections:
            plot_vsection(args, mesh, field, tag, start, end, gauss)
        if 'profile' in args.plot:
            plot_profile(args, mesh, field, resolve_point(args.point or 'max', mesh, center), gauss)

    data_nc.close()
    if grid_nc is not None:
        grid_nc.close()


if __name__ == '__main__':
    main()
