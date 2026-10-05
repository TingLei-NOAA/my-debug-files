#!/usr/bin/env python3
"""
Plot the MPAS-JEDI Dirac test output (B applied to a unit impulse) in a
subdomain centred on the impulse, regional or global mesh.

Follows the conventions of plot_inc.py (horizontal maps via tricontourf on
latCell/lonCell) and plot_inc_CellorEdge_level.py (level-vs-position sections).

The impulse location is given with --center LAT LON LEVEL, the same values as
dirLats / dirLons / dirLevs in the Dirac YAML (LEVEL is the 1-based model level,
1 = lowest, exactly as dirLevs).

One figure per variable, with three panels sharing one color scale:
  left         horizontal map at the impulse level, +-HALFWIDTH km around it
  right top    x (east-west) vertical section through the impulse
  right bottom y (north-south) vertical section through the impulse
Section x axes are the signed distance from the impulse in km (-HW..+HW);
y axes are model level (default) or height (--zcoord height).

Grid fields (latCell, lonCell, latEdge, lonEdge, zgrid) are looked up first in
the Dirac file, then in --grid-file (the MPAS invariant/static file named in
the "invariant" stream of streams.atmosphere).

Optionally (--gauss LX LY LZ) a reference Gaussian centred on the impulse
    G = exp(-0.5*[(dx/LX)^2 + (dy/LY)^2 + (dk/LZ)^2])
is overlaid on all three panels as dashed contours at --gauss-fracs, where
dx/dy are east-west/north-south distances (km) and dk the distance in model
levels (G = 0.61 at one length scale).

Examples
--------
# Dirac at dirLats=37.5, dirLons=-95.0, dirLevs=30, +-500 km box:
python plot_dirac.py mpas.Dirac...nc -v qv theta --center 37.5 -95.0 30 \\
       --grid-file conus.invariant.nc

# zoom to +-300 km, lower 45 levels only, height axis:
python plot_dirac.py mpas.Dirac...nc -v qv --center 37.5 -95.0 30 \\
       --halfwidth 300 --zrange 1 45 --grid-file conus.invariant.nc

# compare with a Gaussian of 150 km x 150 km x 3 levels:
python plot_dirac.py mpas.Dirac...nc -v qv --center 37.5 -95.0 30 \\
       --gauss 150 150 3 --grid-file conus.invariant.nc --no-coast
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
KM_PER_DEG = np.pi * EARTH_RADIUS_KM / 180.0

# MPAS horizontal dimension -> (lat name, lon name)
LOCATION_COORDS = {
    'nCells': ('latCell', 'lonCell'),
    'nEdges': ('latEdge', 'lonEdge'),
    'nVertices': ('latVertex', 'lonVertex'),
}


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
                   help='subdomain half-width in km around the impulse (default 500)')
    p.add_argument('--grid-file', help='MPAS invariant/static file with lat/lon/zgrid')
    p.add_argument('--time', type=int, default=0, help='Time index (default 0)')

    g = p.add_argument_group('vertical sections')
    g.add_argument('--npts', type=int, default=201, help='points along a section')
    g.add_argument('--zcoord', choices=['level', 'height'], default='level',
                   help='vertical axis; "height" needs zgrid (default level)')
    g.add_argument('--zrange', nargs=2, type=float, metavar=('ZMIN', 'ZMAX'),
                   help='vertical axis limits (levels, or km for height)')

    g = p.add_argument_group('Gaussian overlay')
    g.add_argument('--gauss', nargs=3, type=float, metavar=('LX_KM', 'LY_KM', 'LZ_LEV'),
                   help='overlay a Gaussian centred on the impulse with these east-west, '
                        'north-south (km) and vertical (model levels) length scales')
    g.add_argument('--gauss-fracs', nargs='+', type=float, default=[0.25, 0.5, 0.75],
                   help='contour values of the unit Gaussian (default 0.25 0.5 0.75)')

    g = p.add_argument_group('colors and output')
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
            sys.exit('ERROR: "%s" not found; pass the MPAS invariant file with --grid-file' % name)
        return None

    def latlon(self, location):
        """Degrees, lon wrapped to [-180, 180)."""
        if location not in self._latlon:
            latname, lonname = LOCATION_COORDS[location]
            lat = np.degrees(self._read(latname))
            lon = wrap_lon(np.degrees(self._read(lonname)))
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


def wrap_lon(lon):
    return (np.asarray(lon) + 180.0) % 360.0 - 180.0


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
        if self.values.ndim != 2:
            sys.exit('ERROR: "%s" is not a 3D (location x level) field' % name)
        self.nlev = self.values.shape[1]

    def at_level(self, level):
        """level is 1-based, as dirLevs and plot_inc.py."""
        return self.values[:, level - 1]

    def at(self, i, level):
        return self.values[i, level - 1]

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


# ---------------------------------------------------------------- impulse-centred geometry

class Center:
    """The Dirac impulse location and the subdomain around it."""

    def __init__(self, lat, lon, level, halfwidth_km):
        self.lat, self.lon, self.level = lat, float(wrap_lon(lon)), int(level)
        self.hw = halfwidth_km

    def extent(self):
        """[lonmin, lonmax, latmin, latmax] of the +-halfwidth km box."""
        dlat = self.hw / KM_PER_DEG
        dlon = self.hw / (KM_PER_DEG * np.cos(np.radians(self.lat)))
        return [self.lon - dlon, self.lon + dlon, self.lat - dlat, self.lat + dlat]

    def section(self, axis, npts):
        """Points on the x (constant lat) or y (constant lon) line through the
        impulse, and their signed distance (km) from it, -hw..+hw."""
        dist = np.linspace(-self.hw, self.hw, npts)
        if axis == 'x':
            lats = np.full(npts, self.lat)
            lons = self.lon + dist / (KM_PER_DEG * np.cos(np.radians(self.lat)))
        else:
            lats = self.lat + dist / KM_PER_DEG
            lons = np.full(npts, self.lon)
        return lats, wrap_lon(lons), dist


# ---------------------------------------------------------------- Gaussian reference

def gaussian(lat, lon, level, lat0, lon0, level0, lx_km, ly_km, lz_lev):
    """Unit-amplitude Gaussian exp(-0.5*[(dx/Lx)^2 + (dy/Ly)^2 + (dk/Lz)^2]).

    dx, dy: east-west and north-south distances (km) from (lat0, lon0) on the
    local tangent plane; dk = level - level0 in model levels. The value at one
    length scale from the centre is exp(-0.5) = 0.61. lat, lon, level broadcast.
    """
    dx = KM_PER_DEG * np.cos(np.radians(lat0)) * wrap_lon(np.asarray(lon) - lon0)
    dy = KM_PER_DEG * (np.asarray(lat) - lat0)
    dk = np.asarray(level, dtype=float) - level0
    return np.exp(-0.5 * ((dx / lx_km) ** 2 + (dy / ly_km) ** 2 + (dk / lz_lev) ** 2))


class GaussOverlay:
    """The --gauss reference centred on the impulse, and how it is drawn."""

    def __init__(self, args, center):
        self.lx, self.ly, self.lz = args.gauss
        self.center = center
        self.fracs = sorted(args.gauss_fracs)
        print('Gaussian: L = %g km (E-W) x %g km (N-S) x %g levels' % (self.lx, self.ly, self.lz))

    def __call__(self, lat, lon, level):
        c = self.center
        return gaussian(lat, lon, level, c.lat, c.lon, c.level, self.lx, self.ly, self.lz)

    def contour(self, ax, x, y, g, triangulated=False, **kw):
        """Dashed contours of the unit Gaussian g at --gauss-fracs, labelled by fraction."""
        draw = ax.tricontour if triangulated else ax.contour
        cs = draw(x, y, g, levels=self.fracs, colors='k', linestyles='--',
                  linewidths=1.0, **kw)
        ax.clabel(cs, fmt='%.2g', fontsize=7)


# ---------------------------------------------------------------- panel data

class MapData:
    """Field at the impulse level, for the mesh points inside the subdomain."""

    def __init__(self, mesh, field, center):
        lat, lon = mesh.latlon(field.location)
        extent = center.extent()
        # Keep only the points in (a margin around) the subdomain: faster, and
        # avoids triangles spanning the concave edges of a regional mesh.
        margin = 0.1 * (extent[3] - extent[2])
        data = field.at_level(center.level)
        ok = ((lon >= extent[0] - margin) & (lon <= extent[1] + margin) &
              (lat >= extent[2] - margin) & (lat <= extent[3] + margin) & np.isfinite(data))
        if ok.sum() < 3:
            sys.exit('ERROR: no %s points within %g km of the impulse'
                     % (field.location, center.hw))
        self.extent = extent
        self.lat, self.lon, self.values = lat[ok], lon[ok], data[ok]


class SectionData:
    """Field along the x or y line through the impulse: values (npts, nlev)."""

    def __init__(self, args, mesh, field, center, axis):
        self.axis = axis
        self.lat, self.lon, self.dist = center.section(axis, args.npts)
        idx = mesh.nearest(field.location, self.lat, self.lon)
        self.values = field.values[idx, :]
        self.cells = mesh.nearest('nCells', self.lat, self.lon) if args.zcoord == 'height' else idx
        self.zz, self.zlabel = vertical_axis(args, mesh, self.cells, field.nlev)
        self.xx = np.tile(self.dist[:, None], (1, field.nlev))


def vertical_axis(args, mesh, cells, nlev):
    """Vertical coordinate for the columns `cells`: (Z array (len(cells), nlev), label)."""
    if args.zcoord == 'height':
        if mesh.zgrid is None:
            print('WARNING: zgrid not found (use --grid-file); falling back to model level')
        else:
            return mesh.layer_heights_km(cells), 'Height (km)'
    return np.tile(np.arange(1, nlev + 1, dtype=float), (len(cells), 1)), 'Model level'


# ---------------------------------------------------------------- panels

def add_map_axes(fig, slot, args, extent):
    if HAVE_CARTOPY:
        ax = fig.add_subplot(slot, projection=ccrs.PlateCarree())
        ax.set_extent(extent, crs=ccrs.PlateCarree())
        if not args.no_coast:
            ax.coastlines(resolution='50m', linewidth=0.6)
            ax.add_feature(cfeature.BORDERS, linewidth=0.4)
            ax.add_feature(cfeature.STATES, linewidth=0.3, edgecolor='gray')
        gl = ax.gridlines(draw_labels=True, linestyle='--', linewidth=0.4)
        gl.top_labels = gl.right_labels = False
        return ax, {'transform': ccrs.PlateCarree()}
    ax = fig.add_subplot(slot)
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    return ax, {}


def draw_map(ax, tr, m, center, colors, gauss):
    clev, cmap, norm = colors
    cs = ax.tricontourf(m.lon, m.lat, m.values, levels=clev, cmap=cmap, norm=norm,
                        extend='both', **tr)
    if gauss:
        gauss.contour(ax, m.lon, m.lat, gauss(m.lat, m.lon, center.level),
                      triangulated=True, **tr)
    for axis in ('x', 'y'):
        slat, slon, _ = center.section(axis, 2)
        ax.plot(slon, slat, 'k:', linewidth=1.0, **tr)
    ax.plot(center.lon, center.lat, 'kx', markersize=8, **tr)
    ax.set_title('Level %d\n%s' % (center.level, stats_text(m.values)), fontsize=10)
    return cs


def draw_section(ax, s, args, mesh, center, colors, gauss):
    clev, cmap, norm = colors
    ax.contourf(s.xx, s.zz, s.values, levels=clev, cmap=cmap, norm=norm, extend='both')
    if gauss:
        levels = np.arange(1, s.values.shape[1] + 1)
        gauss.contour(ax, s.xx, s.zz, gauss(s.lat[:, None], s.lon[:, None], levels[None, :]))
    if args.zcoord == 'height' and mesh.zgrid is not None:
        terrain = mesh.zgrid[s.cells, 0] / 1000.0
        ax.fill_between(s.dist, 0.0, terrain, color='saddlebrown', zorder=3)
    ax.axvline(0.0, color='k', linestyle=':', linewidth=0.8)
    if args.zcoord == 'level':
        ax.axhline(center.level, color='k', linestyle=':', linewidth=0.8)
    if args.zrange:
        ax.set_ylim(args.zrange)
    ax.set_ylabel(s.zlabel)
    if s.axis == 'x':
        ax.set_xlabel('x: distance from impulse (km), west -> east')
        ax.set_title('x section (lat %.2fN)   %s' % (center.lat, stats_text(s.values)), fontsize=10)
    else:
        ax.set_xlabel('y: distance from impulse (km), south -> north')
        ax.set_title('y section (lon %.2fE)   %s' % (center.lon, stats_text(s.values)), fontsize=10)


def plot_variable(args, mesh, field, center, gauss):
    """One figure: map (left), x section (right top), y section (right bottom)."""
    m = MapData(mesh, field, center)
    sx = SectionData(args, mesh, field, center, 'x')
    sy = SectionData(args, mesh, field, center, 'y')
    # One color scale for all three panels, so they can be compared directly.
    colors = color_setup(np.concatenate([m.values, sx.values.ravel(), sy.values.ravel()]), args)

    fig = plt.figure(figsize=(16, 7.5))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.25], hspace=0.45, wspace=0.15)
    ax_map, tr = add_map_axes(fig, gs[:, 0], args, m.extent)
    ax_x = fig.add_subplot(gs[0, 1])
    ax_y = fig.add_subplot(gs[1, 1])

    cs = draw_map(ax_map, tr, m, center, colors, gauss)
    draw_section(ax_x, sx, args, mesh, center, colors, gauss)
    draw_section(ax_y, sy, args, mesh, center, colors, gauss)

    fig.colorbar(cs, ax=[ax_map, ax_x, ax_y], orientation='horizontal',
                 fraction=0.04, pad=0.1, label=field.label())
    gtxt = ('   Gaussian L = %g km x %g km x %g lev' % (gauss.lx, gauss.ly, gauss.lz)
            if gauss else '')
    fig.suptitle('%s, impulse at (%.2fN, %.2fE) level %d%s'
                 % (field.long_name, center.lat, center.lon, center.level, gtxt))
    fname = os.path.join(args.outdir, '%sdirac_%s.%s' % (args.prefix, field.name, args.fmt))
    fig.savefig(fname, dpi=args.dpi, bbox_inches='tight')
    plt.close(fig)
    print('wrote', fname)


# ---------------------------------------------------------------- main

def report_center(mesh, field, center, i_center):
    """Print where the impulse falls and where the response peaks, as a sanity check."""
    lat, lon = mesh.latlon(field.location)
    print('%s: impulse -> %s %d (%.3fN, %.3fE) level %d, value %.4e'
          % (field.name, field.location, i_center, lat[i_center], lon[i_center],
             center.level, field.at(i_center, center.level)))
    ip, kp = np.unravel_index(np.nanargmax(np.abs(field.values)), field.values.shape)
    print('%s: max |value| %.4e at (%.3fN, %.3fE) level %d'
          % (field.name, abs(field.values[ip, kp]), lat[ip], lon[ip], kp + 1))


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    if not HAVE_CARTOPY:
        print('WARNING: cartopy not available; map drawn without projection or coastlines')

    data_nc = Dataset(args.file, 'r')
    grid_nc = Dataset(args.grid_file, 'r') if args.grid_file else None
    mesh = Mesh(data_nc, grid_nc)
    center = Center(*args.center, halfwidth_km=args.halfwidth)
    gauss = GaussOverlay(args, center) if args.gauss else None

    for name in args.vars:
        field = Field(data_nc, name, args.time)
        if not 1 <= center.level <= field.nlev:
            sys.exit('ERROR: --center level %d outside 1..%d for %s' % (center.level, field.nlev, name))
        i_center = int(mesh.nearest(field.location, center.lat, center.lon)[0])
        report_center(mesh, field, center, i_center)
        plot_variable(args, mesh, field, center, gauss)

    data_nc.close()
    if grid_nc is not None:
        grid_nc.close()


if __name__ == '__main__':
    main()
