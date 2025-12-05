import os
import glob
import re
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

def read_latlon_files(pattern):
    latitudes = []
    longitudes = []

    for filename in glob.glob(pattern):
        with open(filename, 'r') as file:
            content = file.read()
            # Extract the array part from the content using regex
            match = re.search(r'array=\[.*?values:\s*\[(.*?)\]', content, re.DOTALL)
            if match:
                array_str = match.group(1).strip()
                # Convert to a list of floats
                values = [float(x) for x in re.findall(r'-?\d+\.?\d*', array_str)]
                # Separate the values into latitudes and longitudes
                lons = values[0::2]
                lats = values[1::2]
                if len(lons) == len(lats):
                    longitudes.extend(lons)
                    latitudes.extend(lats)
                else:
                    print(f"Warning: Mismatched lat/lon pairs in {filename}")
            else:
                print(f"Warning: No valid array found in {filename}")

    return np.array(longitudes), np.array(latitudes)

def is_float(element):
    try:
        float(element)
        return True
    except ValueError:
        return False

def normalize_longitudes(longitudes):
    """Normalize longitudes to be within the range [-180, 180]."""
    return ((longitudes + 180) % 360) - 180

def plot_latlon(longitudes, latitudes, longitudes2=None, latitudes2=None, output_file=None, padding=5,
    center_lon=-125,point_size=2):
    # Normalize longitudes
    longitudes = normalize_longitudes(longitudes)
  
    if longitudes2 is not None:
        longitudes2 = normalize_longitudes(longitudes2)
    proj = ccrs.PlateCarree(central_longitude=center_lon)              # CHANGE: centered map
#    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
    fig, ax = plt.subplots(subplot_kw={'projection': proj})
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')

    # Calculate the min and max for longitudes and latitudes with padding in degrees
    if longitudes2 is not None and latitudes2 is not None:
        lonmin = min(longitudes.min(), longitudes2.min()) - padding
        lonmax = max(longitudes.max(), longitudes2.max()) + padding
        latmin = min(latitudes.min(), latitudes2.min()) - padding
        latmax = max(latitudes.max(), latitudes2.max()) + padding
    else:
        lonmin = longitudes.min() - padding
        lonmax = longitudes.max() + padding
        latmin = latitudes.min() - padding
        latmax = latitudes.max() + padding

    # Ensure that the ranges form a valid bounding box
#    if lonmin == lonmax:
#        lonmin -= padding
#        lonmax += padding
#    if latmin == latmax:
#        latmin -= padding
#        latmax += padding

    # Ensure the limits are reasonable
#    if lonmin < -180: lonmin = -180
#    if lonmax > 180: lonmax = 180
#    if latmin < -90: latmin = -90
#    if latmax > 90: latmax = 90

    print(f'lonmin/max: {lonmin} / {lonmax}')
    print(f'latmin/max: {latmin} / {latmax}')

#    ax.set_xlim([lonmin, lonmax])
#    ax.set_ylim([latmin, latmax])
 # --- set extent with explicit CRS & avoid set_xlim/ylim ---
    data_crs = ccrs.PlateCarree()                                       # CHANGE: explicit data CRS
    ax.set_extent([lonmin, lonmax, latmin, latmax], crs=data_crs)       # CHANGE: set_extent over xlim/ylim


    # Increase the size of the scattered points
    print("Plotting first set of points:")
    print("Longitudes:", longitudes)
    print("Latitudes:", latitudes)
    ax.scatter(longitudes, latitudes, alpha=0.6,color='blue', s=point_size, label='mgbf filtering grid',
               transform=data_crs)

    if longitudes2 is not None and latitudes2 is not None:
        print("Plotting second set of points:")
        print("Longitudes2:", longitudes2)
        print("Latitudes2:", latitudes2)
        ax.scatter(longitudes2, latitudes2, color='red', alpha=0.6,s=point_size, label='FV3 grid',
        transform=data_crs)
        ax.legend()

    # Add grid lines with labels
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
    gl.top_labels = False
    gl.right_labels = False

    # Optional: format the labels
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}

    if output_file:
        plt.savefig(output_file)
    else:
        plt.show()






def main():
    # Pattern for the first set of files
    pattern1 = 'mgbf_filtering*latlon_*.txt'
    longitudes1, latitudes1 = read_latlon_files(pattern1)

    # Uncomment and set the pattern for the second set of files if needed
    pattern2 = 'model_native_grid_latlon_*.txt'
    longitudes2, latitudes2 = read_latlon_files(pattern2)
    plot_latlon(longitudes1, latitudes1, longitudes2, latitudes2, output_file='halov1-newinter-domains-plot.png')

    # For plotting a single set of points
#    plot_latlon(longitudes1, latitudes1, output_file='latlon_plot.png')

if __name__ == '__main__':
    main()


