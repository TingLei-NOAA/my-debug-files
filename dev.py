import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import math
from pathlib import Path

# List of NetCDF file paths
file_paths = [
            "Data/mgbf_NA/20240527.010000.p484_L30_proc5_2G_35kmv6-mgbf-D1-p64_dirac_SABER_lam.fv_core.res.nc",
            "Data/mgbf_NA/20240527.010000.p484-inho_L30_proc5_2G_35kmv6-mgbf-D1-p64_dirac_SABER_lam.fv_core.res.nc"]

labels = ["old mgbf","new inho mgbf","new-mgbf","new_mgbf_w13_line",
           "xxmgbf_w11_line_normal","mgbf_w13_line_normal","xxxmgbf_w13_line_normal"]

# Specify the horizontal location (adjust these to your desired location)
#X = 143  # Replace with the specific xaxis_1 index you want
Z = 29  # Replace with the specific xaxis_1 index you want
X= 149 
Y = 149  # Replace with the specific yaxis_2 index you want
#Y_list=[142,143,146]
#Z_list=[54,52,52]
Lgv=4.0
Lgh=7.5*13.0/3.0




# Assuming `file_paths` and `labels` are defined somewhere before this code
# Example:
# file_paths = ["file1.nc", "file2.nc"]
# labels = ["File 1", "File 2"]

# Variables for horizontal plotting
Z_list = [0, 1]  # Modify these according to your needs
Y_list = [5, 10]  # Modify these according to your needs

# Set the horizontal and vertical coordinates (replace these with actual values)
#X = 10
#Y = 20

# Placeholder variables to store `zaxis_1`
zaxis_1 = None
xaxis_1 = None
yaxis_1 = None


def build_xy_from_dxdy(dx: np.ndarray, dy: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Build 2D physical coordinates (meters) from dx, dy on the T grid.
    x increases eastward using cumulative dx along the second dimension.
    y increases northward using cumulative dy along the first dimension.
    Origin is (0,0) at [0,0].
    """
    ny, nx = dx.shape
    x = np.zeros((ny, nx), dtype=float)
    y = np.zeros((ny, nx), dtype=float)
    # accumulate x along columns
    for j in range(ny):
        for i in range(1, nx):
            x[j, i] = x[j, i - 1] + dx[j, i - 1]
    # accumulate y along rows
    for j in range(1, ny):
        y[j, :] = y[j - 1, :] + dy[j - 1, :]
    return x, y

# Create a figure with three subplots
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))

# Load dx/dy grid for physical horizontal coordinates and build 2D x/y (meters)
dxdy_path = Path("fv3_grid_dxdy.nc")
if not dxdy_path.exists():
    raise FileNotFoundError(f"{dxdy_path} not found; place fv3_grid_dxdy.nc alongside dev.py")
dxdy_nc = nc.Dataset(dxdy_path, mode="r")
dx_grid = dxdy_nc.variables["dx"][:]  # (grid_yt, grid_xt) meters
dy_grid = dxdy_nc.variables["dy"][:]  # (grid_yt, grid_xt) meters
x_grid, y_grid = build_xy_from_dxdy(dx_grid, dy_grid)

# Subplot 1: Vertical profiles
for index, file_path in enumerate(file_paths):
    # Load the NetCDF file
#  if index == 0 :
  if True :
    dataset = nc.Dataset(file_path, mode='r')
    z_1 = dataset.variables['zaxis_1'][:]
    z_1[:] = np.arange(z_1.size).reshape(z_1.shape)
    print(f"File {index} z1:", z_1)

    # Extract the variables
    if zaxis_1 is None:
        zaxis_1 = dataset.variables['zaxis_1'][:]
        zaxis_1[:] = np.arange(zaxis_1.size).reshape(zaxis_1.shape)
        zstart=max(math.floor(Z-6*Lgv),0)
        zend=min(math.ceil(Z+6*Lgv),len(zaxis_1)-1)
        zslice=slice(zstart,zend)
    if xaxis_1 is None:
        xstart=max(math.floor(X-8*Lgh),0)
        xend=min(math.ceil(X+8*Lgh), x_grid.shape[1]-1)
        xslice=slice(xstart,xend)
        xaxis_1 = x_grid[Y, :]
        print("zslice is ",zslice)
        print("xslice is ",xslice)
    if yaxis_1 is None:
        ystart=max(math.floor(Y-8*Lgh),0)
        yend=min(math.ceil(Y+8*Lgh), y_grid.shape[0]-1)
        yslice=slice(ystart,yend)
        yaxis_1 = y_grid[:, X]
    temperature = dataset.variables['air_temperature'][:]  # Assuming 'T' corresponds to temperature
    print(f"thinkdeb max temp is {np.max(temperature)}")

    # Extract the temperature profile at the specified horizontal location
    temperature_Z_profile = temperature[0, :, Y, X]  # Assuming Time = 0
    temperature_X_profile = temperature[0, Z, Y, :]  # Assuming Time = 0
    temperature_Y_profile = temperature[0, Z, :, Y]  # Assuming Time = 0
    max_X_temperature = np.max(temperature_X_profile)
    max_Y_temperature = np.max(temperature_Y_profile)
    max_Z_temperature = np.max(temperature_Z_profile)
    
    imax_Z = np.argmax(temperature_Z_profile)
    imax_X = np.argmax(temperature_X_profile)
    imax_Y = np.argmax(temperature_Y_profile)
#    if index >0 and 1 > 2 :
    if True :
      print("do nothing now")
      temperature_X_profile=temperature_X_profile/max_X_temperature
      temperature_Y_profile=temperature_Y_profile/max_Y_temperature
      temperature_Z_profile=temperature_Z_profile/max_Z_temperature


    # Update the label to include the maximum value
#    print(f" imax_x is {imax_X}")
#    print(f" xaxis_1 {xaxis_1}")
#    print(f" yaxis_1 {xaxis_1}")
#    if index == 0:
#       label_with_max_X = f"{labels[index]} 250km cutoff radius "
#       label_with_max_Y = f"{labels[index]} 250km cutoff radius "
    else: 
       label_with_max_X = f"{labels[index]} 2D line filter "
    label_with_max_Z = f"{labels[index]} (max: {max_Z_temperature:.2f} at Z={zaxis_1[imax_Z]})"
    label_with_max_X = f"{labels[index]} (max: {max_X_temperature:.2f} at X={xaxis_1[imax_X]})"
    label_with_max_Y = f"{labels[index]} (max: {max_Y_temperature:.2f} at y={yaxis_1[imax_Y]})"
#    if index == 0:
#      label_with_max_Z = f"{labels[index]} 0.3 cut-off radius in sigma"
#    else:
#      label_with_max_Z = f"{labels[index]} 2D line filter"
   
#clt    label_with_max_Z = f"{labels[index]} "
    
    # Plot the temperature profile against the vertical axis in the first subplot (ax1)
    ax1.plot(temperature_Z_profile[zslice], zaxis_1[zslice], label=label_with_max_Z)
    ax2.plot(xaxis_1[xslice],temperature_X_profile[xslice],  label=label_with_max_X)
    ax2.plot(yaxis_1[yslice],temperature_Y_profile[yslice],  label=label_with_max_Y)

    # Horizontal contour around (X,Y)
    temp_slice = temperature[0, Z, yslice, xslice]
    Xmesh, Ymesh = np.meshgrid(xaxis_1[xslice], yaxis_1[yslice])
    cf = ax3.pcolormesh(Xmesh, Ymesh, temp_slice, shading="auto")
    ax3.plot(xaxis_1[X], yaxis_1[Y], marker="o", color="red", markersize=6, label="(X,Y)")
    ax3.set_xlabel("X (m)")
    ax3.set_ylabel("Y (m)")
    ax3.set_title(f"Temperature at Z={Z}")
    cb = plt.colorbar(cf, ax=ax3, orientation="vertical", pad=0.02)
    cb.set_label("Temperature (units)")
    ax3.legend()

    # Close the dataset
    dataset.close()

vfgaussian=np.exp(-((zaxis_1 - zaxis_1[Z]) ** 2) / (Lgv ** 2))
#ax1.plot(vfgaussian[zslice], zaxis_1[zslice], label=f"Gaussian Curve with R={Lgv} grid units")
# Horizontal axes now meters; scale Lgh accordingly if desired (kept as-is)
hfgaussian=np.exp(-((xaxis_1 - xaxis_1[X]) ** 2) / (Lgh ** 2))
hyfgaussian=np.exp(-((yaxis_1 - yaxis_1[X]) ** 2) / (Lgh ** 2))
#ax2.plot(xaxis_1[xslice],hfgaussian[xslice],  label=f"Gaussian Curve approximating 250km cut-off")
#ax3.plot(yaxis_1[yslice],hyfgaussian[yslice],  label=f"Gaussian Curve approximating 250km cut-off")


# Configure the first subplot (Vertical profile plot)
#ax3.set_xlabel('Response Amplitude')
ax2.set_xlabel('Response Amplitude')
ax1.set_ylabel('Height(Grid Units')
ax1.set_title(f'Response Amplitude in the vertical Direction')
ax1.invert_yaxis()  # Invert the y-axis to have height increasing upwards
ax1.legend()  # Add a legend to distinguish between the different files

ax2.set_xlabel('Horizontal distance (m)')
ax2.set_ylabel('Response Amplitude')
ax2.set_title(f'Response Amplitude in the X-Direction')
ax2.legend()  # Add a legend to distinguish between the different files

# Adjust layout and show the combined plot
plt.tight_layout()
# Show the plot
plt.savefig('deb.png')
dxdy_nc.close()
