import math
from typing import Tuple

import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
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
Lgh=110




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

# Create a figure with two subplots
fig, (ax1, ax2,ax3) = plt.subplots(1, 3, figsize=(15, 6))

# Subplot 1: Vertical profiles
for index, file_path in enumerate(file_paths):
    # Load the NetCDF file
#  if index == 0 :
  if True :
    dataset = nc.Dataset(file_path, mode='r')
    z_1 = dataset.variables['zaxis_1'][:]
    print(f"File {index} z1:", z_1)

    # Extract the variables
    if zaxis_1 is None:
        zaxis_1 = dataset.variables['zaxis_1'][:]
        zstart=max(math.floor(Z-6*Lgv),0)
        zend=min(math.ceil(Z+6*Lgv),len(zaxis_1)-1)
        zslice=slice(zstart,zend)
    if xaxis_1 is None:
        xaxis_1 = dataset.variables['xaxis_1'][:]
        x0 = xaxis_1[X]
        xmask = np.where(np.abs(xaxis_1 - x0) <= 8 * Lgh)[0]
        xslice = slice(xmask[0], xmask[-1] + 1) if xmask.size else slice(0, len(xaxis_1))
        print("zslice is ",zslice)
        print("xslice is ",xslice)
    if yaxis_1 is None:
        yaxis_1 = dataset.variables['yaxis_1'][:]
        y0 = yaxis_1[Y]
        ymask = np.where(np.abs(yaxis_1 - y0) <= 8 * Lgh)[0]
        yslice = slice(ymask[0], ymask[-1] + 1) if ymask.size else slice(0, len(yaxis_1))
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
    ax3.plot(yaxis_1[yslice],temperature_Y_profile[yslice],  label=label_with_max_Y)
#    ax2.plot(xaxis_1,temperature_X_profile,  label=label_with_max_X)

    # Close the dataset
    dataset.close()

vfgaussian=np.exp(-((zaxis_1 - zaxis_1[Z]) ** 2) / (Lgv ** 2))
ax1.plot(vfgaussian[zslice], zaxis_1[zslice], label=f"Gaussian Curve with R={Lgv} grid units")
hfgaussian=np.exp(-((xaxis_1 - xaxis_1[X]) ** 2) / (Lgh ** 2))
hyfgaussian=np.exp(-((yaxis_1 - yaxis_1[Y]) ** 2) / (Lgh ** 2))
ax2.plot(xaxis_1[xslice],hfgaussian[xslice],  label=f"Gaussian Curve approximating RF 110km ")
ax3.plot(yaxis_1[yslice],hyfgaussian[yslice],  label=f"Gaussian Curve approximating RF 110km")


# Configure the first subplot (Vertical profile plot)
ax3.set_xlabel('Response Amplitude')
ax2.set_xlabel('Response Amplitude')
ax1.set_ylabel('Height(Grid Units')
ax1.set_title(f'Response Amplitude in the vertical Direction')
ax1.invert_yaxis()  # Invert the y-axis to have height increasing upwards
ax1.legend()  # Add a legend to distinguish between the different files

ax3.set_xlabel('Y (Grid Units)')
ax3.set_ylabel('Response Amplitude')
ax3.set_title(f'Response Amplitude in the Y-Direction')
ax3.legend()  # Add a legend to distinguish between the different files
ax2.set_xlabel('X (Grid Units)')
ax2.set_ylabel('Response Amplitude')
ax2.set_title(f'Response Amplitude in the X-Direction')
ax2.legend()  # Add a legend to distinguish between the different files

# Adjust layout and show the combined plot
plt.tight_layout()
# Show the plot
plt.savefig('deb-xy.png')
quit()
