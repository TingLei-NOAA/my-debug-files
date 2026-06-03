import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import math

# List of NetCDF file paths
file_paths = ["Data/mgbf/20240527.000000.jim_old-dbz-lev2-D1_L30-p196_dirac_SABER_lam.fv_core.res.nc", 
 "Data/mgbf/20240527.000000.jim-dbz-lev2-D1_L30-p196_dirac_SABER_lam.fv_core.res.nc" ]

labels = ["dbz-mgbf","jim-dbz-mgbf","xxmgbf_w9_line_normal",
           "xxmgbf_w11_line_normal","mgbf_w13_line_normal","xxxmgbf_w13_line_normal"]

# Specify the horizontal location (adjust these to your desired location)
#X = 143  # Replace with the specific xaxis_1 index you want
Z = 1  # Replace with the specific xaxis_1 index you want
X= 209
Y = 124  # Replace with the specific yaxis_2 index you want
#Y_list=[142,143,146]
#Z_list=[54,52,52]
Lgv=4.885*np.sqrt(2.0)
Lgh=17.8/13.0*np.sqrt(2.0)




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

# Create a figure with two subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# Subplot 1: Vertical profiles
for index, file_path in enumerate(file_paths):
    # Load the NetCDF file
    dataset = nc.Dataset(file_path, mode='r')
    z_1 = dataset.variables['zaxis_1'][:]
    print(f"File {index} z1:", z_1)

    # Extract the variables
    if zaxis_1 is None:
        zaxis_1 = dataset.variables['zaxis_1'][:]
        nz = len(dataset.variables['zaxis_1'])
        zaxis_1 = np.arange(nz)
    if xaxis_1 is None:
#        xaxis_1 = dataset.variables['xaxis_1'][:]
        nx = len(dataset.variables['xaxis_1'])
        xaxis_1 = np.arange(nx)
        zstart=max(math.floor(Z-6*Lgv),0)
        zend=min(math.ceil(Z+6*Lgv),len(zaxis_1)-1)
        xstart=max(math.floor(X-4*Lgh),0)
        xend=min(math.ceil(X+4*Lgh), len(xaxis_1)-1)
        xslice=slice(xstart,xend)
        zslice=slice(zstart,zend)
        print("zslice is ",zslice)
        print("xslice is ",xslice)
    temperature = dataset.variables['air_temperature'][:]  # Assuming 'T' corresponds to temperature

    # Extract the temperature profile at the specified horizontal location
    temperature_Z_profile = temperature[0, :, Y, X]  # Assuming Time = 0
    temperature_X_profile = temperature[0, Z, Y, :]  # Assuming Time = 0
    max_X_temperature = np.max(temperature_X_profile)
    max_Z_temperature = np.max(temperature_Z_profile)
    
    imax_Z = np.argmax(temperature_Z_profile)
    imax_X = np.argmax(temperature_X_profile)
    if index >0 or 1 > 0 :
#          print("do nothing now")
      temperature_X_profile=temperature_X_profile/max_X_temperature
      temperature_Z_profile=temperature_Z_profile/max_Z_temperature


    # Update the label to include the maximum value
    print(f" imax_x is {imax_X}")
    print(f" xaxis_1 {xaxis_1}")
    print(f" max_x is {max_X_temperature}")
    label_with_max_X = f"{labels[index]} "
    label_with_max_Z = f"{labels[index]} (max: {max_Z_temperature:.2f} at Z={zaxis_1[imax_Z]})"
 #   label_with_max_Z = f"{labels[index]} "
    
    # Plot the temperature profile against the vertical axis in the first subplot (ax1)
    ax1.plot(temperature_Z_profile[zslice], zaxis_1[zslice], label=label_with_max_Z)
    ax2.plot(xaxis_1[xslice],temperature_X_profile[xslice],  label=label_with_max_X)

    # Close the dataset
    dataset.close()

vfgaussian=np.exp(-((zaxis_1 - zaxis_1[Z]) ** 2) / (Lgv ** 2))
ax1.plot(vfgaussian[zslice], zaxis_1[zslice], label=f"gaussian with GSI RF  R={Lgv/np.sqrt(2.0):.2f} grid units")
hfgaussian=np.exp(-((xaxis_1 - xaxis_1[X]) ** 2) / (Lgh ** 2))
ax2.plot(xaxis_1[xslice],hfgaussian[xslice],  label=f"gaussian approximating {Lgh*13/np.sqrt(2.0):.2f}km GSI RF length scale")


# Configure the first subplot (Vertical profile plot)
ax1.set_ylabel('Height (Z)')
ax1.set_title(f'Temperature in vertical direction  at xaxis_1={X}, yaxis_1={Y}')
ax1.invert_yaxis()  # Invert the y-axis to have height increasing upwards
ax1.legend()  # Add a legend to distinguish between the different files

ax2.set_xlabel('X')
ax2.set_ylabel('Temperature (T)')
ax2.set_title(f'Temperature in Horizontal Direction at zaxis_1={Z}, yaxis_1={Y}')
ax2.legend()  # Add a legend to distinguish between the different files

# Adjust layout and show the combined plot
plt.tight_layout()
# Show the plot
plt.savefig('dev-jim-dbz-lev2.png')
quit()
