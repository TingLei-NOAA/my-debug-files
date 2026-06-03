import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt


X = 209
Y = 124

# Python uses zero-based indices: vertical level 2 -> 1, level 30 -> 29.
LEVEL2_INDEX = 1
LEVEL30_INDEX = 29

VARIABLE = "air_temperature"
OUTPUT_FILE = "dev-jim-dbz-vertical-profiles.png"

EXPERIMENTS = {
    "dbz-mgbf": {
        LEVEL2_INDEX: "Data/mgbf/20240527.000000.jim_old-dbz-lev2-D1_L30-p196_dirac_SABER_lam.fv_core.res.nc",
        LEVEL30_INDEX: "Data/mgbf/20240527.000000.dbz-lev30-D1_L30-p196_dirac_SABER_lam.fv_core.res.nc",
    },
    "jim-dbz-mgbf": {
        LEVEL2_INDEX: "Data/mgbf/20240527.000000.jim-dbz-lev2-D1_L30-p196_dirac_SABER_lam.fv_core.res.nc",
        LEVEL30_INDEX: "Data/mgbf/20240527.000000.jim-dbz-lev30-D1_L30-p196_dirac_SABER_lam.fv_core.res.nc",
    },
}


def read_vertical_profile(file_path):
    with nc.Dataset(file_path, mode="r") as dataset:
        nz = len(dataset.variables["zaxis_1"])
        zaxis = np.arange(nz)
        profile = dataset.variables[VARIABLE][0, :, Y, X]
    return zaxis, profile


def plot_experiment(ax, exp_name, files_by_level):
    zaxis = None
    profiles = {}

    for level_index, file_path in files_by_level.items():
        file_zaxis, profile = read_vertical_profile(file_path)
        if zaxis is None:
            zaxis = file_zaxis
        profiles[level_index] = profile

    norm_coeff = np.max(profiles[LEVEL30_INDEX])
    if norm_coeff == 0.0:
        raise ValueError(f"{exp_name} level-30 profile maximum is zero; cannot normalize")

    for level_index in (LEVEL2_INDEX, LEVEL30_INDEX):
        profile = profiles[level_index]
        max_value = np.max(profile)
        imax = np.argmax(profile)
        label = (
            f"Dirac level {level_index + 1} "
            f"(max={max_value:.3g} at z={zaxis[imax]})"
        )
        ax.plot(profile / norm_coeff, zaxis, label=label)

    ax.set_title(exp_name)
    ax.set_xlabel(f"{VARIABLE} / max(level 30 profile)")
    ax.set_ylabel("Vertical index")
    ax.invert_yaxis()
    ax.grid(True, alpha=0.3)
    ax.legend()


fig, axes = plt.subplots(1, 2, figsize=(13, 6), sharey=True)

for ax, (exp_name, files_by_level) in zip(axes, EXPERIMENTS.items()):
    plot_experiment(ax, exp_name, files_by_level)

fig.suptitle(f"Normalized vertical profiles at xaxis_1={X}, yaxis_1={Y}")
plt.tight_layout()
plt.savefig(OUTPUT_FILE, dpi=150)
