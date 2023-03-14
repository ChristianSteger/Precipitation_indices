# Description: Check that the fast script 'precipitation_indices.py' yields
#              identical results (apart from numerical errors) compared to a
#              simpler (but less efficient) method in Python that computes
#              indices from unpartitioned time series (e.g. by applying the
#              function 'numpy.percentiles()' to compute percentiles).
#
# Author: Christian R. Steger, March 2023

# Load modules
import os
import numpy as np
import xarray as xr
import textwrap

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------

# Input names
files_pattern = "/Users/csteger/Desktop/data_in/" \
                + "pr_EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_CLMcom-ETH-" \
                + "COSMO-crCLIM-v1-1_v1_1hr_YYYY01010030-YYYY12312330.nc"
file_res = "/Users/csteger/Desktop/data_out/" \
          + "pr_EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_CLMcom-ETH-COSMO-"\
          + "crCLIM-v1-1_v1_1hr_1979-1988_all_day_perc.nc"

# Other settings
years = np.arange(1979, 1988 + 1)  # considered years 1979 - 1988
subdom = {"rlon": slice(170, 230), "rlat": slice(150, 200)}

# -----------------------------------------------------------------------------
# Load data (for sub-domain)
# -----------------------------------------------------------------------------

# Check if all files exist
print("\n".join(textwrap.wrap(("Process files "
                               + files_pattern.split("/")[-1]), 79)))
files = [files_pattern.replace("YYYY", str(i)) for i in years]
if not all([os.path.isfile(i) for i in files]):
    raise ValueError("Not all input files exist")

# Load precipitation data for subdomain
prec = np.empty((0, 50, 60), dtype=np.float32)
prec_yr_max = np.empty((len(files), 50, 60), dtype=np.float32)
for ind, i in enumerate(files):
    ds = xr.open_dataset(i)
    ds = ds.isel(rlon=subdom["rlon"], rlat=subdom["rlat"])
    pr = ds["pr"].values * 3600.0  # unit conversion to [mm hour-1]
    ds.close()
    prec_yr_max[ind, :, :] = np.max(pr, axis=0)
    prec = np.concatenate((prec, pr), axis=0)
    print("File loaded: " + str(ind + 1) + "/" + str(len(files)))

# Load computed precipitation indices
var = ("mean", "max", "wet_day_freq", "intensity",
       "perc_90.00", "perc_95.00", "perc_99.00", "perc_99.90")
ds = xr.open_dataset(file_res)
ds = ds.isel(rlon=subdom["rlon"], rlat=subdom["rlat"])
res = {i: ds[i].values for i in var}
ds.close()

# -----------------------------------------------------------------------------
# Recompute precipitation indices and print maximal (and mean) abs. deviations
# -----------------------------------------------------------------------------

# Temporal mean
prec_mean = prec.mean(axis=0)
dev_max = np.abs(prec_mean - res["mean"]).max()
print("Max. abs. deviation in temporal mean: %.15f" % dev_max + " mm hour-1")

# Hourly maximum
prec_max = prec_yr_max.mean(axis=0)
dev_max = np.abs(prec_max - res["max"]).max()
print("Max. abs. deviation in hourly maximum: %.15f" % dev_max + " mm hour-1")

# Wet day frequency
mask = (prec > 0.1)
wet_day_freq = mask.sum(axis=0).astype(np.float32) / float(mask.shape[0])
dev_max = np.abs(wet_day_freq - res["wet_day_freq"]).max()
print("Max. absolute deviation in wet day frequency: %.15f" % dev_max)

# Intensity
prec_cp = prec.copy()
prec_cp[~mask] = np.nan
intensity = np.nanmean(prec_cp, axis=0)
dev_max = np.abs(intensity - res["intensity"]).max()
print("Max. absolute deviation in intensity: %.10f" % dev_max + " mm hour-1")

# Percentile
qs = np.array([90.0, 95.0, 99.0, 99.9])
prec_per = np.percentile(prec, q=qs, axis=0, method="linear")
print("Absolute deviation in percentiles [mm hour-1]:")
print("       maximum       mean")
for ind, i in enumerate(qs):
    dev_max = np.abs(prec_per[ind, :, :] - res["perc_" + "%.2f" % i]).max()
    dev_mean = np.abs(prec_per[ind, :, :] - res["perc_" + "%.2f" % i]).mean()
    print("%.2f" % i + ": %.10f" % dev_max + ", %.10f" % dev_mean)
