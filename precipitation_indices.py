# Description: Compute multiple precipitation indices (temporal mean, daily/
#              hourly maximum, wet day frequency and percentiles) from NetCDF
#              climate data with high spatiotemporal resolution. Statistics are
#              computed from yearly blocks and subsequently merged. For per-
#              centiles, this is achieved by keeping the largest occurring
#              precipitation values per grid cell (and updating them during the
#              iteration through the yearly blocks). All specified percentiles
#              can then be computed from the maximal values kept in memory.
#
# Literature referenced:
# - Ban et al. (2021): The first multi‐model ensemble of regional climate
#   simulations at kilometer‐scale resolution, part I: evaluation of
#   precipitation, https://doi.org/10.1007/s00382-021-05708-w
# - Schär et al. (2016): Percentile indices for assessing changes in heavy
#   precipitation events, https://doi.org/10.1007/s10584-016-1669-2
#
# Author: Christian R. Steger, March 2023

# Load modules
import sys
import os
import numpy as np
import time
import xarray as xr
import textwrap
import warnings

# Load required functions
sys.path.append("/Users/csteger/Downloads/Precipitation_indices/")
# set to path in which file 'auxiliary.py' is stored
import auxiliary as aux

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------

# Input/Output
files_pattern = "/Users/csteger/Desktop/data_in/pr_EUR-11_ECMWF-ERAINT_"\
                + "evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_1hr_" \
                + "YYYY01010030-YYYY12312330.nc"
# input file pattern -> important: label year(s) with 'YYYY'
path_out = "/Users/csteger/Desktop/data_out/"  # output directory
file_out_fp = "pr_EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_CLMcom-ETH-" \
              + "COSMO-crCLIM-v1-1_v1_1hr"
# first part of output name

# Other settings
time_freq_in = "1hr"  # temporal input frequency ("1hr" or "day")
years = np.arange(1979, 1988 + 1)  # considered years 1979 - 1988
percentile_method = "all"
# percentile method according to Schär et al. (2016) ("all", "wet")
qs = np.array([90.0, 95.0, 99.0, 99.9])  # percentiles [0.0, 100.0]
year_subsel = "yearly"
# temporal subselection for year ("yearly", "JJA", "SON", "DJF", "MAM")
prec_thresh = {"1hr": 0.1, "day": 1.0}  # [mm]
# threshold for wet day/hour according to Ban et al. (2021)
unit_con_fac = {"1hr": 3600.0, "day": 1.0}
# factor to convert input units to precipitation flux per temporal input
# frequency. With the current setting, hourly data must have input units
# [kg m-2 s-1] and daily data [mm day-1]
block_size_max_in = 2.5  # None, 2.5, 5.0
# maximal data volume that is read during an opening call (None or value) [GB]

# Important notes:
# - Opening larger NetCDF files (ca. > 5 GB) can cause the xarray/NetCDF
#   function to freeze. To avoid this, a maximal value can be specified with
#   'block_size_max_in'. With 'block_size_max_in = None', all data is loaded in
#   one step.

# -----------------------------------------------------------------------------
# Preprocessing steps
# -----------------------------------------------------------------------------
print(" Preprocessing steps " .center(79, "#"))

# Check input settings
if time_freq_in not in ("1hr", "day"):
    raise ValueError("Unknown selected temporal granularity")
if percentile_method not in ("all", "wet"):
    raise ValueError("Unknown value for 'percentile_method'")
if (qs.min() < 85.0) or (qs.max() > 100.0):
    raise ValueError("Allowed range for qs of [85.0, 100.0] is exceeded")
if year_subsel not in ("yearly", "JJA", "SON", "DJF", "MAM"):
    raise ValueError("Unknown value for 'year_subsel'")

# Check if all files exist
print("\n".join(textwrap.wrap(("Process files "
                               + files_pattern.split("/")[-1]), 79)))
files = [files_pattern.replace("YYYY", str(i)) for i in years]
if not all([os.path.isfile(i) for i in files]):
    raise ValueError("Not all input files exist")

# Load ad check metadata from NetCDF
ds = xr.open_dataset(files[0], decode_times=False)
if "calendar" in list(ds["time"].attrs):
    mod_cal = ds["time"].calendar
elif "calender" in list(ds["time"].attrs):
    mod_cal = ds["time"].calender
else:
    raise ValueError("Calendar attribute not found")
if mod_cal not in ("standard", "gregorian", "proleptic_gregorian", "360_day"):
    raise ValueError("Unknown calendar")
if "rlon" in list(ds.coords):
    len_x = ds.coords["rlon"].size
    len_y = ds.coords["rlat"].size
    out_dim = ("rlat", "rlon")
elif "x" in list(ds.coords):
    len_x = ds.coords["x"].size
    len_y = ds.coords["y"].size
    out_dim = ("y", "x")
else:
    raise ValueError("Unknown spatial coordinates")
ds.close()

# Compute total number of time steps
map_freq = {"1hr": "H", "day": "D"}
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    time_axis = xr.date_range(start=str(years[0]), end=str(years[-1] + 1),
                              freq=map_freq[time_freq_in], calendar=mod_cal,
                              closed="left")
da = xr.DataArray(time_axis, [("time", time_axis)])
if year_subsel != "yearly":
    da = da.sel(time=da["time.season"] == year_subsel)
num_ts_tot = da["time"].size
print("Total number of time steps: " + str(num_ts_tot))

# -----------------------------------------------------------------------------
# Compute temporal mean, daily/hourly maximum, wet day frequency and
# percentiles for year or season
# -----------------------------------------------------------------------------
print(" Compute precipitation indices " .center(79, "#"))

t_beg_tot = time.time()

# Allocate arrays
prec_mean_yr = np.empty((len(years), len_y, len_x), dtype=np.float32)
prec_max_yr = np.empty((len(years), len_y, len_x), dtype=np.float32)
ts_above_thresh_yr = np.empty((len(years), len_y, len_x), dtype=np.int32)
prec_int_yr = np.empty((len(years), len_y, len_x), dtype=np.float32)
num_keep = (np.ceil(num_ts_tot * (1.0 - qs.min() / 100)) + 1) \
    .astype(np.int32)
prec_keep = np.empty((len_y, len_x, num_keep), dtype=np.float32)
prec_keep.fill(-999.0)
print("Size of array 'prec_keep': %.1f" % (prec_keep.nbytes / (10 ** 9))
      + " GB")

# Conversion of input precipitation unit to [mm h-1] or [mm day-1]
out_unit = {"1hr": "mm h-1", "day": "mm day-1"}

# Loop through years
num_ts_yr = np.empty(len(years), dtype=np.int32)
for ind, year in enumerate(years):

    print((" Process year " + str(year) + " ").center(79, "-"))

    # Load data
    t_beg = time.time()
    file_in = files[ind]
    print(textwrap.fill("Load file " + file_in.split("/")[-1], 79))
    # -------------------------------------------------------------------------
    # Load data in one block
    # -------------------------------------------------------------------------
    if block_size_max_in is None:
        ds = xr.open_dataset(file_in)
        if year_subsel != "yearly":
            ds = ds.sel(time=da["time.season"] == year_subsel)
        len_t = ds.coords["time"].size
        prec = ds["pr"].values
        ds.close()
    # -------------------------------------------------------------------------
    # Load data in multiple blocks
    # -------------------------------------------------------------------------
    else:
        ds = xr.open_dataset(file_in)
        if year_subsel != "yearly":
            ds = ds.sel(time=da["time.season"] == year_subsel)
        len_t = ds.coords["time"].size
        ds.close()
        prec = np.empty((len_t, len_y, len_x), dtype=np.float32)
        prec.fill(np.nan)
        num_blocks = int(np.ceil((prec.nbytes / (10 ** 9))
                                 / block_size_max_in))
        lim = np.linspace(0, len_t, num_blocks + 1, dtype=np.int32)
        for i in range(num_blocks):
            ds = xr.open_dataset(file_in)
            if year_subsel != "yearly":
                ds = ds.sel(time=da["time.season"] == year_subsel)
            slice_t = slice(lim[i], lim[i + 1])
            prec[slice_t, :, :] = ds["pr"][slice_t, :, :].values
            ds.close()
            print("Data blocks loaded: " + str(i + 1) + "/" + str(num_blocks))
    # -------------------------------------------------------------------------
    prec *= unit_con_fac[time_freq_in]  # conversion to [mm h-1] or [mm day-1]
    print("Data loaded (" + "%.1f" % (time.time() - t_beg) + " s)")
    num_ts_yr[ind] = len_t

    # Check range of precipitation data
    prec_min = prec.min()
    if prec_min < 0.0:
        print("Warning: negative precipitation value(s) detected")
        print("Minimum value: %.5f" % prec_min + " " + out_unit[time_freq_in])
        if prec_min < -1.0:
            raise ValueError("Negative precipitation smaller than -1.0 "
                             + out_unit[time_freq_in] + " found")
    if np.isnan(prec_min):
        raise ValueError("NaN-value(s) detected")

    # Compute temporal mean
    t_beg = time.time()
    prec_mean_yr[ind, :, :] = prec.mean(axis=0)
    print("Mean computed (" + "%.1f" % (time.time() - t_beg) + " s)")

    # Compute daily/hourly maximum
    t_beg = time.time()
    prec_max_yr[ind, :, :] = prec.max(axis=0)
    print("Daily/hourly maximum (" + "%.1f" % (time.time() - t_beg) + " s)")

    # Compute wet day frequency
    t_beg = time.time()
    mask_inc = (prec > prec_thresh[time_freq_in])
    ts_above_thresh_yr[ind, :, :] = mask_inc.sum(axis=0)
    print("Wet day frequency computed (" + "%.1f" % (time.time() - t_beg)
          + " s)")

    # Update maximal values to compute percentiles later
    t_beg = time.time()
    if percentile_method == "all":
        aux.update_max_values_all_day(prec_keep, prec, len_y, len_x, num_keep)
    else:
        aux.update_max_values_wet_day(prec_keep, prec, len_y, len_x, num_keep,
                                      prec_thresh[time_freq_in])
    print("Maximal values updated (" + "%.1f" % (time.time() - t_beg) + " s)")

    # Intensity
    t_beg = time.time()
    prec[~mask_inc] = 0.0
    num_ts_cons = ts_above_thresh_yr[ind, :, :].astype(np.float32)
    mask_issue = (num_ts_cons == 0.0)
    if np.any(mask_issue):
        print("Warning: some grid cells never exceed the threshold value")
        print("Set affected grid cells to NaN")
        num_ts_cons[mask_issue] = np.nan
        # avoid division by 0.0 by setting the grid cells to NaN
    prec_int_yr[ind, :, :] = np.mean(prec, axis=0) \
        * (prec.shape[0] / num_ts_cons)
    print("Intensity computed (" + "%.1f" % (time.time() - t_beg) + " s)")

print((" Aggregate statistics over " + str(len(years)) + " years ")
      .center(79, "-"))

# Compute temporal mean, daily/hourly maximum and wet day frequency for entire
# period
prec_mean = np.average(prec_mean_yr, weights=num_ts_yr, axis=0) \
    .astype(np.float32)
prec_max = prec_max_yr.mean(axis=0)
if num_ts_tot != num_ts_yr.sum():
    raise ValueError("Inconsistency in total number of time steps")
ts_above_thresh = ts_above_thresh_yr.sum(axis=0).astype(np.int32)
prec_wet_day_freq = (ts_above_thresh / num_ts_tot).astype(np.float32)

# Compute intensity for entire period
prec_int = aux.weighted_average(prec_int_yr, weights=ts_above_thresh_yr,
                                weight_zero="ignore")
# -> with 'weight_zero="ignore"', the intensity of a grid cell is only set to
#    NaN when all years have invalid values. With 'weight_zero="consider"',
#    a single invalid value in a year causes the grid cell to be NaN.
num_gc = np.any(ts_above_thresh_yr == 0, axis=0).sum()
num_gc_all = np.all(ts_above_thresh_yr == 0, axis=0).sum()
print("Number of grid cells with no values above threshold for some / all")
print("years: " + str(num_gc) + " / " + str(num_gc_all))
print("Number of grid cells with intensity of NaN: "
      + str(np.isnan(prec_int).sum()))

# Compute percentiles for entire period
prec_per = np.empty((len(qs), len_y, len_x), dtype=np.float32)
prec_per.fill(np.nan)
if percentile_method == "all":
    print("Compute all day precipitation percentiles")
    x = np.linspace(0.0, 100.0, num_ts_tot, dtype=np.float32)
    if qs.min() < x[-num_keep]:
        raise ValueError("x-position for interpolation is out of range")
    for i in range(len_y):
        for j in range(len_x):
            prec_per[:, i, j] = np.interp(qs, x[-num_keep:],
                                          prec_keep[i, j, :])
else:
    print("Compute wet day precipitation percentiles")
    for i in range(len_y):
        for j in range(len_x):
            if ts_above_thresh[i, j] == 0:
                continue
            x = np.linspace(0.0, 100.0, ts_above_thresh[i, j],
                            dtype=np.float32)
            if ts_above_thresh[i, j] >= num_keep:
                if qs.min() < x[-num_keep]:
                    raise ValueError("x-position for interpolation is out of "
                                     + "range")
                prec_per[:, i, j] = np.interp(qs, x[-num_keep:],
                                              prec_keep[i, j, :])
            else:
                prec_per[:, i, j] = np.interp(qs, x, prec_keep[i, j, -len(x):])

# -----------------------------------------------------------------------------
# Save precipitation indices to NetCDF file
# -----------------------------------------------------------------------------
print(" Save statistics to NetCDF file ".center(79, "#"))

# Processing information and addition to output file name
if year_subsel == "yearly":
    info = "Considered years: " + str(years[0]) + " - " + str(years[-1]) \
           + ", threshold for wet day: %.2f" % prec_thresh[time_freq_in] \
           + " mm"
    fn_add = str(years[0]) + "-" + str(years[-1]) + "_" + percentile_method \
        + "_day_perc"
else:
    info = "Considered years: " + str(years[0]) + " - " + str(years[-1]) \
           + ", sub-yearly period: " + str(year_subsel) + ", threshold for " \
           + "wet day: %.2f" % prec_thresh[time_freq_in] + " mm"
    fn_add = str(years[0]) + "-" + str(years[-1]) + "_" + str(year_subsel) \
        + "_" + percentile_method + "_day_perc"

# Save to NetCDF file
nan_val = -999.0
ds = xr.open_mfdataset(files[0])
ds = ds.drop_dims("time")
ds.attrs["precipitation_indices"] = info
ds["mean"] = (out_dim, prec_mean)
ds["mean"].attrs["units"] = out_unit[time_freq_in]
ds["mean"].attrs["long_name"] = year_subsel + " mean averaged over all years"
ds["max"] = (out_dim, prec_max)
ds["max"].attrs["units"] = out_unit[time_freq_in]
ds["max"].attrs["long_name"] = year_subsel + " maximum averaged over all years"
ds["wet_day_freq"] = (out_dim, prec_wet_day_freq)
ds["wet_day_freq"].attrs["units"] = "-"
ds["wet_day_freq"].attrs["long_name"] = "wet day frequency"
ds["intensity"] = (out_dim, np.nan_to_num(prec_int, nan=nan_val))
ds["intensity"].attrs["units"] = out_unit[time_freq_in]
ds["intensity"].attrs["long_name"] = "intensity"
for ind, q in enumerate(qs):
    name = "perc_%.2f" % q
    ds[name] = (out_dim, np.nan_to_num(prec_per[ind, :, :], nan=nan_val))
    ds[name].attrs["units"] = out_unit[time_freq_in]
    ds[name].attrs["long_name"] = "%.2f" % q + " " + percentile_method \
                                  + "-day percentile"
encoding_nan = {i: {"_FillValue": nan_val} for i in ["intensity"]
                + ["perc_%.2f" % i for i in qs]}
encoding_no_nan = {i: {"_FillValue": None} for i in
                   ["mean", "max", "wet_day_freq"]
                   + list({"rlon", "rlat", "x", "y", "lon", "lat"}
                          .intersection(set(ds.variables)))}
ds.to_netcdf(path_out + file_out_fp + "_" + fn_add + ".nc",
             encoding=(encoding_nan | encoding_no_nan))

print("Total elapsed time: %.1f" % (time.time() - t_beg_tot) + " s")
