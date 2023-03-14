# Description: Compare output from 'precipitation_indices.py' with results from
#              the shell script applying CDO.
#
# Author: Christian R. Steger, March 2023

# Load modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
from cmcrameri import cm

mpl.style.use("classic")

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------

# Input names
file_py = "/Users/csteger/Desktop/data_out/" \
          + "pr_EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_CLMcom-ETH-COSMO-" \
          + "crCLIM-v1-1_v1_1hr_1979-1988_all_day_perc.nc"
files_cdo = ("/Users/csteger/Desktop/data_in/",
             {"mean": "10ymean_year_mean.nc",
              "max": "10ymean_year_max.nc",
              "wet_day_freq": "10ymean_year_fre.nc",
              "intensity": "10ymean_year_int.nc",
              "perc_90.00": "10ymean_year_q90.nc",
              "perc_95.00": "10ymean_year_q95.nc",
              "perc_99.00": "10ymean_year_q99.nc",
              "perc_99.90": "10ymean_year_q99.9.nc"})
path_plot = "/Users/csteger/Desktop/"

# -----------------------------------------------------------------------------
# Compare results
# -----------------------------------------------------------------------------

# Load precipitation indices from Python script
var = ("mean", "max", "wet_day_freq", "intensity",
       "perc_90.00", "perc_95.00", "perc_99.00", "perc_99.90")
ds = xr.open_dataset(file_py)
res_py = {i: ds[i].values for i in var}
rlon = ds["rlon"].values
rlat = ds["rlat"].values
crs_rot = ccrs.RotatedPole(
    pole_latitude=ds["rotated_pole"].grid_north_pole_latitude,
    pole_longitude=ds["rotated_pole"].grid_north_pole_longitude
)
ds.close()

# Load precipitation indices from CDO shell script
res_cdo = {}
for i in var:
    ds = xr.open_dataset(files_cdo[0] + files_cdo[1][i], decode_times=False)
    if i == "wet_day_freq":
        values = ds["precipitation_days_index_per_time_period"].values
    elif i == "intensity":
        values = ds["simple_daily_intensity_index_per_time_period"].values
    else:
        values = ds["pr"].values
    res_cdo[i] = values.squeeze()
    ds.close()

# Plot comparison
titles = {"mean": "Temporal mean [mm/hour]",
          "max": "Yearly maximum [mm/hour]",
          "wet_day_freq": "Wet day frequency [-]",
          "intensity": "Intensity [mm/hour]",
          "perc_90.00": "90% percentile [mm/hour]",
          "perc_95.00": "95% percentile [mm/hour]",
          "perc_99.00": "99% percentile [mm/hour]",
          "perc_99.90": "99.9% percentile [mm/hour]"}
fig = plt.figure(figsize=(16, 18.5))
gs = gridspec.GridSpec(4 * 3, 4, left=0.1, bottom=0.1, right=0.9, top=0.9,
                       hspace=0.12, wspace=0.1,
                       height_ratios=[1.0, 0.08, 0.08] * 4)
for i in range(2):
    for j in range(4):
        ind = (i * 4) + j
        # ---------------------------------------------------------------------
        ax = plt.subplot(gs[j * 3, (i * 2)], projection=crs_rot)
        data = res_py[var[ind]]
        levels = MaxNLocator(nbins=10, steps=[1, 2, 5, 10], symmetric=False) \
            .tick_values(np.nanpercentile(data, 10.0),
                         np.nanpercentile(data, 90.0))
        cmap = cm.davos_r
        norm = mpl.colors.BoundaryNorm(levels, ncolors=cmap.N, extend="both")
        plt.pcolormesh(rlon, rlat, data, cmap=cmap, norm=norm)
        ax.set_aspect("auto")
        ax.add_feature(cfeature.COASTLINE, linewidth=1.0, color="black")
        ax.add_feature(cfeature.BORDERS, linewidth=1.0, color="black")
        plt.title(titles[var[ind]], fontsize=12, fontweight="bold")
        # ---------------------------------------------------------------------
        ax = plt.subplot(gs[j * 3 + 1, (i * 2)])
        cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm,
                                       orientation="horizontal")
        cb.ax.tick_params(labelsize=10)
        # ---------------------------------------------------------------------
        ax = plt.subplot(gs[j * 3, (i * 2) + 1], projection=crs_rot)
        data = np.abs(res_cdo[var[ind]] - res_py[var[ind]])  # absolute error
        levels = MaxNLocator(nbins=6, steps=[1, 2, 5, 10], symmetric=False) \
            .tick_values(0.0, np.nanpercentile(data, 95.0))
        cmap = cm.lajolla
        norm = mpl.colors.BoundaryNorm(levels, ncolors=cmap.N, extend="max")
        plt.pcolormesh(rlon, rlat, data, cmap=cmap, norm=norm)
        ax.set_aspect("auto")
        ax.add_feature(cfeature.COASTLINE, linewidth=1.0, color="black")
        ax.add_feature(cfeature.BORDERS, linewidth=1.0, color="black")
        plt.title("CDO - Python", fontsize=12, fontweight="bold")
        t = plt.text(0.86, 0.92, "{:.1e}".format(np.nanmax(data)),
                     horizontalalignment="center",
                     verticalalignment="center",
                     transform=ax.transAxes, fontsize=10)
        t.set_bbox(dict(facecolor="white", alpha=1.0, edgecolor="black"))
        # ---------------------------------------------------------------------
        ax = plt.subplot(gs[j * 3 + 1, (i * 2) + 1])
        cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm,
                                       orientation="horizontal", format="%.0e")
        cb.ax.tick_params(labelsize=10)
        # ---------------------------------------------------------------------
fig.savefig(path_plot + "Precipitation_indices_Python_vs_CDO.png", dpi=300,
            bbox_inches="tight")
plt.close(fig)
