# Description: Compute multiple precipitation indices with CDO
#
# Author: Christian R. Steger, March 2023

# -----------------------------------------------------------------------------

# Concatenate files and convert to [mm h-1] (~70s)
files_pattern="pr_EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_CLMcom-ETH-COSMO-"\
"crCLIM-v1-1_v1_1hr_????01010030-????12312330.nc"
cdo -L -mulc,3600 -cat ${files_pattern} 10y_tmp.nc

# Temporal mean (~40s)
cdo -s -L -timmean -yearmean 10y_tmp.nc 10ymean_YEAR_mean.nc

# Wet day frequency (R>0.1mm/day) (~50s)
cdo -s eca_pd,0.1 10y_tmp.nc 10yfre.nc
cdo -s divc,87672 10yfre.nc 10ymean_YEAR_fre.nc

# Intensity (mean precipitation amount on wet days) (~50s)
cdo -s eca_sdii,0.1 10y_tmp.nc 10ymean_YEAR_int.nc

# Hourly maximum (~40s)
cdo -s -L -timmean -yearmax 10y_tmp.nc 10ymean_YEAR_max.nc

# Percentiles (2*~40s + 4*~110s = ~520s)
cdo timmin 10y_tmp.nc 10y_tmp_timmin.nc
cdo timmax 10y_tmp.nc 10y_tmp_timmax.nc
export CDO_PCTL_NBINS=1001
cdo -L -timpctl,90 10y_tmp.nc 10y_tmp_timmin.nc 10y_tmp_timmax.nc \
10ymean_YEAR_q90.nc
cdo -L -timpctl,95 10y_tmp.nc 10y_tmp_timmin.nc 10y_tmp_timmax.nc \
10ymean_YEAR_q95.nc
cdo -L -timpctl,99 10y_tmp.nc 10y_tmp_timmin.nc 10y_tmp_timmax.nc \
10ymean_YEAR_q99.nc
cdo -L -timpctl,99.9 10y_tmp.nc 10y_tmp_timmin.nc 10y_tmp_timmax.nc \
10ymean_YEAR_q99.9.nc

# -----------------------------------------------------------------------------
