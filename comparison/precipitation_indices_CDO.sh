# Description: Compute multiple precipitation indices with CDO
#
# Author: Christian R. Steger, March 2023

# -----------------------------------------------------------------------------

# Settings (files pattern and number of time steps in concatenate files)
files_pattern="pr_EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_CLMcom-ETH-COSMO-"\
"crCLIM-v1-1_v1_1hr_????01010030-????12312330.nc"
num_ts=17544

t_beg=$(date +%s)

# Concatenate files and convert to [mm h-1]
cdo -L -mulc,3600 -cat ${files_pattern} 10y_tmp.nc

# Temporal meam
cdo -s -L -timmean -yearmean 10y_tmp.nc 10ymean_year_mean.nc

# Wet day frequency (>0.1mm/hour)
cdo -s eca_pd,0.1 10y_tmp.nc 10yfre.nc
cdo -s divc,${num_ts} 10yfre.nc 10ymean_year_fre.nc

# Intensity (mean precipitation amount on wet days)
cdo -s eca_sdii,0.1 10y_tmp.nc 10ymean_year_int.nc

# Hourly maximum
cdo -s -L -timmean -yearmax 10y_tmp.nc 10ymean_year_max.nc

# Percentiles
cdo timmin 10y_tmp.nc 10y_tmp_timmin.nc
cdo timmax 10y_tmp.nc 10y_tmp_timmax.nc
export CDO_PCTL_NBINS=1001
cdo -L -timpctl,90 10y_tmp.nc 10y_tmp_timmin.nc 10y_tmp_timmax.nc \
10ymean_year_q90.nc
cdo -L -timpctl,95 10y_tmp.nc 10y_tmp_timmin.nc 10y_tmp_timmax.nc \
10ymean_year_q95.nc
cdo -L -timpctl,99 10y_tmp.nc 10y_tmp_timmin.nc 10y_tmp_timmax.nc \
10ymean_year_q99.nc
cdo -L -timpctl,99.9 10y_tmp.nc 10y_tmp_timmin.nc 10y_tmp_timmax.nc \
10ymean_year_q99.9.nc

t_end=$(date +%s)
t_elapsed=$(( t_end - t_beg ))
echo "Total elapsed time: "${t_elapsed}"s"

# Remove temporary files
rm 10y_tmp.nc 10yfre.nc 10y_tmp_timmax.nc 10y_tmp_timmin.nc

# -----------------------------------------------------------------------------
