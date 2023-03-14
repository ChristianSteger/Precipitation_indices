Accurate and fast way to compute the **precipitation indices**
- temporal mean
- maximum yearly daily/hourly precipitation
- wet day frequency
- percentiles

from climate data organised in yearly blocks of NetCDF files. Daily and hourly temporal frequency can be processed and percentiles can be computed as **all-day/hour** or **wet-day/hour** according to Schär et al. (2016).
Indices can be derived for **entire years** or seasons (**JJA**, **SON**, **DJF**, **MAM**).

To reduce memory consumption, yearly blocks are processed sequentially and the relevant data is kept in memory. To e.g. compute the 90% percentile, the largest 10% of data of a grid cell is kept in memory
and sequentially updated. After all yearly blocks have been processed, the percentile is computed via interpolation. This method is identical to apply the function *numpy.percentile()*. Part of the percentile
calculation is parallelised with [Numba](http://numba.pydata.org) – the performance of the script depends thus on the available number of CPUs (respectively threads).

# Dependencies and installation

The following Python packages are required to use this repository: NumPy, SciPy, netCDF4, Xarray, Numba.
First, create a Conda environment with all required Python packages:

```bash
conda create -n prec_percentiles -c conda-forge numpy scipy matplotlib netcdf4 xarray numba
```

and **activate this environment**. In the script **precipitation_indices.py**, the path to the script **auxiliary.py** must then be adjusted. Afterwards, **precipitation_indices.py** can be used.

The script **python_versus_CDO.py** additional requires the packages Matplotlib and cartopy for plotting. The script **submit_daint_multithread.sh** allows to submit **precipitation_indices.py** as a job on the CSCS cluster Daint.

# Example data

The example data for the script can be downloaded from [ESGF](https://esgf-data.dkrz.de/search/cordex-dkrz/).
Select the following options: Project &rarr; CORDEX, Domain &rarr; EUR-11, Experiment &rarr; evaluation, RCM Model &rarr; COSMO-crCLIM-v1-1, Time Frequency &rarr; 1hr, Variable &rarr; pr. Then download the first 10 files:

*pr_EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_1hr_197901010030-197912312330.nc*\
..\
*pr_EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_1hr_198801010030-198812312330.nc*\

# Comparison with different calculation methods

The folder **comparison** contains Python scripts to compare obtained results with other methods:

- **precipitation_indices_CDO.sh**: Compute the same precipitation indices with [CDO](https://code.mpimet.mpg.de/projects/cdo/). This methods is approximately ~7 times slower than applying **precipitation_indices.py**
(depending on the available CPUs/threads). With CDO, percentiles are estimated by sorting the data into bins (see [CDO Manual](https://code.mpimet.mpg.de/projects/cdo/embedded/cdo.pdf)) to reduce memory consumption.
By default, 101 bins are applied. To obtain more accurate results, the bin size can be increase with e.g. *export CDO_PCTL_NBINS=1001*.

- **python_versus_CDO.py**: Compare output from **precipitation_indices.py** with results from the above shell script applying CDO.
The below plot shows the result for the eight precipitation indices where the output from **precipitation_indices.py** is shown in the left and the difference to CDO in the right column.
The comparably larger (and non-random) deviations for *temporal mean* is caused by the slight difference in computing averages over the 10 years (in the Python script, the average is weighted by the number of days in a year).
For the *yearly maximum*, the *wet day frequency* and the *intensity*, the deviations are very small (and random).
For *percentiles*, the deviations are considerably larger and particularly occur in the southern part of the domain (at least for the percentiles in the range of 90% - 99%).
The deviations in the percentiles are caused by the different methods applied to compute them (Python: consider all individual values, CDO: bin data in 1001 bins).

![Alt text](https://github.com/ChristianSteger/Media/blob/master/Precipitation_indices_Python_vs_CDO.png?raw=true "Output from python_versus_CDO.py")

- **python_evaluate.py**: Check that the fast script **precipitation_indices.py** yields identical results to a different (but more memory-intensive) method implemented in Python that computes indices for a time series in one step.

# References
- Ban, N., Caillaud, C., Coppola, E. et al. (2021): The first multi-model ensemble of regional climate simulations at kilometer-scale resolution, part I: evaluation of precipitation. Clim Dyn 57, 275–302. https://doi.org/10.1007/s00382-021-05708-w
- Chinita, M. J., Richardson, .M, Teixeira, J. and Miranda, P. M. A. (2021): Global mean frequency increases of daily and sub-daily heavy precipitation in ERA 5. Environ Res Lett 16(7):74035. https://doi.org/10.1088/1748-9326/ac0caa
- Schär, C., Ban, N., Fischer, E. M. et al. (2016): Percentile indices for assessing changes in heavy precipitation events. Climatic Change 137, 201–216. https://doi.org/10.1007/s10584-016-1669-2
