# Precipitation_indices
Accurate and fast way to compute precipitation indices (temporal mean, maximum daily/hourly precipitation, wet day frequency and percentiles) from climate data organised in yearly blocks of NetCDF files.
Part of the percentile calculation is parallelised with Numba - the performance of the script depends thus on the available number of CPUs (respectively threads).

# Package dependencies and installation

The following Python packages are required to use this repository: NumPy, SciPy, netCDF4, Xarray, Numba.
First, create a Conda environment with all required Python packages:

```bash
conda create -n prec_percentiles -c conda-forge numpy scipy matplotlib netcdf4 xarray numba
```

and **activate this environment**. In the script **precipitation_indices.py**, the path to the script **auxiliary.py** must then be adjusted. Afterwards, **precipitation_indices.py** can be used.

The script **python_versus_CDO.py** additional requires the packages Matplotlib and cartopy for plotting. The script **submit_daint_multithread.sh** allows to submit **precipitation_indices.py** as a job on the CSCS cluster Daint.

# Comparison with different calculation methods.

The folder **comparison** contains Python scripts to compare obtained results with other methods:
- **precipitation_indices_CDO.sh**: Compute the same precipitation indices with [CDO](https://code.mpimet.mpg.de/projects/cdo/). This methods is approximately ~10 times slower than applying **precipitation_indices.py**
(depending on the available CPUs/threads). With CDO, the accuracy of the percentile computation can be improved by setting *CDO_PCTL_NBINS* to a higher value than its default one -- e.g. *export CDO_PCTL_NBINS=1001*.
- **python_versus_CDO.py**: Compare output from **precipitation_indices.py** with results from the above shell script applying CDO.
- **python_evaluate.py**: Check that the fast script **precipitation_indices.py** yields identical results to a different (but slower) method implemented in Python.

# References
- Ban, N., Caillaud, C., Coppola, E. et al. (2021): The first multi-model ensemble of regional climate simulations at kilometer-scale resolution, part I: evaluation of precipitation. Clim Dyn 57, 275–302. https://doi.org/10.1007/s00382-021-05708-w
- Chinita, M. J., Richardson, .M, Teixeira, J. and Miranda, P. M. A. (2021): Global mean frequency increases of daily and sub-daily heavy precipitation in ERA 5. Environ Res Lett 16(7):74035. https://doi.org/10.1088/1748-9326/ac0caa
- Schär, C., Ban, N., Fischer, E. M. et al. (2016): Percentile indices for assessing changes in heavy precipitation events. Climatic Change 137, 201–216. https://doi.org/10.1007/s10584-016-1669-2
