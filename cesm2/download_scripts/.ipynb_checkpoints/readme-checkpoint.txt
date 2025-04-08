# Data download
1. Navigate to https://www.earthsystemgrid.org/dataset/ucar.cgd.cesm2le.output.html
2. Select the child dataset of interest (we used monthly averages in all cases)
3. Select the variable of interest (TS, QFLX, PRECT, etc...)
4. Select `Download Options`
5. In Filter by File Name, enter `*cmip6*1231.001*h0`; which restricts to the first ensemble member (001) of the runs initialized in year 1231 that follow CMIP6 (rather than the biomass burning cases).
5a. to repeat for all ensemble members starting at year 1231 use `*cmip6*1231*h0`
6. Select all --> Download Options --> Python script. 