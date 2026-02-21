# ------------------------------------------------------
# 
# functions to process era5 climate data and generate
# climate input files for SCEPTER 
# 
# ------------------------------------------------------
# Need the following 1d arrays --------------
# 
# [ time ]: decimal year 
# [ temperature ]: degrees C
# [ soil moisture ]: mm / m
# 

# %% 
import os 

import icechunk 
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.style as mplstyle
import numpy as np
import pandas as pd
import xarray as xr

mplstyle.use('fast')



# %%
def add_decimal_year_coord(
    ds: xr.Dataset,
    roundtime_to: int,
    native_timecoord_name: str='time'
)->xr.Dataset:
    '''
    Compute decimal year values from datetime coordinate
    and add the values as a new coordinate in the xr
    dataset

    Parameters
    ----------
    ds : xr.Dataset
        dataset with datetime coordinates (or year coordinates)
    roundtime_to : int
        number of decimal points in decimal_year (5 usually works for hourly)
    native_timecoord_name : str
        name of the `time` coord in `ds`
    
    Returns
    -------
    xr.Dataset
        input dataset with decimal year time coordinate added
    '''
    # --- add time year coordinate
    t0 = ds[native_timecoord_name][0].values  # first time point
    if native_timecoord_name == "year":
        t0 = pd.to_datetime(f"{t0}")
        # convert time to pandas datetime index
        time_dt = pd.to_datetime(ds[native_timecoord_name].values, format='%Y')
    else:
        t0 = pd.to_datetime(t0)
        # convert time to pandas datetime index
        time_dt = pd.to_datetime(ds[native_timecoord_name].values)

    # compute elapsed time in fractional years
    years_since_start = (time_dt - t0).total_seconds() / (365.25 * 24 * 3600)
    # round 
    years_since_start = np.round(years_since_start, roundtime_to)

    # update dataset 
    ds = ds.assign_coords(time_years=(native_timecoord_name, years_since_start))

    return ds


# %%
def create_climvars_ds(
    rtds: xr.Dataset,
    site_dict: dict,
    mintime: str,
    maxtime: str,
    soilwater_layers: list,
    soiltemperature_layers: list,
    runoff_vars: list,
    roundtime_to: int=5,
    era5_unit_conversions_on: bool=True,
)-> xr.Dataset:
    '''
    Select point for each site, extract just the necessary variables, and 
    convert units. Return the resulting dataset, and the number of years 
    that the default data spans (used as an input for the 
    `save_all_case_climfiles_as_txt` function later.)

    Parameters
    ----------
    rtds : xr.Dataset
        era5 icechunked dataset (or analogous input)
    site_dict : dict
        each key is a site (e.g., 'site_1') and each value is a dict with elements for
        xlat, xlon, and name. Note that lon is 0-360.
    mintime : str
        Minimum time to pull out of data (e.g., `2000-01-01`)
    maxtime : str
        Maximum time to pull out of data (e.g., `2010-12-31`)
    soilwater_layers : list
        list of variables to use for soilwater input (generally *_layer_1, *_layer_2, etc.)
    soiltemperature_layers : list 
        list of variables to use for soil temperature input (generally *_layer_1, *_layer_2, etc.)
    runoff_vars : list
        list of variables to use for runoff input (often just `sub_surface_runoff`)
    roundtime_to : int
        number of decimal points in decimal_year (5 usually works for hourly)
    era5_unit_conversions_on : bool
        variables in `soilwater_layers`, `soiltemperature_layers`, and `runoff_vars`
        are converted following era5 to scepter conversion (see `era5_unit_conversions` function)


    '''
    # --- generate ds of these points 
    dslist = []
    for site, coords in site_dict.items():
        tmpds = rtds.sel(latitude=coords['xlat'], longitude=coords['xlon'], method='nearest').sel(time=slice(mintime, maxtime))
        tmpds = tmpds.expand_dims(site=[coords['name']]).copy()
        dslist.append(tmpds)

    dsx = xr.concat(dslist, dim='site')

    # --- pull out just the variables we want 
    allvars = [soilwater_layers, soiltemperature_layers, runoff_vars]
    allvars_list = sum(allvars, [])
    dsvar = dsx[allvars_list].copy()

    # get number of years to pass on to repeat the ltm timeseries
    dsvar = add_decimal_year_coord(dsvar, roundtime_to=roundtime_to)
    nyears_data = int(np.round(dsvar.time_years.max().values, 0))

    # --- convert units if necessary 
    if era5_unit_conversions_on:
        dsvar = era5_unit_conversions(
            dsvar,
            soilwater_layers,
            soiltemperature_layers,
            runoff_vars,
        )

    # --- create average vars across lists
    dsvar['soilwater_mm_meter'] = dsvar[soilwater_layers].to_array(dim='var').mean(dim='var')
    dsvar['runoff_mm_month'] = dsvar[runoff_vars].to_array(dim='var').mean(dim='var')
    dsvar['temperature_c'] = dsvar[soiltemperature_layers].to_array(dim='var').mean(dim='var')

    # --- return result
    return dsvar, nyears_data


def era5_unit_conversions(
    dsvar: xr.Dataset,
    soilwater_layers: list,
    soiltemperature_layers: list,
    runoff_vars: list,
)-> xr.Dataset:
    '''
    Convert era5 variable units to SCEPTER units. Soil water 
    gets converted from [m3 m-3] to [mm m-1]; temperature is
    converted from [K] to [C]; and runoff is converted from 
    [m per native time step] to [mm month-1]

    Parameters
    ----------
    dsvar : xr.Dataset
        era5 dataset with variables in the soilwater, temp, and runoff lists.
    soilwater_layers : list
        list of variables to use for soilwater input (generally *_layer_1, *_layer_2, etc.)
    soiltemperature_layers : list 
        list of variables to use for soil temperature input (generally *_layer_1, *_layer_2, etc.)
    runoff_vars : list
        list of variables to use for runoff input (often just `sub_surface_runoff`)

    Returns
    -------
    xr.Dataset 
        The input dataset with the set variables converted to new units

    '''
    # --- unit conversions
    # [ era5 soil water ]
    # ... native: m3 m-3
    # ... target: mm m-1 (note, collapsing from 3d to 1d doesn't require any operation)
    mm_per_m = 1000
    dsvar[soilwater_layers] = dsvar[soilwater_layers] * mm_per_m

    # [ era5 temperature ]
    # ... native: K
    # ... target: C
    K_to_C = -273.15
    dsvar[soiltemperature_layers] = dsvar[soiltemperature_layers] + K_to_C

    # [ era5 runoff (or infiltration) ]
    # ... native: m per native time period
    # ... target: mm month-1
    # ( note, we assume 1 month = 30.417 days )
    # get time step duration in months
    alltime = pd.to_datetime(dsvar.time.values)
    dt_days = np.diff(alltime) / np.timedelta64(1, 'D') # time step in days
    dt_months = dt_days / 30.417
    # pad the last value (repeat it) to get the original length
    dt_months = np.append(dt_months, dt_months[-1])
    # convert to data array
    dt_months_da = xr.DataArray(dt_months, coords={'time': dsvar.time}, dims='time')

    # scale runoff to month-1
    for var in runoff_vars: # must loop because we're scaling with a data array
        dsvar[var] = dsvar[var] / dt_months_da

    # --- return result
    return dsvar 



def create_dsdict_across_timeResolutions(
    dsvar: xr.Dataset,
    default_res_name: str,
    default_resolution: bool,
    daily_mean: bool, 
    daily_ltm: bool,
    monthly_mean: bool,
    monthly_ltm: bool,
    annual_mean: bool,
    annual_ltm: bool,
    roundtime_to: int = 5,
    start_at_zero: bool=True,
)->dict:
    '''
    Compute climate for different time resolutions and put 
    each resulting dataset in a dictionary. There are 5 different
    time-resolution options [default_resolution, monthly_mean, monthly_ltm,
    annual_mean, annual_ltm]

    Parameters
    ----------
    dsvar : xr.Dataset
        Input climate dataset. Time coord should be in `datetime`. No restrictions
        on variables. 
    default_res_name : str
        Name for the native time resolution (used as the key for the default_resolution
        ds in the output dictionary). This is ignored if default_resolution=False.
    default_resolution : bool
        True to keep the data at the default resolution.
    daily_mean : bool
        True to keep the data averaged at each day 
    daily_ltm : bool
        True to keep the long-term daily mean of the data (length = 365)
    monthly_mean : bool
        True to keep the data averaged at each month (length = 12 * nyears)
    monthly_ltm : bool
        True to keep the long-term monthly mean of the data (length = 12)
    annual_mean : bool
        True to keep the data averaged at each year (length = nyears)
    annual_ltm : bool
        True to keep the long-term annual mean of the data (length = 1)
    rountime_to : int
        number of decimal points in decimal_year (5 usually works for hourly)
    start_at_zero : bool
        whether to make the `time_years` start at zero. If false,
        time_years starts at 0+dt
    
    Returns
    -------
    dict 
        dictionary of xr.Datasets, each dataset being the same climate data at 
        some different time resolution. 
    '''
    # --- create datasets for each time resolution 
    ds_dict = {}

    # [ default resolution ]
    if default_resolution:
        if not start_at_zero: # shift values by one dt
            dt = np.diff(dsvar.time_years.values)[0]
            dsvar = dsvar.assign_coords(time_years=dsvar['time_years'] + dt)
        ds_dict[default_res_name] = dsvar.copy()


    # [ daily ]
    if daily_mean:
        tmpds_daily = dsvar.resample(time='D').mean()
        # update time_years and save
        if not start_at_zero:
            tmpds_daily2 = add_decimal_year_coord(tmpds_daily, roundtime_to=roundtime_to)
            dt = np.diff(tmpds_daily2.time_years.values)[0]
            ds_dict['daily'] = tmpds_daily2.assign_coords(time_years=tmpds_daily2['time_years'] + dt)
        else:
            ds_dict['daily'] = add_decimal_year_coord(tmpds_daily, roundtime_to=roundtime_to)

    # [ daily ltm ]
    if daily_ltm:
        tmpds_daily = dsvar.resample(time='D').mean()
        # then get ltm
        tmpds_daily_ltm = tmpds_daily.groupby("time.dayofyear").mean("time").copy()
        # update time_years (which was lost from the mean)
        time_years_daily_ltm = tmpds_daily_ltm.dayofyear.values/365 - (1/365)  # (time starts at zero)
        tmpds_daily_ltm = tmpds_daily_ltm.assign_coords(time_years=('dayofyear', time_years_daily_ltm))  # assign coords
        if not start_at_zero:
            dt = np.diff(tmpds_daily_ltm.time_years.values)[0]
            tmpds_daily_ltm = tmpds_daily_ltm.assign_coords(time_years=tmpds_daily_ltm['time_years'] + dt)
        # add to dict
        ds_dict['daily_ltm'] = tmpds_daily_ltm

    # [ monthly ]
    if monthly_mean:
        tmpds_monthly = dsvar.resample(time='1MS').mean()
        # update time_years and save (lost from the mean)
        if not start_at_zero:
            tmpds_monthly2 = add_decimal_year_coord(tmpds_monthly, roundtime_to=roundtime_to)
            dt = np.diff(tmpds_monthly2.time_years.values)[0]
            ds_dict['monthly'] = tmpds_monthly2.assign_coords(time_years=tmpds_monthly2['time_years'] + dt)
        else: # (make decimal year outside of `add_decimal_year_coord` functionality because we want even spacing)
            time = tmpds_monthly['time']
            time_years_in = (time.dt.year - time.dt.year[0] + (time.dt.month-1)/12)
            ds_dict['monthly'] = tmpds_monthly.assign_coords(time_years=time_years_in)

    # [ monthly ltm ]
    # this appears to return evenly spaced timesteps even though we start from months
    if monthly_ltm:
        # clim_repeated = xr.DataArray(      # to plot 1-12 over and over for the dsx['time'] period
        #     data=dsvar['time.month'].values,
        #     coords={'time': dsvar['time']},
        #     dims='time'
        # )
        tmpds_monthly_ltm = dsvar.groupby('time.month').mean(dim='time').copy()
        # update time_years (which was lost from the mean)
        time_years_monthly_ltm = (tmpds_monthly_ltm.month.values / 12) - (1/12)  # (time starts at zero)
        tmpds_monthly_ltm = tmpds_monthly_ltm.assign_coords(time_years=('month', time_years_monthly_ltm))  # assign coords
        if not start_at_zero:
            dt = np.diff(tmpds_monthly_ltm.time_years.values)[0]
            tmpds_monthly_ltm = tmpds_monthly_ltm.assign_coords(time_years=tmpds_monthly_ltm['time_years'] + dt)
        ds_dict['monthly_ltm'] = tmpds_monthly_ltm

    # [ annual ]
    if annual_mean:
        tmpds_yearly = dsvar.groupby('time.year').mean(dim='time')
        # update time_years and save (lost from the mean)
        if not start_at_zero:
            tmpds_yearly2 = add_decimal_year_coord(tmpds_yearly, roundtime_to=roundtime_to, native_timecoord_name='year')
            dt = 0.08333
            ds_dict['yearly'] = tmpds_yearly2.assign_coords(time_years=tmpds_yearly2['time_years'] + dt)
        else: # (do it outside of `add_decimal_year_coord` because we want evenly spaced decimal year grid)
            time = tmpds_yearly['year']
            time_years_in = xr.DataArray(np.arange(time.size), dims="year", coords={"year": time})
            ds_dict['yearly'] = tmpds_yearly.assign_coords(time_years=time_years_in)
            # ds_dict['yearly'] = add_decimal_year_coord(tmpds_yearly, roundtime_to=roundtime_to, native_timecoord_name='year')

    # [ annual ltm ]
    if annual_ltm:
        tmpds_yearly_ltm = dsvar.groupby('time.year').mean(dim='time').mean(dim='year')
        # add time_year with monthly resolution (just the annual mean repeated each month)
        ds_dict['yearly_ltm'] = tmpds_yearly_ltm.expand_dims(time_years = list(np.linspace(0, 1-0.08333,12)))
    
    return ds_dict 


# %% 
# ---- [ FUNCTION TO WRITE OUTPUT FILE GIVEN ARRAYS ]
def write_output(outdir, var, arrvar, arrtime):
    # file_path = os.path.join(outdir, var.filename)
    # detect if writing to s3 or local
    if outdir.startswith("s3://"):
        import s3fs
        fs = s3fs.S3FileSystem()
        file_path = f"{outdir.rstrip('/')}/{var.filename}"
        f = fs.open(file_path, "w")
    else:
        os.makedirs(outdir, exist_ok=True)
        file_path = os.path.join(outdir, var.filename)
        f = open(file_path, "w")

    # open the file in write mode
    with f as file:
        # write header
        h1 = "# " + var.colname_time
        h2 = var.colname_var
        file.write(f"{h1}\t{h2}\n")
        # write each pair of values from the arrays to the file
        for value1, value2 in zip(arrtime, arrvar):
            formatted_value1 = f"{value1:.7f}"
            formatted_value2 = f"{value2:.7f}"
            file.write(
                f"{formatted_value1}\t{formatted_value2}\n"
            )  # Separate values by tab and end line with newline character

# ---- [ DEFINE VAR CLASS ]
class Var:
    def __init__(self, varname, filename, dataset_loc, colname_var, colname_time):
        self.varname = varname
        self.filename = filename
        self.dataset_loc = dataset_loc     # not used for era5
        self.colname_var = colname_var
        self.colname_time = colname_time


def save_clim_plot(
    climpath: str,
    savehere: str,
    savename: str='clim-timeseries.png',
):
    '''
    Read in the three clim `.in` files and plot the timeseries as a 
    three-panel figure. Save the figure. Climate files should be named
    `q_temp.in`; `T_temp.in`; and `Wet_temp.in` following SCEPTER convention.

    Parameters
    ----------
    climpath : str
        Path to the three climate files
    savehere : str
        Path to where the figure gets saved
    savename : str
        Name of the figure that gets saved
    
    Returns 
    -------
    None
    '''
    # ----------------------
    # plot the results and save in the dir we just created 
    # get file names
    f1, f2, f3 = "q_temp.in", "T_temp.in", "Wet_temp.in"

    # Read the tab-delimited .txt file into a DataFrame
    df1 = pd.read_csv(os.path.join(climpath, f1), delimiter="\t")
    df2 = pd.read_csv(os.path.join(climpath, f2), delimiter="\t")
    df3 = pd.read_csv(os.path.join(climpath, f3), delimiter="\t")

    # --- plot 
    fig = plt.figure(figsize = (8,12))
    gs = gridspec.GridSpec(nrows=3, ncols=1, figure=fig)

    # initialize subplots
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[2, 0])

    # [ plot 1 ]
    ax1.plot(df1.iloc[:, 0], df1.iloc[:, 1], lw=2, c='steelblue')
    ax1.set_xlabel(df1.columns[0])
    ax1.set_ylabel(df1.columns[1])

    ax2.plot(df2.iloc[:, 0], df2.iloc[:, 1], lw=2, c='crimson')
    ax2.set_xlabel(df2.columns[0])
    ax2.set_ylabel(df2.columns[1])

    ax3.plot(df3.iloc[:, 0], df3.iloc[:, 1], lw=2, c='navy')
    ax3.set_xlabel(df3.columns[0])
    ax3.set_ylabel(df3.columns[1])

    # save figure 
    if savehere.startswith("s3://"):
        import s3fs
        fs = s3fs.S3FileSystem()
        with fs.open(os.path.join(savehere, savename), "wb") as f:
            plt.savefig(f, dpi=300, format="png")  # specify format since using file handle
    else:
        plt.savefig(os.path.join(savehere, savename), dpi=300)
    plt.close();


def save_all_case_climfiles_as_txt(
    ds_dict: dict,
    nyears_data: int,
    save_maindir: str,
    inputvar_details_fn: str,
    save_clim_fig: bool=True,
    subdir_rule: str="combined",
    ltm_repeat: bool=False,
):
    '''
    Save SCEPTER climate files for a series of cases, looping over 
    each site in every dataset, and over each time resolution in the
    dictionary of datasets (ds_dict). Files get saved at `save_maindir/
    sitename/key/` or `save_maindir/sitename_key` depending on subdir_rule. 

    Parameters
    ----------
    ds_dict : dict
        dictionary of datasets used to generate the files. Each must have column 
        names listed in the `varnames` column of the `inputvar_details.txt` file.
    nyears_data : int
        number of years of the default resolution data (used to repeat the long-
        term mean data by the same amount)
    save_maindir : str
        Where to save the resulting data (excluding the sitename and key subdirs).
        Often `savepath/era5_{mintime}_{maxtime}`.
    inputvar_details_fn : str
        Path to the `inputvar_details.txt` file relating the dataset variable 
        names to the respective climate filenames / headers / etc. 
    saveplots : bool
        Whether to plot and save the climate data as we create it 
    subdir_rule : str
        [separate | combined] separate means the subdirs are split 
        (`save_maindir/sitename/key`) combined gets you: (`save_maindir/sitename_key`)
    ltm_repeat : bool
        [True | False] if True, then repeat the ltm data for however many years of input 
        data we use (then multiplied by `nyears_repeat`). 
        If False, then repeat the ltm data once (then multiplied by `nyears_repeat`).
    
    Returns
    -------
    None
    
    '''
    # --- loop over ds_dict 
    for key, tmpds in ds_dict.items():

        # --- loop over sites 
        for sitename, tmpds_site in tmpds.groupby('site'):
            # [ track progress ]
            print(f'Now saving: {sitename} -- {key}')
            # [ remove the site dimension ]
            tmpds_site = tmpds_site.squeeze('site')
            # [ save sub sub dir ]
            if subdir_rule == 'combined':
                save_subdir = f'{sitename}_{key}'
            else:
                save_subdir = os.path.join(sitename, key)
            savehere = os.path.join(save_maindir, save_subdir)
            if not os.path.exists(savehere):
                os.makedirs(savehere)
            
            # [ loop through variables ]
            with open(inputvar_details_fn) as file:
                # skip the first line
                next(file)
                # loop through each line in the file
                for line in file:
                    # create var object
                    thisvar = Var(*line.strip().split(","))

                    # [ make arrays ]
                    # if key == "yearly_ltm":
                    #     arrtime = np.array([0])
                    # else:
                    arrtime = tmpds_site.time_years.values
                    arrvar = tmpds_site[thisvar.varname].values
                    # repeat arrays if `ltm`
                    if (ltm_repeat) and ('_ltm' in key):
                        arrtime = np.concatenate([arrtime + i for i in range(nyears_data)])
                        arrvar = np.tile(arrvar, nyears_data)

                    # [ write output ]
                    write_output(savehere, thisvar, arrvar, arrtime)

            # [ plot + save ]
            if save_clim_fig:
                save_clim_plot(
                    climpath = savehere,
                    savehere = savehere,
                    savename = 'clim-timeseries.png',
                )