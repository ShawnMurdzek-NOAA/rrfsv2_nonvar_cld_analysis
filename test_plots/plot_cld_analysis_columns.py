"""
Plot Columns of Processed Observations, Model Background, and Model Analysis

Passed Arguments
----------------
sys.argv[1] : YAML with input parameters

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import sys
import yaml
import pandas as pd
import numpy as np
import copy
import matplotlib.pyplot as plt
import datetime as dt
import xarray as xr
import metpy.calc as mc
from metpy.units import units


#---------------------------------------------------------------------------------------------------
# Functions
#---------------------------------------------------------------------------------------------------

def read_param(fname):
    """
    Read input YAML file
    """

    # Read input parameters
    with open(fname, 'r') as fptr:
        param = yaml.safe_load(fptr)
    
    return param


def read_nasa_larc(fname):
    """
    Read NASA LaRC observations and output as a DataFrame
    """

    df = pd.read_csv(fname, sep="\s+")
    df.replace(99999., value=np.nan, inplace=True)
    df.loc[np.isclose(df['ptop'], -20), 'ptop'] = 1013.

    return df


def read_metar_ceil(fname):
    """
    Read METAR ceilometer observations and output as a DataFrame
    """

    df = pd.read_csv(fname, sep="\s+")
    df.replace(-99999., value=np.nan, inplace=True)

    # Decode CLAM field into oktas
    # -1 indicates missing obs
    # https://www.emc.ncep.noaa.gov/mmb/data_processing/table_20.htm#0-20-011
    for i in range(1, 4):
        clam = df[f"cover{i}"].values
        oktas = -1 * np.ones(len(clam))
        oktas[clam < 9] = clam[clam < 9]
        oktas[np.isclose(clam, 9)] = 8
        oktas[np.isclose(clam, 10)] = 3
        oktas[np.isclose(clam, 11)] = 3
        oktas[np.isclose(clam, 12)] = 6
        oktas[np.isclose(clam, 13)] = 1
        df[f"okta{i}"] = oktas

    return df


def read_refl(fname):
    """
    Read MRMS reflectivty observations and output as a DataSet
    """

    ds = xr.open_dataset(fname)

    # Set missing obs to NaN for easier plotting
    ds['reflectivity'].values[np.isclose(ds['reflectivity'], -999)] = np.nan

    return ds


def read_lght(fname):
    """
    Read lightning observations and output as a DataFrame
    """

    df = pd.read_csv(fname, sep="\s+")

    return df


def extract_col(lat, lon, larc_df, ceil_df, refl_ds, lght_df, mesh_ds, bgd_ds, ana_ds):
    """
    Extract a column at the MPAS cell closest to the input (lat, lon) coordinate
    """

    # Determine MPAS cell index closest to the (lat, lon) coordinate
    # Note that the MPAS cell index from the processed ob files will be idx+1 (b/c Python
    # uses 0 as the first index)
    mpas_lat = np.rad2deg(mesh_ds['latCell'].values)
    mpas_lon = np.rad2deg(mesh_ds['lonCell'].values)
    mpas_lon[mpas_lon > 180] = mpas_lon[mpas_lon > 180] - 360
    idx = np.argmin((mpas_lon - lon)**2 + (mpas_lat - lat)**2)

    # Extract observational fields
    out_dict = {}
    out_dict['larc'] = larc_df.iloc[idx, :]
    out_dict['ceil'] = ceil_df.loc[ceil_df['icell'] == (idx+1), :]
    out_dict['refl'] = refl_ds.sel(cell=idx)
    z = out_dict['refl']['height'].values - mesh_ds['ter'][idx].values
    out_dict['refl']['z_agl'] = xr.DataArray(z, 
                                             dims=['height'], 
                                             coords={'height':out_dict['refl']['height'].values},
                                             attrs={'units':'m AGL'})
    out_dict['lght'] = lght_df.loc[lght_df['cell_id'] == (idx+1), :]

    # Extract MPAS fields
    out_dict['lat'] = mpas_lat[idx]
    out_dict['lon'] = mpas_lon[idx]
    z = mesh_ds['zgrid'][idx, :].values - mesh_ds['ter'][idx].values
    out_dict['hgt'] = 0.5*(z[1:] + z[:-1])
    out_dict['bgd'] = bgd_ds.sel(nCells=idx)
    out_dict['ana'] = ana_ds.sel(nCells=idx)
    out_dict['bgd'] = compute_mpas_rh(out_dict['bgd'])
    out_dict['ana'] = compute_mpas_rh(out_dict['ana'])

    return out_dict


def compute_mpas_rh(ds):
    """
    Compute relative humidity for a MPAS DataSet using theta, qv, pressure_p, and pressure_base
    """

    # Extract required fields
    theta = ds['theta'].values * units.K
    qv = ds['qv'].values
    p = (ds['pressure_p'].values + ds['pressure_base'].values) * units.Pa

    # Compute RH using MetPy
    T = mc.temperature_from_potential_temperature(p, theta)
    rh = mc.relative_humidity_from_mixing_ratio(p, T, qv).to('percent')

    # Save to DataSet
    ds['RH'] = copy.deepcopy(ds['theta'])
    ds['RH'].values = rh.magnitude
    ds['RH'].attrs['units'] = '%'
    ds['RH'].attrs['long_name'] = 'Relative humidity'

    return ds


def plot_col(in_dict):
    """
    Plot columns
    """

    # Configure plot
    fig, axes = plt.subplots(nrows=2, ncols=5, figsize=(12, 8), sharey=True)
    plt.subplots_adjust(hspace=0.35, wspace=0.1, left=0.07, bottom=0.09, right=0.99, top=0.89)
    lsize = 11

    # Ceilometer observations
    ax = axes[0, 0]
    for i in range(1, 4):
        if ~np.isnan(in_dict['ceil'][f"base{i}"].values[0]):
            ax.plot(in_dict['ceil'][f"okta{i}"], in_dict['ceil'][f"base{i}"], 'bo')
        if np.isclose(in_dict['ceil']['cover1'].values[0], 0):
            ax.plot(0, 0, 'bo', ms=10)
    ax.set_xlim([-1.5, 8.5])
    ax.set_xlabel('ceilometer okta', size=lsize)

    # Add height of ptop from NASA LaRC to ceilometer plot
    ptop = in_dict['larc']['ptop'] * 1e2
    p1d = in_dict['bgd']['pressure_p'].values + in_dict['bgd']['pressure_base'].values
    p1d = np.squeeze(p1d)
    ztop = np.interp([ptop], np.flip(p1d), np.flip(in_dict['hgt']))
    ax.axhline(ztop, ls='--', c='r', label='LaRC ptop')

    # Add cld_cover_3d (if it exists) to ceilometer plot
    if 'cld_cover_3d' in in_dict['ana']:
        cld_cover_3d = in_dict['ana']['cld_cover_3d'].values.T
        cld_cover_3d[np.where(cld_cover_3d < 0)] = np.nan
        cld_cover_3d = 8 * cld_cover_3d
        ax.plot(cld_cover_3d, in_dict['hgt'], 'k-', label='cld_cover_3d')
    ax.legend()

    # Reflectivity observations
    ax = axes[1, 0]
    ax.plot(in_dict['refl']['reflectivity'], in_dict['refl']['z_agl'], 'b-')
    ax.set_xlim([-5, 80])
    ax.set_xlabel(f"reflectivity ({in_dict['refl']['reflectivity'].attrs['units']})", size=lsize)

    # MPAS background, analysis, and increments: Thermodynamic fields
    fields = ['theta', 'qv', 'RH']
    for i, f in enumerate(fields):

        ax_raw = axes[0, i+1]
        ax_raw.plot(in_dict['bgd'][f].T, in_dict['hgt'], 'b-', label='bgd')
        ax_raw.plot(in_dict['ana'][f].T, in_dict['hgt'], 'r-', label='ana')
        ax_raw.set_xlabel(f"{in_dict['bgd'][f].attrs['long_name']}\n({in_dict['bgd'][f].attrs['units']})", size=lsize)
        ax_raw.legend()

        ax_inc = axes[1, i+1]
        ax_inc.plot(in_dict['ana'][f].T - in_dict['bgd'][f].T, in_dict['hgt'], 'b-')
        ax_inc.set_xlabel(f"diff {in_dict['bgd'][f].attrs['long_name']}\n({in_dict['bgd'][f].attrs['units']})", size=lsize)

    axes[0, 1].set_xlim([250, 325])

    # MPAS background, analysis, and increments: Hydrometeor fields
    colors = ['#66CCEE', '#AA3377', '#228833', '#CCBB44', '#4477AA']
    fields = ['qc', 'qi', 'qr', 'qs', 'qg']
    ax_raw = axes[0, 4]
    ax_inc = axes[1, 4]
    for i, (f, c) in enumerate(zip(fields, colors)):
        ax_raw.plot(1e3*in_dict['bgd'][f].T, in_dict['hgt'], c=c, ls='--', lw=2, label=f"{f}_b")
        ax_raw.plot(1e3*in_dict['ana'][f].T, in_dict['hgt'], c=c, ls='-', lw=1, label=f"{f}_a")
        ax_raw.set_xlabel(f"hydrometeor fields\n(g kg-1)", size=lsize)
        ax_raw.legend(ncol=2)

        ax_inc.plot(1e3*(in_dict['ana'][f].T - in_dict['bgd'][f].T), in_dict['hgt'], c=c, ls='-', label=f)
        ax_inc.set_xlabel(f"diff hydrometeor fields\n(g kg-1)", size=lsize)
        ax_inc.legend()

    # Other annotations
    larc_str = "NASA LaRC: "
    for f in ['ptop', 'teff', 'lwp', 'frac']:
        larc_str = larc_str + f"{f} = {in_dict['larc'][f]}, "
    if len(in_dict['lght']) == 0:
        lght_str = ''
    else:
        lght_str = f"nstrikes = {in_dict['lght']['nstrikes'].values[0]}"
    plt.suptitle(f"{in_dict['lat']:.3f} deg N, {in_dict['lon']:.3f} deg E\n{larc_str}\n{lght_str}", size=14)
    for i in range(2):
        axes[i, 0].set_ylabel('height (m AGL)', size=lsize)
        for j in range(5):
            axes[i, j].grid()
            axes[i, j].set_ylim([0, 4000])

    plt.savefig(f"cld_col_{in_dict['lat']:.3f}_{in_dict['lon']:.3f}.png")
    plt.close(fig)


#---------------------------------------------------------------------------------------------------
# Main Program
#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    start = dt.datetime.now()
    print('Starting plot_cld_analysis_columns.py')
    print(f"Time = {start.strftime('%Y%m%d %H:%M:%S')}")

    # Read input YAML file
    param = read_param(sys.argv[1])

    # Read observations
    print('Reading observations...')
    larc_df = read_nasa_larc(param['larc_fname'])
    ceil_df = read_metar_ceil(param['ceil_fname'])
    refl_ds = read_refl(param['refl_fname'])
    lght_df = read_lght(param['lght_fname'])

    # Read MPAS fields
    print('Reading MPAS fields...')
    mesh_ds = xr.open_dataset(param['invariant_fname'])
    bgd_ds = xr.open_dataset(param['bgd_fname'])
    ana_ds = xr.open_dataset(param['ana_fname'])

    for lat, lon in zip(param['lat'], param['lon']):
        print(f"Making plot for ({lat}, {lon})")

        # Extract column
        col_dict = extract_col(lat, lon, 
                               larc_df, ceil_df, refl_ds, lght_df, 
                               mesh_ds, bgd_ds, ana_ds)

        # Make plot
        plot_col(col_dict)

    print('Program finished!')
    print(f"Elapsed time = {(dt.datetime.now() - start).total_seconds()} s")
    

"""
End plot_cld_analysis_columns.py
"""
