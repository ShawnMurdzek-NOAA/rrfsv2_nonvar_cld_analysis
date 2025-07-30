"""
Plot Processed NSSL Mosaic from a NetCDF File

Passed Arguments
----------------
sys.argv[1] : Input netCDF file name
sys.argv[2] : Domain name

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import sys
import xarray as xr
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import datetime as dt


#---------------------------------------------------------------------------------------------------
# Parameters
#---------------------------------------------------------------------------------------------------

# Input text file
infile = sys.argv[1]

# Plotting parameters
param = {'reflectivity': {'cmap': 'plasma',
                          'vmin': -15,
                          'vmax': 70,
                          'hgt': 3000.},
         'cref': {'cmap': 'plasma',
                  'vmin': -15,
                  'vmax': 70}}

# Plotting domain
domain = sys.argv[2]
minlat = 18
minlon = -135
maxlat = 55
maxlon = -60
markersize = 0.1
if domain == 'iowa':
    minlat = 39
    minlon = -98
    maxlat = 45
    maxlon = -88
    markersize = 3
elif domain != 'full':
    print(f"invalid domain {domain}, switching to full")


#---------------------------------------------------------------------------------------------------
# Program
#---------------------------------------------------------------------------------------------------

start = dt.datetime.now()
print('Starting NSSL Mosaic Plotting Program')
print(f"time = {start.strftime('%Y%m%d %H:%M:%S')}")
print()

# Open file
ds = xr.open_dataset(infile)

# Compute composite reflectivity
cref = np.amax(ds['reflectivity'].values, axis=0)
ds['cref'] = xr.DataArray(data=cref,
                          dims=['cell'],
                          coords={'cell':ds['cell'].values},
                          attrs={'units':ds['reflectivity'].attrs['units']})

# Create state borders
borders = cfeature.NaturalEarthFeature(category='cultural',
                                       scale='50m',
                                       facecolor='none',
                                       name='admin_1_states_provinces')

# Plot data
for c in param.keys():
    print(f"Plotting {c}")
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal())

    if 'hgt' in param[c]:
        zind = np.where(np.isclose(ds['height'].values, param[c]['hgt']))[0][0]
        data = ds[c][zind, :].values
        cbar_label = f"{param[c]['hgt']}-{ds['height'].attrs['units']} {c} ({ds[c].attrs['units']})"
        save_fname = f"{c}_{param[c]['hgt']}_{domain}_mosaic.png"
    else:
        data = ds[c].values
        cbar_label = f"{c} ({ds[c].attrs['units']})"
        save_fname = f"{c}_{domain}_mosaic.png"

    cax = ax.scatter(ds['longitude'], ds['latitude'], c=data, s=markersize,
                     transform=ccrs.PlateCarree(),
                     cmap=param[c]['cmap'], vmin=param[c]['vmin'], vmax=param[c]['vmax'])

    cbar = plt.colorbar(cax, ax=ax, orientation='horizontal') 
    cbar.set_label(cbar_label, size=14)

    ax.set_extent([minlon, maxlon, minlat, maxlat])
    ax.coastlines('50m', linewidth=1, edgecolor='k')
    ax.add_feature(borders, linewidth=0.75, edgecolor='gray')
    ax.set_title('Processed NSSL Mosaic', size=18)

    plt.savefig(save_fname)
    plt.close()

print()
print('Done!')
print(f"elapsed time = {(dt.datetime.now() - start).total_seconds()} s")
print()


"""
End plot_NSSL_mosaic_ref.py
"""
