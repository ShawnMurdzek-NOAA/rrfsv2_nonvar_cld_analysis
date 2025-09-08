"""
Plot Processed METAR Cloud Observations Dumped to a Text File

Passed Arguments
----------------
sys.argv[1] : Input text file name
sys.argv[2] : MPAS invariant file with 'zgrid' and 'ter' fields
sys.argv[3] : Domain name

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import sys
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import datetime as dt
import xarray as xr


#---------------------------------------------------------------------------------------------------
# Create Custom Colorbar for Max Okta
#---------------------------------------------------------------------------------------------------

plasma = mpl.colormaps['plasma_r'].resampled(9)
okta_colors = np.zeros([10, 4])
okta_colors[1:, :] = plasma(np.linspace(0, 1, 9))
okta_colors[0, :] = [120/256, 120/256, 120/256, 1]
okta_cmap = mpl.colors.ListedColormap(okta_colors)


#---------------------------------------------------------------------------------------------------
# Parameters
#---------------------------------------------------------------------------------------------------

# Input text file
infile = sys.argv[1]

# MPAS invariant file
mpas_fname = sys.argv[2]

# Plotting parameters
param = {'ceil': {'cmap': 'plasma',
                  'vmin': 0,
                  'vmax': 2500,
                  'units': 'm'},
         'max_okta': {'cmap': okta_cmap,
                      'vmin': -1.5,
                      'vmax': 8.5,
                      'units': '-1 = missing'},
         'okta': {'cmap': okta_cmap,
                  'vmin': -1.5,
                  'vmax': 8.5,
                  'units': '-1 = missing'}}

# Plotting domain
domain = sys.argv[3]
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
elif domain == 'KYTN':
    minlat = 36
    minlon = -87
    maxlat = 37
    maxlon = -85
    markersize = 5
elif domain == 'south':
    minlat = 25
    minlon = -105
    maxlat = 42
    maxlon = -79
    markersize = 1
elif domain != 'full':
    print(f"invalid domain {domain}, switching to full")


#---------------------------------------------------------------------------------------------------
# Helper Function to Create Plots
#---------------------------------------------------------------------------------------------------

def create_plot(lat, lon, c, name, param, bounds, domain, markersize):
    """
    Create a scatterplot of ceilometer obs mapped onto the MPAS mesh

    Parameters
    ----------
    lat : array
        Latitude coordinates of interpolated ceilometer obs (deg N)
    lon : array
        Longitude coordinates of interpolated ceilometer obs (deg E)
    c : array
        Field used to color the scatterplot points
    name : string
        Name of field being plotting
    param : dictionary
        Additional plotting parameters. One of the keys should match 'name'
    bounds : list
        [minlon, maxlon, minlat, maxlat]
    domain : string
        Domain name
    markersize : float
        Scatterplot marker size

    Returns
    -------
    fig : matplotlib.pyplot.figure
        Figure with scatterplot

    """

    # Create state borders
    borders = cfeature.NaturalEarthFeature(category='cultural',
                                           scale='50m',
                                           facecolor='none',
                                           name='admin_1_states_provinces')

    # Create figure and axes
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal())

    cax = ax.scatter(lon, lat, c=c, s=markersize,
                     transform=ccrs.PlateCarree(),
                     cmap=param[name]['cmap'], vmin=param[name]['vmin'], vmax=param[name]['vmax'])

    # Add colorbar
    cbar = plt.colorbar(cax, ax=ax, orientation='horizontal')
    if 'units' in param[name]:
        cbar_label = f"{name} ({param[name]['units']})"
    else:
        cbar_label = name
    cbar.set_label(cbar_label, size=14)

    # Set tick locations for oktas
    if name in ['max_okta', 'okta']:
        cbar.ax.set_xticks(np.arange(-1, 9))

    ax.set_extent(bounds)
    ax.coastlines('50m', linewidth=1, edgecolor='k')
    ax.add_feature(borders, linewidth=0.75, edgecolor='gray')
    ax.set_title('Processed METAR Cloud Obs', size=18)

    return fig


#---------------------------------------------------------------------------------------------------
# Program
#---------------------------------------------------------------------------------------------------

start = dt.datetime.now()
print('Starting METAR Cloud Plotting Program')
print(f"time = {start.strftime('%Y%m%d %H:%M:%S')}")
print()

# Open METAR text file
df = pd.read_csv(infile, sep="\s+")
df.replace(-99999., value=np.nan, inplace=True)
print(f"Total entries in input file = {len(df)}")
print()

# Open MPAS invariant file and determin vertical grid
mpas_ds = xr.open_dataset(mpas_fname)
mpas_z = np.mean(mpas_ds['zgrid'].values - mpas_ds['ter'].values[:, np.newaxis], axis=0)
nbin = len(mpas_z) - 1

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

# Compute max oktas
df['max_okta'] = np.amax(df.loc[:, [f"okta{i}" for i in range(1, 4)]].values, axis=1)

# Diagnose ceilings (CLAM = 4, 5, 6, 7, 8, 9)
ceil = -1 * np.ones(len(df))
for i in range(1, 4):
    clam = df[f"cover{i}"].values
    okta = df[f"okta{i}"].values
    hocb = df[f"base{i}"].values
    cond = (np.isclose(ceil, -1) & (okta >= 4))
    ceil[cond] = hocb[cond]
ceil[np.isclose(ceil, -1)] = np.nan
df['ceil'] = ceil

# Plot data
for name in param.keys():

    if name == 'okta':
        for j in range(nbin):
            cond = np.zeros(len(df))
            for k in range(1, 4):
                cond = np.logical_or(cond, np.logical_and(df[f'base{k}'] >= mpas_z[j], 
                                                          df[f'base{k}'] < mpas_z[j+1]))
            subset = df.loc[cond]
            max_okta = np.amax(subset.loc[:, [f"okta{i}" for i in range(1, 4)]].values, axis=1)
            if len(subset) > 0:
                print(f"Plotting {name} {j}")
                fig = create_plot(subset['MPAS_lat'], subset['MPAS_lon'], max_okta, name, param, 
                                  [minlon, maxlon, minlat, maxlat], domain, markersize)
                plt.suptitle(f'{mpas_z[j]:.1f} <= z < {mpas_z[j+1]:.1f} m', size=18)
                plt.savefig(f"{name}_{j}_{domain}_metarcld.png")
                plt.close()

    else:
        print(f"Plotting {name}")
        fig = create_plot(df['MPAS_lat'], df['MPAS_lon'], df[name], name, param, 
                          [minlon, maxlon, minlat, maxlat], domain, markersize)
        plt.savefig(f"{name}_{domain}_metarcld.png")
        plt.close()

print()
print('Done!')
print(f"elapsed time = {(dt.datetime.now() - start).total_seconds()} s")
print()


"""
End plot_metar_cld.py
"""
