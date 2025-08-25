"""
Plot Processed METAR Cloud Observations Dumped to a Text File

Passed Arguments
----------------
sys.argv[1] : Input text file name
sys.argv[2] : Domain name

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

# Plotting parameters
param = {'ceil': {'cmap': 'plasma',
                  'vmin': 0,
                  'vmax': 7000,
                  'units': 'm'},
         'max_okta': {'cmap': okta_cmap,
                      'vmin': -1.5,
                      'vmax': 8.5,
                      'units': '-1 = missing'}}

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
# Program
#---------------------------------------------------------------------------------------------------

start = dt.datetime.now()
print('Starting METAR Cloud Plotting Program')
print(f"time = {start.strftime('%Y%m%d %H:%M:%S')}")
print()

# Open file
df = pd.read_csv(infile, sep="\s+")
df.replace(-99999., value=np.nan, inplace=True)
print(f"Total entries in input file = {len(df)}")
print()

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

    cax = ax.scatter(df['MPAS_lon'], df['MPAS_lat'], c=df[c], s=markersize,
                     transform=ccrs.PlateCarree(),
                     cmap=param[c]['cmap'], vmin=param[c]['vmin'], vmax=param[c]['vmax'])

    cbar = plt.colorbar(cax, ax=ax, orientation='horizontal')
    if 'units' in param[c]:
        cbar_label = f"{c} ({param[c]['units']})"
    else:
        cbar_label = c
    cbar.set_label(cbar_label, size=14)

    # Set tick locations for max_oktas
    if c == 'max_okta':
        cbar.ax.set_xticks(np.arange(-1, 9))

    ax.set_extent([minlon, maxlon, minlat, maxlat])
    ax.coastlines('50m', linewidth=1, edgecolor='k')
    ax.add_feature(borders, linewidth=0.75, edgecolor='gray')
    ax.set_title('Processed METAR Cloud Obs', size=18)

    plt.savefig(f"{c}_{domain}_metarcld.png")
    plt.close()

print()
print('Done!')
print(f"elapsed time = {(dt.datetime.now() - start).total_seconds()} s")
print()


"""
End plot_metar_cld.py
"""
