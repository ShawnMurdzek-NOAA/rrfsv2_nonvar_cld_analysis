"""
Plot NASA LaRC Obs Dumped in a Text File

Passed Arguments
----------------
sys.argv[1] : Input text file name dumped by larccld.fd
sys.argv[2] : Domain name

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import datetime as dt


#---------------------------------------------------------------------------------------------------
# Parameters
#---------------------------------------------------------------------------------------------------

# Input text file
infile = sys.argv[1]

# Output file tag
out_tag = sys.argv[1].split('.')[0]

# Plotting parameters
param = {'ptop': {'cmap': 'plasma',
                  'vmin': 0,
                  'vmax': 1020},
         'teff': {'cmap': 'plasma',
                  'vmin': 200,
                  'vmax': 300},
         'lwp': {'cmap': 'plasma',
                 'vmin': 0,
                 'vmax': 5},
         'phase': {'cmap': 'plasma',
                   'vmin': 0,
                   'vmax': 5},
         'frac': {'cmap': 'plasma',
                  'vmin': 0,
                  'vmax': 1},
         'nlev_cld': {'cmap': 'plasma',
                      'vmin': 0,
                      'vmax': 1}}

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
print('Starting NASA LaRC Plotting Program')
print(f"time = {start.strftime('%Y%m%d %H:%M:%S')}")
print()

# Open file
df = pd.read_csv(infile, sep="\s+")
df.replace(99999., value=np.nan, inplace=True)
df.loc[np.isclose(df['ptop'], -20), 'ptop'] = 1013.
col = list(df.columns)
col.remove('lat')
col.remove('lon')
print(f"Total entries in input file = {len(df)}")
print('Fields to plot =', col)
print()

# Create state borders
borders = cfeature.NaturalEarthFeature(category='cultural',
                                       scale='50m',
                                       facecolor='none',
                                       name='admin_1_states_provinces')

# Plot data
for c in col:
    print(f"Plotting {c}")
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal())

    cax = ax.scatter(df['lon'], df['lat'], c=df[c], s=markersize, 
                     transform=ccrs.PlateCarree(),
                     cmap=param[c]['cmap'], vmin=param[c]['vmin'], vmax=param[c]['vmax'])
    
    # Too slow
    #cax = ax.tripcolor(df['lon'], df['lat'], df[c],
    #                   transform=ccrs.PlateCarree(),
    #                   cmap=param[c]['cmap'], vmin=param[c]['vmin'], vmax=param[c]['vmax'])

    # Too slow
    #cax = ax.tricontourf(df['lon'], df['lat'], df[c], 
    #                    np.linspace(param[c]['vmin'], param[c]['vmax'], 16),
    #                    transform=ccrs.PlateCarree(),
    #                    cmap=param[c]['cmap'], extend='both')

    cbar = plt.colorbar(cax, ax=ax, orientation='horizontal')
    cbar.set_label(c, size=14)

    ax.set_extent([minlon, maxlon, minlat, maxlat])
    ax.coastlines('50m', linewidth=1, edgecolor='k')
    ax.add_feature(borders, linewidth=0.75, edgecolor='gray')
    ax.set_title(out_tag, size=18)

    plt.savefig(f"{c}_{domain}_{out_tag}.png")
    plt.close()

print()
print('Done!')
print(f"elapsed time = {(dt.datetime.now() - start).total_seconds()} s")
print()


"""
End plot_nasa_larc.py
"""
