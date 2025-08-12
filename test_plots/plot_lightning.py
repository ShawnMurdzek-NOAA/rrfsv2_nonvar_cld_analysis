"""
Plot Lightning Obs Dumped in a Text File

Passed Arguments
----------------
sys.argv[1] : Text file with raw lightning obs
sys.argv[2] : Text file with interpolated lightning obs
sys.argv[3] : Domain name

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

# Input text files
in_raw = sys.argv[1]
in_interp = sys.argv[2]

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
elif domain != 'full':
    print(f"invalid domain {domain}, switching to full")


#---------------------------------------------------------------------------------------------------
# Program
#---------------------------------------------------------------------------------------------------

start = dt.datetime.now()
print('Starting Lightning Plotting Program')
print(f"time = {start.strftime('%Y%m%d %H:%M:%S')}")
print()

# Open file with raw lightning obs
df_raw = pd.read_csv(in_raw, sep="\s+")
df_raw = df_raw.loc[df_raw['quality'] == 0, :]
df_raw['nstrikes'] = np.ones(len(df_raw))

# Open file with interpolated lightning obs
df_interp = pd.read_csv(in_interp, sep="\s+")

# Create state borders
borders = cfeature.NaturalEarthFeature(category='cultural',
                                       scale='50m',
                                       facecolor='none',
                                       name='admin_1_states_provinces')

# Determine bin edges for 2D histogram
#delta = 0.25
#xbin = np.arange(minlon, maxlon + 2*delta, delta)
#ybin = np.arange(minlat, maxlat + 2*delta, delta)

# Plot data
fig = plt.figure(figsize=(8, 10))
for i, (label, df) in enumerate(zip(['raw obs', 'interpolated to MPAS'], [df_raw, df_interp])):
    print(f"Plotting {label} (nlightning strikes = {df['nstrikes'].sum()})")
    ax = fig.add_subplot(2, 1, i+1, projection=ccrs.LambertConformal())

    # Create 2D histogram for lightning strikes
    #v, x, y = np.histogram2d(df['lon'], df['lat'], bins=[xbin, ybin], weights=df['nstrikes'])
    #v[np.isclose(v, 0)] = np.nan

    #cax = ax.pcolormesh(x, y, v.T, transform=ccrs.PlateCarree(),
    #                    cmap='plasma', vmin=0, vmax=4000)
    
    # Plot lightning locations with a scatter plot
    cax = plt.scatter(df['lon'], df['lat'], c=df['nstrikes'], transform=ccrs.PlateCarree(),
                      cmap='plasma', vmin=0, vmax=100, s=3, alpha=1)

    cbar = plt.colorbar(cax, ax=ax, orientation='horizontal', aspect=30, pad=0.05)
    cbar.set_label('number of lightning strikes', size=14)

    ax.set_title(label, size=16)
    ax.set_extent([minlon, maxlon, minlat, maxlat])
    ax.coastlines('50m', linewidth=1, edgecolor='k')
    ax.add_feature(borders, linewidth=0.75, edgecolor='gray')

plt.subplots_adjust(left=0.05, bottom=0.07, right=0.95, top=0.93, hspace=0.25)
plt.savefig(f"lightning_strikes_{domain}.png")
plt.close()

print()
print('Done!')
print(f"elapsed time = {(dt.datetime.now() - start).total_seconds()} s")
print()


"""
End plot_lightning.py
"""
