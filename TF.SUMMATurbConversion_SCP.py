# ------------------------------------------------------------------------------
# netcdf/numpy/xray/stats
import numpy as np
from datetime import datetime, timedelta
import pandas as pd
import xarray as xray

# ------------------------------------------------------------------------------
# OS interaction
import sys
import pickle
import os

# ------------------------------------------------------------------------------
# plotting packages
import seaborn as sns
sns.set_style("whitegrid")
sns.set_context('talk')
import matplotlib
from matplotlib.pyplot import subplots
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.basemap import Basemap

# ------------------------------------------------------------------------------
# Custom packages
import turbpy
import turbpy.multiConst as mc

# ------------------------------------------------------------------------------
# Directory Lists
# ------------------------------------------------------------------------------
# Unix - This isn't my desktop! These paths are for Jessica's unix box
if 'linux' in sys.platform:
    dir_pre = '/home/lapok/'
# Mac
elif 'darwin' in sys.platform:
    dir_pre = '/Users/karllapo/gdrive/'
# Paths - uniform between platforms
dirProj = dir_pre + 'SnowHydrology/proj/TurbulenceFeedbacks/'
dirData = dir_pre + 'SnowHydrology/proj/TurbulenceFeedbacks/data/SCP'

# Test that all paths work
os.chdir(dirProj)
os.chdir(dirData)

# ------------------------------------------------------------------------------
# Open netcdfs
# ------------------------------------------------------------------------------
os.chdir(dirData)
Atower = xray.open_dataset('SCP.Atower.netcdf')
Ctower = xray.open_dataset('SCP.Ctower.netcdf')
Mtower = xray.open_dataset('SCP.Mtower.netcdf')
soilObs = xray.open_dataset('SCP.soil_obs.netcdf')
radObs = xray.open_dataset('SCP.rad_obs.netcdf')
fluxObs = xray.open_dataset('SCP.flux_obs.netcdf')

# ------------------------------------------------------------------------------
# Load 30min averages
# ------------------------------------------------------------------------------
os.chdir(dirData)
fluxObs_30m = xray.open_dataset('SCP.flux_obs_30m.netcdf')
radObs_30m = xray.open_dataset('SCP.rad_obs_30m.netcdf')
soilObs_30m = xray.open_dataset('SCP.soil_obs_30m.netcdf')

# ------------------------------------------------------------------------------
# Surface temperature from upwelling longwave
# ------------------------------------------------------------------------------
sigma = 5.67 * 10 ** (-8.)                   # Stefan-Boltzmann constant
Tsfc = (radObs.Rlw_out / sigma) ** (1. / 4.)
# Insert into Mtower xarray.Dataset
Mtower['Tsfc'] = (('time'), Tsfc - 273.15)

# ------------------------------------------------------------------------------
# Bulk Richardson Number
# ------------------------------------------------------------------------------
Mtower['UBar_15m'] = (('time'), (Mtower.U_15m_M ** 2.
                                 + Mtower.V_15m_M**2) ** (1. / 2.)
                      )
Mtower.UBar_15m[Mtower.UBar_15m > 10] = np.nan

Mtower['UBar_2m'] = (('time'), (Mtower.u_2m_M ** 2.
                                + Mtower.v_2m_M ** 2.) ** (1. / 2.)                     )
Mtower.UBar_15m[Mtower.UBar_2m > 10] = np.nan

RiBulk_15m, _, _ = turbpy.bulkRichardson(Mtower.T_15m_M + 273.15,
                                         Mtower.Tsfc + 273.15,
                                         Mtower.UBar_15m,
                                         15.)
Mtower['RiBulk_15m'] = RiBulk_15m
RiBulk_2m, _, _ = turbpy.bulkRichardson(Mtower.T_1m_M + 273.15,
                                        Mtower.Tsfc + 273.15,
                                        Mtower.UBar_2m,
                                        15.)
Mtower['RiBulk_2m'] = RiBulk_2m

# Offline Turbulence
# Met variables
scalarGroundSnowFraction = 1.
soilRelHumidity = 1.
airPress = 101000
# Turbulence parameters
z0Ground = .005
# Control variables
#ixStability = ('standard',
#               'louisInversePower',
#               'mahrtExponential',
#               'moninObukhov')

ixStability = ('moninObukhov',
				'standard',
               	'louisInversePower',
               	'mahrtExponential')
# Select stable conditions only for offline comparison
# 15m and 1m Ri_Bulk values have strongly overlapping stable values
ind = np.nonzero((Mtower.RiBulk_15m > 0.).values)

sensible_1m = xray.Dataset()
sensible_1m.coords['time'] = Mtower.time[ind]
latent_1m = xray.Dataset()
latent_1m.coords['time'] = Mtower.time[ind]

sensible_15m = xray.Dataset()
sensible_15m.coords['time'] = Mtower.time[ind]
latent_15m = xray.Dataset()
latent_15m.coords['time'] = Mtower.time[ind]

for stab in ixStability:
    senHeatGround_15m = np.ones(Mtower.time[ind].size) * np.nan
    latHeatGround_15m = np.ones(Mtower.time[ind].size) * np.nan
    senHeatGround_1m = np.ones(Mtower.time[ind].size) * np.nan
    latHeatGround_1m = np.ones(Mtower.time[ind].size) * np.nan

    for n, d in enumerate(Mtower.time[ind]):
        if (d['time.hour'].values == 0) & (d['time.minute'].values == 2):
            print(stab + ": " +
                  str(d['time.month'].values) + "/" +
                  str(d['time.day'].values) + "/" +
                  str(d['time.year'].values) +
                  " (" + str(int(n / Mtower.time.size * 100.)) + "%)")
        ds = Mtower.sel(time=d)
        snowDepth = 0.  # (m)

        # "1"m variables
        airTemp_1m = ds.T_1m_M + 273.15  # (C) -> (K)
        svp, _ = turbpy.conversionTools.satVapPress(airTemp_1m)  # Saturation vapor pressure
        airVaporPress_1m = svp * ds.RH_1m_M / 100.  # Vapor pressure (Pa)
        windspd_1m = ds.UBar_2m  # (m/s)

        # "10"m variables
        airTemp_15m = ds.T_15m_M + 273.15  # (C) -> (K)
        svp, _ = turbpy.conversionTools.satVapPress(airTemp_15m)  # Saturation vapor pressure
        airVaporPress_15m = svp * ds.RH_8m_M / 100.  # Vapor pressure (Pa)
        windspd_15m = ds.UBar_15m  # (m/s)

        # Surface variables
        sfcTemp = ds.Tsfc + 273.15  # (C) -> (K)
        soilRelHumidity = soilObs.Qsoil_c.sel(time=d) / 100.  # (fraction)
        sfcVaporPress, _ = turbpy.conversionTools.satVapPress(sfcTemp)
        sfcVaporPress = sfcVaporPress * soilRelHumidity

        # 1m
        if np.any(np.isnan([snowDepth, airTemp_1m, sfcTemp,
                            airVaporPress_1m, sfcVaporPress,
                            windspd_1m
                            ]
                           )
                  ):
            continue
        else:
            # Offline turbulence - 1m
            mHeight = 1.
            (_, _, senHeatGround_1m[n], latHeatGround_1m[n],
             _, _, _, _
             ) = turbpy.turbFluxes(airTemp_1m, airPress, airVaporPress_1m,
                                   windspd_1m, sfcTemp, sfcVaporPress,
                                   snowDepth, mHeight, groundSnowFraction=1,
                                   ixStability=stab, z0Ground=.005,
                                   )

        # 15m
        if np.any(np.isnan([snowDepth, airTemp_15m, sfcTemp,
                            airVaporPress_15m, sfcVaporPress, windspd_15m
                            ]
                           )
                  ):
            continue
        else:
            mHeight = 15.
            # Offline turbulence - 15m
            (_, _, senHeatGround_15m[n], latHeatGround_15m[n],
             _, _, _, _
             ) = turbpy.turbFluxes(airTemp_15m, airPress, airVaporPress_15m,
                                   windspd_15m, sfcTemp, sfcVaporPress,
                                   snowDepth, mHeight, groundSnowFraction=1,
                                   ixStability=stab, z0Ground=.005
                                   )

    sensible_15m[stab] = (('time'), senHeatGround_15m)
    latent_15m[stab] = (('time'), latHeatGround_15m)

    sensible_1m[stab] = (('time'), senHeatGround_1m)
    latent_1m[stab] = (('time'), latHeatGround_1m)

# Save offline sensible and latent heat fluxes
# All is appended to preserve previous version (5% of total data) in case
# of potential flaws
os.chdir(dirProj)
sensible_15m.to_netcdf('OfflineTurb.SCP.sensible_15m.ALL.nc')
sensible_1m.to_netcdf('OfflineTurb.SCP.sensible_1m.ALL.nc')
latent_15m.to_netcdf('OfflineTurb.SCP.latent_15m.ALL.nc')
latent_1m.to_netcdf('OfflineTurb.SCP.latent_1m.ALL.nc')
