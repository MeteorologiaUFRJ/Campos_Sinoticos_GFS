#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 23:11:12 2022

@author: coqueiro
"""
from datetime import datetime
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from netCDF4 import num2date
import numpy as np
from siphon.catalog import TDSCatalog
import metpy.calc as mpcalc
from metpy.units import units

best_gfs = TDSCatalog(
    'https://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/catalog.html?dataset=grib/NCEP/GFS/Global_0p25deg/Best')
print(list(best_gfs.datasets))
best_ds = best_gfs.datasets[0]
ncss = best_ds.subset()
query = ncss.query()

query.lonlat_box(north=10, south=-70, east=-10, west=-120).time(datetime.utcnow())
query.accept('netcdf4')


u = query.variables('u-component_of_wind_isobaric')
v = query.variables('v-component_of_wind_isobaric')
msl = query.variables('Pressure_reduced_to_MSL_msl')
data = ncss.get_data(query)
print(list(data.variables))

u = data.variables['u-component_of_wind_isobaric']
v = data.variables['v-component_of_wind_isobaric']
msl = data.variables['Pressure_reduced_to_MSL_msl']

# Time variables can be renamed in GRIB collections. Best to just pull it out of the
# coordinates attribute on temperature
time_name = u.coordinates.split()[1]
time_var = data.variables[time_name]
lat_var = data.variables['latitude']
lon_var = data.variables['longitude']

# Get the actual data values and remove any size 1 dimensions
u = u[:]
v = v[:]
msl = msl[:]

lat_vals = lat_var[:].squeeze()
lon_vals = lon_var[:].squeeze()

# # Convert the number of hours since the reference time to an actual date
# u = num2date(u[:].squeeze(), time_var.units)
# v = num2date(v[:].squeeze(), time_var.units)
# msl = num2date(msl[:].squeeze(), time_var.units)

# Combine 1D latitude and longitudes into a 2D grid of locations
lon_2d, lat_2d = np.meshgrid(lon_vals, lat_vals)

#dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)
   
#vorticidade = mpcalc.vorticity(u, v, dx=lon_2d, dy=lat_2d, x_dim=- 1, y_dim=- 2)*10**5