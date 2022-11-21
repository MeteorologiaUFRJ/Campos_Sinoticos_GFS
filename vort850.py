#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 12:52:28 2022

@author: ladsin
"""
from datetime import datetime
import metpy.calc as mpcalc
from metpy.units import units
import numpy as np
import xarray as xr

data = xr.open_mfdataset('/work/archive/Everson/Coqueiro/CENPES/DADOS/Dados netcdf/Previsões /D/*.grib2.nc')

#seleciona as fatias de lat e lon
lon_slice = slice(-180., 180.)
lat_slice = slice(0., -50.)

lats = data.variables['lat'][:]
lons = data.variables['lon'][:]

# Criar um np.array de zeros
vort_850 = np.zeros([len(data['time']),len(data['lat']),len(data['lon'])])

for i in range(len(data['time'])):
    u_wind = units('m/s') * data['U_GRD_L100'][i,:,:].squeeze()
    v_wind = units('m/s') * data['V_GRD_L100'][i,:,:].squeeze()
    
    dx, dy =  mpcalc.lat_lon_grid_deltas(data['lon'], data['lat'])
        
    # Vai preenchendo o np.array...
    vort_850[i,:,:] = mpcalc.vorticity(u_wind, v_wind, dx=dx, dy=dy)*10**5
    
# Cria o datarray dentro da estrutura do data 
data['vort_850'] = xr.DataArray(vort_850, coords=[data['time'], data['lat'], data['lon']], dims = ['time','lat','lon'])   
data['vort_850'].attrs['name'] = 'vo'
data['vort_850'].attrs['long_name'] = 'Vorticidade 850 hPa'
data['vort_850'].attrs['units'] = '1/second'
data['lat'].attrs['units'] = 'degreees_north'
data['lon'].attrs['units'] = 'degrees_east'

# ds = xarray.Dataset({'Vorticidade':(('x', 'y'),speed)},
#                     coords = {'latitude': (('x', 'y'), lats),
#                               'longitude': (('x', 'y'), lons)},
#                     attrs={'variable':'Wind Speed'})

# Deixa somente o vort_850 dentro do data
data = data[['vort_850']]

# Salva no arquivo netcdf
data.to_netcdf('/work/archive/Everson/Coqueiro/CENPES/DADOS/dados_novos/D.nc')

print ('finished saving')