#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 16:52:10 2022

@author: coqueiro
"""

#importando bibliotecas
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import metpy.calc as mpcalc
from metpy.units import units
import numpy as np
import xarray as xr
import cartopy.io.shapereader as shpreader # Import shapefiles
from datetime import datetime, timedelta  # basicas datas e tipos de tempo
import cmocean

#dataset
file_1 = xr.open_dataset(
    '/home/coqueiro/Downloads/GFS_Global_0p25deg_20220903_1200.grib2.nc4'
    ).metpy.parse_cf()

file_1 = file_1.assign_coords(dict(
    longitude = (((file_1.longitude.values + 180) % 360) - 180))
    ).sortby('longitude')

#seleciona as fatias de lat e lon
lon_slice = slice(-120., 10.)
lat_slice = slice(10., -70.)

#pega as lat/lon
lats = file_1.latitude.sel(latitude=lat_slice).values
lons = file_1.longitude.sel(longitude=lon_slice).values

#seleciona o nivel em hPa
level = 1000 * units('hPa')

for i in range(len(file_1.variables['time1'])):
    
    
    u = file_1['u-component_of_wind_isobaric'].metpy.sel(
        time = file_1.time1[i], 
        vertical=level, 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze()
    
    v = file_1['v-component_of_wind_isobaric'].metpy.sel(
        time = file_1.time1[i], 
        vertical=level, 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze()
    
    pnmm = file_1.Pressure_reduced_to_MSL_msl.metpy.sel(
        time = file_1.time1[i], 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze()*0.01*units.hPa/units.Pa
    
    #data
    vtime = file_1.time1.data[i].astype('datetime64[ms]').astype('O')
    
    dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)
    
    vorticidade = mpcalc.vorticity(u, v, dx=dx, dy=dy, x_dim=- 1, y_dim=- 2)*10**5
    
    # escolha o tamanho do plot em polegadas (largura x altura)
    plt.figure(figsize=(25,25))
    
    # usando a projeção da coordenada cilindrica equidistante 
    ax = plt.axes(projection=ccrs.PlateCarree())
    
    
    shapefile = list(
        shpreader.Reader(
        '/home/coqueiro/ufrj/Estagio_supervisionado/shapefile/unidades_federativas/Brasil/UFEBRASIL.shp'
        ).geometries()
        )
    
    ax.add_geometries(
        shapefile, 
        ccrs.PlateCarree(), 
        edgecolor = 'black', 
        facecolor='none', 
        linewidth=0.5
        )
    
    # adiciona continente e bordas
    ax.coastlines(resolution='50m', color='black', linewidth=1)
    ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=1)
    gl = ax.gridlines(crs=ccrs.PlateCarree(),
                      color='gray',
                      alpha=1.0, 
                      linestyle='--', 
                      linewidth=0.5,
                      xlocs=np.arange(-180, 180, 10), 
                      ylocs=np.arange(-90, 90, 10), 
                      draw_labels=True
                      )
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 29, 'color': 'black'}
    gl.ylabel_style = {'size': 29, 'color': 'black'}
    
    # intevalos da pnmm
    intervalo_min2 = np.amin(np.array(pnmm))
    intervalo_max2 = np.amax(np.array(pnmm))
    interval_2 = 2              # de quanto em quanto voce quer que varie
    levels_2 = np.arange(intervalo_min2, intervalo_max2, interval_2)
    
    # intevalos da divergencia - umidade
    intervalo_min3 = -30
    intervalo_max3 = 1
    interval_3 = 2             # de quanto em quanto voce quer que varie
    levels_3 = np.arange(intervalo_min3, intervalo_max3, interval_3)
    
    # adiciona mascara de terra
    ax.add_feature(cfeature.LAND)
    
    # plota a imagem geopotencial
    sombreado = ax.contourf(lons, 
                            lats, 
                            vorticidade, 
                            cmap=cmocean.cm.dense_r, 
                            levels = levels_3, 
                            extend = 'both'
                            )
    
    # plota a imagem pressao
    contorno = ax.contour(lons,
                          lats, 
                          pnmm, 
                          colors='black', 
                          linewidths=0.8, 
                          levels=levels_2
                          )
    
    ax.clabel(contorno, 
              inline = 1, 
              inline_spacing = 1, 
              fontsize=20, 
              fmt = '%3.0f', 
              colors= 'black'
              )
    
    # adiciona legenda 
    barra_de_cores = plt.colorbar(sombreado, 
                                  orientation = 'horizontal', 
                                  pad=0.04, 
                                  fraction=0.04
                                  )
    font_size = 20 # Adjust as appropriate.
    barra_de_cores.ax.tick_params(labelsize=font_size)
    
    # Add a title
    plt.title('Vorticidade relativa (1/s) - 1000 hPa',
              fontweight='bold', 
              fontsize=35, 
              loc='left'
              )
    #previsao
    #plt.title('Valid Time: {}'.format(vtime), fontsize=35, loc='right')
    #analise
    plt.title('Análise: {}'.format(vtime), fontsize=35, loc='right')
    
    
    #--------------------------------------------------------------------------
    # Salva imagem
    plt.savefig(f'/home/coqueiro/ufrj/Estagio_supervisionado/imagens/vorticidade/vorticidade_relativa_{format(vtime)}.png', bbox_inches='tight')