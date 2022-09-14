#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 01:12:18 2022

@author: coqueiro
"""

#importando bibliotecas
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import metpy.calc as mpcalc
from metpy.units import units
import numpy as np
import xarray as xr
import cartopy.io.shapereader as shpreader # Import shapefiles
from datetime import datetime, timedelta  # basicas datas e tipos de tempo
import cmocean

#dataset
file_1 = xr.open_dataset(
    '/home/coqueiro/ufrj/lapt/dados/lapt_10-12092022.grib2.nc4'
    ).metpy.parse_cf()

file_1 = file_1.assign_coords(dict(
    longitude = (((file_1.longitude.values + 180) % 360) - 180))
    ).sortby('longitude')

file_2 = xr.open_dataset(
    '/home/coqueiro/ufrj/lapt/dados/prec_rate.nc4'
    ).metpy.parse_cf()

file_2 = file_2.assign_coords(dict(
    longitude = (((file_2.longitude.values + 180) % 360) - 180))
    ).sortby('longitude')


#extent
lon_slice = slice(-120., 10.)
lat_slice = slice(10., -70.)

#pega as lat/lon
lats = file_1.latitude.sel(latitude=lat_slice).values
lons = file_1.longitude.sel(longitude=lon_slice).values

#seta as variaveis
level = 1000 * units('hPa')

for i in range(len(file_1.variables['time'])):
    
    
    u = file_1['u-component_of_wind_isobaric'].metpy.sel(
        time = file_1.time[i], 
        vertical=level, 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze()
    
    v = file_1['v-component_of_wind_isobaric'].metpy.sel(
        time = file_1.time[i], 
        vertical=level, 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze()
    
    
    p = file_2.Precipitation_rate_surface.metpy.sel(
        time = file_2.time[i], 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze()*1e3
    
    
    q = file_1.Specific_humidity_isobaric.metpy.sel(
        time = file_1.time[i],
        vertical=level,
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze()
    
    rm = mpcalc.mixing_ratio_from_specific_humidity(q)*1e3
    
    #data
    vtime = file_1.time.data[i].astype('datetime64[ms]').astype('O')
    
    # escolha o tamanho do plot em polegadas (largura x altura)
    plt.figure(figsize=(25,25))
    
    # usando a projeção da coordenada cilindrica equidistante 
    ax = plt.axes(projection=ccrs.PlateCarree())
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
    
    # intevalos da agua precipitavel
    intervalo_min2 = 0
    intervalo_max2 = 9
    interval_2 = 0.1           # de quanto em quanto voce quer que varie
    levels_2 = np.arange(intervalo_min2, intervalo_max2, interval_2)
    
    # intevalos da corrente de jato
    intervalo_min3 = np.amin(np.array(rm))
    intervalo_max3 = np.amax(np.array(rm))
    interval_3 = 1          # de quanto em quanto voce quer que varie
    levels_3 = np.arange(intervalo_min3, intervalo_max3, interval_3)
    
    # corrente de jato
    img = plt.contourf(lons,
                       lats, 
                       p,
                       cmap='Blues', 
                       extend='max',
                       levels=levels_2,
                       transform=ccrs.PlateCarree())
    
    img2 = ax.contour(lons, 
                      lats, 
                      rm, 
                      colors='black', 
                      linewidths=1, 
                      levels=levels_3, 
                      transform=ccrs.PlateCarree())
    
    ax.clabel(img2, 
              inline = 1, 
              inline_spacing = 1, 
              fontsize=25, 
              fmt = '%3.0f', 
              colors= 'black'
              )
    
    # ax.streamplot(lons, 
    #               lats, 
    #               u, 
    #               v, 
    #               density=[4,4], 
    #               linewidth=2, 
    #               arrowsize=2.5,
    #               color='black', 
    #               transform=ccrs.PlateCarree())
    
    #adicionando shapefile
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
    ax.coastlines(resolution='10m', color='black', linewidth=3)
    ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=3)
    # adiciona mascara de terra
    ax.add_feature(cfeature.LAND)
    
    
    # adiciona legenda 
    barra_de_cores = plt.colorbar(img, 
                                  orientation = 'horizontal', 
                                  pad=0.04, 
                                  fraction=0.04
                                  )
    font_size = 20 # Adjust as appropriate.
    barra_de_cores.ax.tick_params(labelsize=font_size)
    
    # Add a title
    plt.title('Razão de Mistura e \n'
              'Taxa de precipitação (Kg.s/m²) 1000 hPa',
              fontweight='bold', 
              fontsize=35, 
              loc='left'
              )
    #previsao
    #plt.title('Valid Time: {}'.format(vtime), fontsize=35, loc='right')
    #analise
    plt.title('Análise: {}'.format(vtime), fontsize=20, loc='right')
    
    
    #--------------------------------------------------------------------------
    # Salva imagem
    plt.savefig(f'/home/coqueiro/ufrj/Estagio_supervisionado/imagens/taxa_de_prec/taxa_de_precipitacao_{format(vtime)}.png', bbox_inches='tight')