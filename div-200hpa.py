#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 17:47:00 2022

@author: bmiranda
"""
#importando bibliotecas
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.colors 
import metpy.calc as mpcalc
from metpy.units import units
import numpy as np
import xarray as xr
import cartopy.io.shapereader as shpreader 
from datetime import datetime, timedelta 
import cmocean
import matplotlib.colors as mcolors

file_1 = xr.open_dataset(
    '/home/bmiranda/Desktop/ES2/bia-isa/dados/GFS_Global_0p25deg_ana_20221011_1200.grib2.nc4'
    ).metpy.parse_cf()

file_1 = file_1.assign_coords(dict(
    longitude = (((file_1.longitude.values + 180) % 360) - 180))
    ).sortby('longitude')

#extent
lon_slice = slice(-90., -10.)
lat_slice = slice(10., -70.)

#pega as lat/lon
lats = file_1.latitude.sel(latitude=lat_slice).values
lons = file_1.longitude.sel(longitude=lon_slice).values

# variaveis repetidas em cada loop
dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)

level = 200 * units('hPa')

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
    
    divergencia = mpcalc.divergence(u, v, dx=dx, dy=dy, x_dim=- 1, y_dim=- 2) *  1e6
    # mag = np.sqrt(u**2+v**2)
    
    #time
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
                      xlocs=np.arange(-180, 180, 5), 
                      ylocs=np.arange(-90, 90, 5), 
                      draw_labels=True
                      )
    gl.top_labels = False
    gl.right_labels = False
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 29, 'color': 'black'}
    gl.ylabel_style = {'size': 29, 'color': 'black'}
    
    # # intevalos da corrente de jato
    # intervalo_min3 = 30
    # intervalo_max3 = 110
    # interval_3 = 10            # de quanto em quanto voce quer que varie
    # levels_3 = np.arange(intervalo_min3, intervalo_max3, interval_3)
   
    # intevalos da divergencia
    divq_min = 0
    divq_max = 200.
    n_levs = 20 # numero de intervalos
    # divlevs = np.round(np.linspace(divq_min, divq_max, n_levs), 1)
    divlevs = np.round(np.linspace(divq_min, divq_max, n_levs))
    # plot
    
    # corrente de jato
    # jato = plt.contourf(lons,
    #                    lats, 
    #                    mag, 
    #                    cmap=cmocean.cm.amp, 
    #                    levels = levels_3, 
    #                    extend='both')

    sombreado = ax.contourf(lons, 
                        lats, 
                        divergencia, 
                        cmap = 'BuPu', 
                        levels = divlevs, 
                        extend = 'min'
                        )

    ax.streamplot(lons, lats, u, v, density=[3,3], linewidth=1.5, color='black', transform=ccrs.PlateCarree())

    #adicionando shapefile
    shapefile = list(
        shpreader.Reader(
        '/home/bmiranda/Desktop/ES2/bia-isa/shapefiles/BR_UF_2021/BR_UF_2021.shp'
        ).geometries()
        )
    
    ax.add_geometries(
        shapefile, ccrs.PlateCarree(), 
        edgecolor = 'black', 
        facecolor='none', 
        linewidth=0.5
        )
    
    # adiciona continente e bordas
    ax.coastlines(resolution='10m', color='black', linewidth=3)
    ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=3)
    
    #adiciona legenda 
    barra_de_cores = plt.colorbar(sombreado, 
                                  orientation = 'horizontal', 
                                  pad=0.04, 
                                  fraction=0.04
                                  )
    font_size = 20 # Adjust as appropriate.
    barra_de_cores.ax.tick_params(labelsize=font_size)
    

    
    # Add a title
    plt.title('Divergência em 200 hPa',
              fontweight='bold', 
              fontsize=35, 
              loc='left'
              )
    
    #previsao
    #plt.title('Valid time: {}'.format(vtime), fontsize=30, loc='right')
    #analise
    plt.title('Análise: {}'.format(vtime), fontsize=35, loc='right')
    # plt.savefig(f'mov_vert-linhas-corrente-500hpa_{vtime}.png', bbox_inches='tight')