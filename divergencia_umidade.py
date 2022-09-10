#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  4 11:10:10 2022

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
    '/home/coqueiro/Downloads/GFS_Global_0p25deg_20220910_0600.grib2.nc4'
    ).metpy.parse_cf()

file_1 = file_1.assign_coords(dict(
    longitude = (((file_1.longitude.values + 180) % 360) - 180))
    ).sortby('longitude')

#extent
lon_slice = slice(-120., 10.)
lat_slice = slice(20., -20.)

#pega as lat/lon
lats = file_1.latitude.sel(latitude=lat_slice).values
lons = file_1.longitude.sel(longitude=lon_slice).values

#seta as variaveis
level = 1000 * units('hPa')

# CRIANDO O CMAP E NORM PARA A COLORBAR

# intevalos da divergencia - umidade
divq_min = -2.
divq_max = 2.
n_levs = 10 # numero de intervalos
divlevs = np.round(np.linspace(divq_min, divq_max, n_levs), 1)

# lista de cores, em ordem crescete. RGBA
colors = ['mediumseagreen', 'mediumaquamarine', 'palegreen', 'white', 'wheat', 'gold', 'goldenrod']

# cria um novo cmap a partir do pre-existente
cmap = mcolors.LinearSegmentedColormap.from_list(
    'Custom cmap', colors, divlevs.shape[0] - 1)
cmap.set_over('darkgoldenrod')
cmap.set_under("seagreen")

# nromaliza com base nos intervalos
norm = mcolors.BoundaryNorm(divlevs, cmap.N) # util para o PCOLORMESH, CONTOURF nao usa

# variaveis repetidas em cada loop
dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)


for i in range(len(file_1.variables['time1'])):
    args = dict(
        time = file_1.time1[i] ,
        vertical=level,
        latitude=lat_slice,
        longitude=lon_slice
    )
    
    u = file_1['u-component_of_wind_isobaric'].metpy.sel(**args).metpy.unit_array.squeeze()
    v = file_1['v-component_of_wind_isobaric'].metpy.sel(**args).metpy.unit_array.squeeze()
    q = file_1.Specific_humidity_isobaric.metpy.sel(**args).metpy.unit_array.squeeze()
    
    #data
    vtime1 = file_1.time1.data[i].astype('datetime64[ms]').astype('O')

    divergencia = mpcalc.divergence(u, v, dx=dx, dy=dy, x_dim=- 1, y_dim=- 2)
    divergencia_umidade = divergencia * q * 1e6
    
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
    
    # plota a imagem divergencia
    sombreado = ax.contourf(lons, 
                            lats, 
                            divergencia_umidade, 
                            cmap = cmap, 
                            levels = divlevs, 
                            extend = 'neither'
                            )
    

    ax.streamplot(lons, lats, u, v, density=[3,3], linewidth=1.5, color='black', transform=ccrs.PlateCarree())
    
    
    shapefile = list(
        shpreader.Reader(
        '/home/coqueiro/Downloads/br_unidades_da_federacao/BR_UF_2019.shp'
        ).geometries()
        )
    
    ax.add_geometries(
        shapefile, ccrs.PlateCarree(), 
        edgecolor = 'black', 
        facecolor='none', 
        linewidth=0.5
        )
    
    # adiciona continente e bordas
    ax.coastlines(resolution='10m', color='black', linewidth=1)
    ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=1)
    # adiciona mascara de terra
    ax.add_feature(cfeature.LAND)
    
    
    # adiciona legenda 
    barra_de_cores = plt.colorbar(sombreado, 
                                  orientation = 'horizontal', 
                                  pad=0.04, 
                                  fraction=0.04
                                  )
    font_size = 20 # Adjust as appropriate.
    barra_de_cores.ax.tick_params(labelsize=font_size)
    barra_de_cores.ax.set_xticks(divlevs)
    
    # Add a title
    plt.title('Divergencia de umidade 1000hPa',
              fontweight='bold', 
              fontsize=35, 
              loc='left'
              )
    
    #previsao
    #plt.title('Valid time1: {}'.format(vtime1), fontsize=35, loc='right')
    #analise
    plt.title('Análise: {}'.format(vtime1), fontsize=35, loc='right')
    
    #--------------------------------------------------------------------------
    # Salva imagem
    plt.savefig(f'/home/coqueiro/ufrj/Estagio_supervisionado/imagens/divergencia_umidade/divergencia_umidade_{vtime1}.png', bbox_inches='tight')

