#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 10 11:26:22 2022

@author: coqueiro
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
import cartopy.io.shapereader as shpreader # Import shapefiles
from datetime import datetime, timedelta  # basicas datas e tipos de tempo
import cmocean
import matplotlib.colors as mcolors

#dataset

file_1 = xr.open_dataset(
    '/home/ladsin/Downloads/GFS_analise_11_13.nc4'
    ).metpy.parse_cf()

file_1 = file_1.assign_coords(dict(
    longitude = (((file_1.longitude.values + 180) % 360) - 180))
    ).sortby('longitude')

#extent
lon_0 = -120.
lon_1 = -20.
lat_0 = 10.
lat_1 = -55.

lon_slice = slice(lon_0, lon_1)
lat_slice = slice(lat_0, lat_1)

#pega as lat/lon
lats = file_1.latitude.sel(latitude=lat_slice).values
lons = file_1.longitude.sel(longitude=lon_slice).values

level = 850*units('hPa')
for i in range(len(file_1.variables['time'])):

    # selecionando variaveis
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
    
    q = file_1.Specific_humidity_isobaric.metpy.sel(
        time = file_1.time[i], 
        vertical=level, 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze()* 1e3
    
    #data
    vtime = file_1.time.data[i].astype('datetime64[ms]').astype('O')

    # ============================================================================ #
    # Colorbar
    # ============================================================================ #
    
    # define os intervalos da legenda
    clevs = np.array([2, 4, 6, 8, 10, 12, 14, 16, 18])
    
    # lista de cores, em ordem crescete. RGBA
    colors = np.array([ # (R, G, B, A)
        [242, 98, 0, 255],
        [249, 155, 77, 255],
        [254, 217, 118, 255],
        [255, 247, 188, 255],
        [190, 220, 230, 255],
        [156, 194, 255, 255],
        [59, 118, 255, 255],
        [0, 77, 182, 255]
    ]) / 255 # divide por 255
    
    # cria um novo cmap a partir do pre-existente
    cmap = mcolors.LinearSegmentedColormap.from_list(
        'specific humidity cmap', colors, clevs.shape[0] - 1)
    cmap.set_over(np.array([0, 37, 89, 255])/255)
    cmap.set_under('white')
    
    # nromaliza com base nos intervalos
    norm = mcolors.BoundaryNorm(clevs, cmap.N) # usado no PColorMesh
    
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
    
    # corrente de jato
    sombreado = ax.contourf(lons, lats, q, cmap = cmap, levels = clevs, extend='both')
    ax.streamplot(lons, lats, u, v, density=[4,4], linewidth=2, arrowsize=2.5, color='black', transform=ccrs.PlateCarree())
    
    #adicionando shapefile
    shapefile = list(
        shpreader.Reader(
        '/work/archive/Everson/Coqueiro/script_gfs/GFS-analysis_and_forecast-main/shapefiles/BR_UF_2021/BR_UF_2021.shp'
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
    
    # adiciona legenda 
    barra_de_cores = plt.colorbar(sombreado, 
                                  orientation = 'horizontal', 
                                  pad=0.04, 
                                  fraction=0.04
                                  )
    font_size = 20 # Adjust as appropriate.
    barra_de_cores.ax.tick_params(labelsize=font_size)
    
    	
    # Add a title
    plt.title('Umidade específica 850 hPa',
              fontweight='bold', 
              fontsize=35, 
              loc='left'
              )
    
    #previsao
    #plt.title('Valid time: {}'.format(vtime), fontsize=35, loc='right')
    #analise
    plt.title('Análise: {}'.format(vtime), fontsize=35, loc='right')
    
    #--------------------------------------------------------------------------
    # Salva imagem
    plt.savefig(f'/work/archive/Everson/Coqueiro/Estagio/plots/umidade/umidade_{format(vtime)}.png', bbox_inches='tight')
   