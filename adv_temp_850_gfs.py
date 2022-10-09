#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 17:39:26 2022

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

#dataset
file_1 = xr.open_dataset(
    '/home/coqueiro/Downloads/prev.nc4'
    ).metpy.parse_cf()

file_1 = file_1.assign_coords(dict(
    longitude = (((file_1.longitude.values + 180) % 360) - 180))
    ).sortby('longitude')

#
#extent
lon_slice = slice(-90., 10.)
lat_slice = slice(15., -70.)
lon_slice_2 = slice(-90., 10.,15)
lat_slice_2 = slice(15., -70.,15)

#pega as lat/lon
lats = file_1.latitude.sel(latitude=lat_slice).values
lons = file_1.longitude.sel(longitude=lon_slice).values

lats_2 = file_1.latitude.sel(latitude=lat_slice_2).values
lons_2 = file_1.longitude.sel(longitude=lon_slice_2).values

#seta as variaveis
level_1 = 850 * units('hPa')
level_2 = 500 * units('hPa')
level_3 = 1000 * units('hPa')

# cria uma escala de cores:
colors = ["#0000CD","#FF0000"]
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)
cmap.set_over("#FF0000")
cmap.set_under("#0000CD")


for i in range(len(file_1.variables['time3'])):
    
    geopotencial = file_1.Geopotential_height_isobaric.metpy.sel(
        time3 = file_1.time3[i], 
        vertical=level_2, 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze()
    
    geopotencial_1000 = file_1.Geopotential_height_isobaric.metpy.sel(
        time3 = file_1.time3[i], 
        vertical=level_3, 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze()
    
    geopotencial_500 = file_1.Geopotential_height_isobaric.metpy.sel(
        time3 = file_1.time3[i], 
        vertical=level_2, 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze()
    
    u = file_1['u-component_of_wind_isobaric'].metpy.sel(
        time3 = file_1.time3[i], 
        vertical=level_1, 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze()
    
    v = file_1['v-component_of_wind_isobaric'].metpy.sel(
        time3 = file_1.time3[i], 
        vertical=level_1, 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze()
    
    u_2 = file_1['u-component_of_wind_isobaric'].metpy.sel(
        time3 = file_1.time3[i], 
        vertical=level_1, 
        latitude=lat_slice_2, 
        longitude=lon_slice_2
        ).metpy.unit_array.squeeze().to('kt')
    
    v_2 = file_1['v-component_of_wind_isobaric'].metpy.sel(
        time3 = file_1.time3[i], 
        vertical=level_1, 
        latitude=lat_slice_2, 
        longitude=lon_slice_2
        ).metpy.unit_array.squeeze().to('kt')
    
    pnmm = file_1.Pressure_reduced_to_MSL_msl.metpy.sel(
        time3 = file_1.time3[i], 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze()* 0.01 * units.hPa/units.Pa
    
    t = file_1.Temperature_isobaric.metpy.sel(
        time3 = file_1.time3[i],  
        vertical=level_1, 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze()
    
    #data
    vtime3 = file_1.time3.data[i].astype('datetime64[ms]').astype('O')
    
    dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)
    
    adv_temp = mpcalc.advection(t, u=u, v=v, dx=dx, dy=dy, x_dim=- 1, y_dim=- 2)*1000
    
    espessura = geopotencial_500 - geopotencial_1000
    
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
    
    
    # intevalos da pnmm
    intervalo_min2 = np.amin(np.array(pnmm))
    intervalo_max2 = np.amax(np.array(pnmm))
    interval_2 = 2              # de quanto em quanto voce quer que varie
    levels_2 = np.arange(intervalo_min2, intervalo_max2, interval_2)
    
    # intevalos da adv de temp
    intervalo_min3 = -1.5
    intervalo_max3 = 1.6
    interval_3 = 0.1              # de quanto em quanto voce quer que varie
    levels_3 = np.arange(intervalo_min3, intervalo_max3, interval_3)
    
    # intevalos da geopotencial
    intervalo_min1 = np.amin(np.array(geopotencial_500))
    intervalo_max1 = np.amax(np.array(geopotencial_500))
    interval_1 = 50              # de quanto em quanto voce quer que varie
    levels_1 = np.arange(intervalo_min1, intervalo_max1, interval_1)

    
    # plota a imagem adv de temp
    sombreado = ax.contourf(lons, 
                            lats, 
                            adv_temp, 
                            cmap='seismic', 
                            levels = levels_3, 
                            extend = 'neither'
                            )
    
    # #plota a imagem pressao
    # contorno_1 = ax.contour(lons,
    #                       lats, 
    #                       pnmm, 
    #                       colors='black', 
    #                       linewidths=1, 
    #                       levels=levels_2
    #                       )
    
    # ax.clabel(contorno_1, 
    #           inline = 1, 
    #           inline_spacing = 1, 
    #           fontsize=20, 
    #           fmt = '%3.0f', 
    #           colors= 'black'
    #           )
    
    # plota a imagem geopotencial
    contorno_2 = ax.contour(lons,
                          lats, 
                          geopotencial_500,
                          cmap=cmap, 
                          linewidths=3, 
                          linestyles='dashed',
                          levels=levels_1
                          )
    
    ax.clabel(contorno_2, 
              inline = 1, 
              inline_spacing = 1, 
              fontsize=24,
              fmt = '%3.0f', 
              colors= 'black'
              )
    
    plt.barbs(
      lons_2,
      lats_2,
      u_2,
      v_2,
      fill_empty=True,
      length=7,
      sizes=dict(emptybarb=0.1, height=0.6),
      barbcolor="black",
      barb_increments=dict(flag=50),)
    
    #adicionando shapefile
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
    ax.add_feature(cfeature.LAND)
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
    plt.title('Adveccao de temperatura (K/s) - 850 hPa',
              fontweight='bold', 
              fontsize=35, 
              loc='left'
              )
    
    #previsao
    plt.title('Valid Time: {}'.format(vtime3), fontsize=25, loc='right')
    #analise
    #plt.title('Análise: {}'.format(vtime), fontsize=35, loc='right')
    
    #--------------------------------------------------------------------------
    # Salva imagem
    plt.savefig(f'/home/coqueiro/ufrj/Estagio_supervisionado/imagens/adv_temperatura/adveccao_de_temperatura_{vtime3}.png', bbox_inches='tight')
