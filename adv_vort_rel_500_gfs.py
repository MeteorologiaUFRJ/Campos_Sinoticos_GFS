#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 18:20:34 2022

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
step = 15

lon_slice = slice(lon_0, lon_1)
lat_slice = slice(lat_0, lat_1)
lon_slice_2 = slice(lon_0, lon_1, step)
lat_slice_2 = slice(lat_0, lat_1, step)

#pega as lat/lon
lats = file_1.latitude.sel(latitude=lat_slice).values
lons = file_1.longitude.sel(longitude=lon_slice).values

lats_2 = file_1.latitude.sel(latitude=lat_slice_2).values
lons_2 = file_1.longitude.sel(longitude=lon_slice_2).values


#seta as variaveis
level = 500 * units('hPa')
level_3 = 1000 * units('hPa')

# cria uma escala de cores:
colors = ["#0000CD","#FF0000"]
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)
cmap.set_over("#FF0000")
cmap.set_under("#0000CD")


for i in range(len(file_1.variables['time'])):
    
    geopotencial = file_1.Geopotential_height_isobaric.metpy.sel(
        time = file_1.time[i], 
        vertical=level, 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze()
    
    geopotencial_1000 = file_1.Geopotential_height_isobaric.metpy.sel(
        time = file_1.time[i], 
        vertical=level_3, 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze()
    
    geopotencial_500 = file_1.Geopotential_height_isobaric.metpy.sel(
        time = file_1.time[i], 
        vertical=level, 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze()
    
    
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
    
    u_2 = file_1['u-component_of_wind_isobaric'].metpy.sel(
        time = file_1.time[i], 
        vertical=level_3, 
        latitude=lat_slice_2, 
        longitude=lon_slice_2
        ).metpy.unit_array.squeeze().to('kt')
    
    v_2 = file_1['v-component_of_wind_isobaric'].metpy.sel(
        time = file_1.time[i], 
        vertical=level_3, 
        latitude=lat_slice_2, 
        longitude=lon_slice_2
        ).metpy.unit_array.squeeze().to('kt')
    
    
    pnmm = file_1.Pressure_reduced_to_MSL_msl.metpy.sel(
        time = file_1.time[i], 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze()* 0.01 * units.hPa/units.Pa
    
    #data
    vtime = file_1.time.data[i].astype('datetime64[ms]').astype('O')
    
    dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)
    vorticidade = mpcalc.vorticity(u, v, dx=dx, dy=dy, x_dim=- 1, y_dim=- 2)*10**5
    adv_vort = mpcalc.advection(vorticidade, u=u, v=v, dx=dx, dy=dy, x_dim=- 1, y_dim=- 2)*100
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
    
    # intevalos da adv vorticidade
    intervalo_min3 = -1.6
    intervalo_max3 = 0
    interval_3 = 0.10              # de quanto em quanto voce quer que varie
    levels_3 = np.arange(intervalo_min3, intervalo_max3, interval_3)
    
    # intevalos da geopotencial
    intervalo_min1 = np.amin(np.array(geopotencial_500))
    intervalo_max1 = np.amax(np.array(geopotencial_500))
    interval_1 = 50              # de quanto em quanto voce quer que varie
    levels_1 = np.arange(intervalo_min1, intervalo_max1, interval_1)
    
    
    # adiciona mascara de terra
    ax.add_feature(cfeature.LAND)
    
    # plota a imagem adv vorticidade
    sombreado = ax.contourf(lons, 
                            lats, 
                            adv_vort, 
                            cmap='gist_heat', 
                            levels = levels_3, 
                            extend = 'min'
                            )
    
    
    # # plota a imagem pressao
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
                          geopotencial, 
                          cmap=cmap, 
                          linewidths=3, 
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
        '/work/archive/Everson/Coqueiro/script_gfs/GFS-analysis_and_forecast-main/shapefiles/BR_UF_2021/BR_UF_2021.shp'
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
    
    # adiciona legenda 
    barra_de_cores = plt.colorbar(sombreado, 
                                  orientation = 'horizontal', 
                                  pad=0.04, 
                                  fraction=0.04
                                  )
    font_size = 20 # Adjust as appropriate.
    barra_de_cores.ax.tick_params(labelsize=font_size)
    
    	
    # Add a title
    plt.title('Adv de vort relativa (1/s²) - 500 hPa',
              fontweight='bold', 
              fontsize=30, 
              loc='left'
              )
    
    #previsao
    #plt.title('Valid Time: {}'.format(vtime), fontsize=20, loc='right')
    #analise
    plt.title('Análise: {}'.format(vtime), fontsize=25, loc='right')
    
    #--------------------------------------------------------------------------
    # Salva imagem
    plt.savefig(f'/work/archive/Everson/Coqueiro/Estagio/plots/adv_vort/adv_vort_{vtime}.png', bbox_inches='tight')