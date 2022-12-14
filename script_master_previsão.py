#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 23:34:26 2022

@author: everson
"""

###################### Importando bibliotecas ############################
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
import os

############################ Dataset ###################################
file_1 = xr.open_dataset(
    '/home/everson/Downloads/GFS_Global_0p25deg_20221123_1800.grib2.nc4'
    ).metpy.parse_cf()

# Converte longitude
file_1 = file_1.assign_coords(dict(
    longitude = (((file_1.longitude.values + 180) % 360) - 180))
    ).sortby('longitude')

# Extensão (Extent)
lon_slice = slice(-120., -20.)
lat_slice = slice(15., -65.)
lon_slice_2 = slice(-120., -20.,15)
lat_slice_2 = slice(15., -65.,15)

# Seleciona as latitudes e longitudes
lats = file_1.latitude.sel(latitude=lat_slice).values
lons = file_1.longitude.sel(longitude=lon_slice).values
lats_2 = file_1.latitude.sel(latitude=lat_slice_2).values
lons_2 = file_1.longitude.sel(longitude=lon_slice_2).values

#################### Seleciona os niveis em Z #########################
level_1 = 1000 * units('hPa')
level_2 = 500 * units('hPa')
level_3 = 850 * units('hPa')
level_4 = 250 * units('hPa') 
level_5 = 200 * units('hPa')

#################### Atribui DX e DY ###################################
dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)

######################### Escalas de cores ##############################

# Linhas geopotencial em 500 hPa
colors = ["#0000CD","#FF0000"]
cmap_1 = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)
cmap_1.set_over("#FF0000")
cmap_1.set_under("#0000CD")

# Espessura
colors = ["#2d001c", "#5b0351", "#780777", "#480a5e", "#1e1552", 
          "#1f337d", "#214c9f", "#2776c6", "#2fa5f1", "#1bad1d", 
          "#8ad900", "#ffec00", "#ffab00", "#f46300", "#de3b00", 
          "#ab1900", "#6b0200", '#3c0000']

cmap_2 = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)
cmap_2.set_over('#3c0000')
cmap_2.set_under('#28000a')

################### Looping do calculos ###############################
for i in range(len(file_1.variables['time2'])):
    
    args_0 = dict(
        time = file_1.time[i], 
        latitude=lat_slice, 
        longitude=lon_slice
        )
    
    args_1 = dict(
        time2 = file_1.time2[i], 
        latitude=lat_slice, 
        longitude=lon_slice
        )
    
    args_2 = dict(
        time2 = file_1.time2[i], 
        vertical=level_2, 
        latitude=lat_slice, 
        longitude=lon_slice
        )
    
    args_3 = dict(
        time2 = file_1.time2[i], 
        vertical = level_3, 
        latitude = lat_slice, 
        longitude = lon_slice
        )
    
    args_4 = dict(
        time2 = file_1.time2[i], 
        vertical=level_3, 
        latitude=lat_slice_2, 
        longitude=lon_slice_2
        )
    
    args_5 = dict(
        time2 = file_1.time2[i], 
        vertical=level_4, 
        latitude=lat_slice, 
        longitude=lon_slice)
    
    args_6 = dict(
        time2 = file_1.time2[i], 
        vertical=level_1,
        latitude=lat_slice, 
        longitude=lon_slice
        )
    
    args_7 = dict(
        time2 = file_1.time2[i], 
        vertical=level_4, 
        latitude=lat_slice, 
        longitude=lon_slice)
    
    # Geopotencial em 1000 hPa (Adv. temp e vort)
    geopotencial_1000 = file_1.Geopotential_height_isobaric.metpy.sel(**args_3).metpy.unit_array.squeeze()
    
    # Geopotencial em 500 hPa (Adv. temp e vort)
    geopotencial_500 = file_1.Geopotential_height_isobaric.metpy.sel(**args_2).metpy.unit_array.squeeze()
    
    # Componente U em 850 hPa (Adv. temp e vort)
    u_500 = file_1['u-component_of_wind_isobaric'].metpy.sel(**args_2).metpy.unit_array.squeeze().to('kt')
    
    # Componente V em 850 hPa (Adv. temp e vort)
    v_500 = file_1['v-component_of_wind_isobaric'].metpy.sel(**args_2).metpy.unit_array.squeeze().to('kt')
    
    # Componente U em 850 hPa (Adv. temp e vort)
    u_850 = file_1['u-component_of_wind_isobaric'].metpy.sel(**args_3).metpy.unit_array.squeeze().to('kt')
    
    # Componente V em 850 hPa (Adv. temp e vort)
    v_850 = file_1['v-component_of_wind_isobaric'].metpy.sel(**args_3).metpy.unit_array.squeeze().to('kt')
    
    # Componete U_2 em 850 hPa (Adv. temp e vort)
    u_850_2 = file_1['u-component_of_wind_isobaric'].metpy.sel(**args_4).metpy.unit_array.squeeze().to('kt')
    
    # Componente V_2 em 850 hPa (Adv. temp e vort)
    v_850_2 = file_1['v-component_of_wind_isobaric'].metpy.sel(**args_4).metpy.unit_array.squeeze().to('kt')
    
    # Componente U em 1000 hPa (Vort rel 1000)
    u_1000 = file_1['u-component_of_wind_isobaric'].metpy.sel(**args_6).metpy.unit_array.squeeze()
    
    # Componente V em 1000 hPa (Vort rel 1000)
    v_1000 = file_1['v-component_of_wind_isobaric'].metpy.sel(**args_6).metpy.unit_array.squeeze()
    
    # Componente U em 250 hPa (Corrente de jato 250)
    u_250 = file_1['u-component_of_wind_isobaric'].metpy.sel(**args_5).metpy.unit_array.squeeze()
    
    # Componente V em 250 hPa (Corrente de jato 250)
    v_250 = file_1['v-component_of_wind_isobaric'].metpy.sel(**args_5).metpy.unit_array.squeeze()
    
    # Componente U em 250 hPa (Corrente de jato 250)
    u_200 = file_1['u-component_of_wind_isobaric'].metpy.sel(**args_7).metpy.unit_array.squeeze()
    
    # Componente V em 250 hPa (Corrente de jato 250)
    v_200 = file_1['v-component_of_wind_isobaric'].metpy.sel(**args_7).metpy.unit_array.squeeze()
    
    # Pressão ao nivel medio do mar
    pnmm = file_1.Pressure_reduced_to_MSL_msl.metpy.sel(**args_1).metpy.unit_array.squeeze()* 0.01 * units.hPa/units.Pa
    
    # Temperatura em 850 hPa
    temp = file_1.Temperature_isobaric.metpy.sel(**args_3).metpy.unit_array.squeeze().to('degC')
    
    # ROLE
    rol = file_1['Upward_Long-Wave_Radp_Flux_atmosphere_top_Mixed_intervals_Average'].metpy.sel(**args_0).metpy.unit_array.squeeze()
    
    # Umidade especifica
    q = file_1.Specific_humidity_isobaric.metpy.sel(**args_3).metpy.unit_array.squeeze()* 1e3
    
    ############################### DATA ######################################
    vtempo = file_1.time2.data[i].astype('datetime64[ms]').astype('O')
    
    ######################## Advecçõa de Temperatura ##########################
    
    ###################### Especificações do Plot #############################
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
    
    # Calculo Adv. de temperatura
    adv_temp = mpcalc.advection(temp, u = u_850, v = v_850, dx = dx, dy = dy, x_dim = -1, y_dim = -2) * 1000
    
    # Intevalos da adv de temp
    intervalo_min_1 = -2.5
    intervalo_max_1 = 2.6
    interval_1 = 0.1              # de quanto em quanto voce quer que varie
    levels_1 = np.arange(intervalo_min_1, intervalo_max_1, interval_1)
    
    # intevalos da geopotencial
    intervalo_min_2 = np.amin(np.array(geopotencial_500))
    intervalo_max_2 = np.amax(np.array(geopotencial_500))
    interval_2 = 50              # de quanto em quanto voce quer que varie
    levels_2 = np.arange(intervalo_min_2, intervalo_max_2, interval_2)

    # Plota a imagem adv de temp
    adv_temp_sombreado = ax.contourf(lons, 
                            lats, 
                            adv_temp, 
                            cmap='seismic', 
                            levels = levels_1, 
                            extend = 'neither'
                            )
    
    # Plota a imagem geopotencial
    geopotencial_500_contorno = ax.contour(lons,
                          lats, 
                          geopotencial_500,
                          cmap = cmap_1,
                          linewidths=3, 
                          linestyles='dashed',
                          levels=levels_2
                          )
    
    ax.clabel(geopotencial_500_contorno, 
              inline = 1, 
              inline_spacing = 1, 
              fontsize=24,
              fmt = '%3.0f', 
              colors= 'black'
              )
    
    plt.barbs(
      lons_2,
      lats_2,
      u_850_2,
      v_850_2,
      fill_empty=True,
      length=7,
      sizes=dict(emptybarb=0.1, height=0.8),
      barbcolor="black",
      barb_increments=dict(flag=50),)
    
    #adicionando shapefile
    shapefile = list(
        shpreader.Reader(
        '/home/everson/Downloads/GFS-analysis_and_forecast-main/shapefiles/BR_UF_2021/BR_UF_2021.shp'
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
    barra_de_cores = plt.colorbar(adv_temp_sombreado, 
                                  orientation = 'horizontal', 
                                  pad=0.04, 
                                  fraction=0.04
                                  )
    font_size = 20 # Adjust as appropriate.
    barra_de_cores.ax.tick_params(labelsize=font_size)
    
   
    # Add a title
    plt.title('Adveccao de temperatura (K/s) - 850 hPa',
              fontweight='bold', 
              fontsize=30, 
              loc='left'
              )
    
    # Titulo da previsao
    plt.title('Valid time: {}'.format(vtempo), fontsize=20, loc='right')
    
    #--------------------------------------------------------------------------
    # Salva imagem
    plt.savefig(f'/home/everson/Estagio_supervisionado/imagens/adv_temperatura/adveccao_de_temperatura_{vtempo}.png', bbox_inches='tight')
    plt.clf()
    ######################################################################################################################################
    ######################################################################################################################################

    ######################### Advecção de vorticidade relativa negativa ##############################
    
    ###################### Especificações do Plot #############################
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
    
    # Calculo Adv. vorticidade relativa negativa
    vorticidade = mpcalc.vorticity(u_1000, v_1000, dx = dx, dy = dy, x_dim = -1, y_dim = -2)*10**5
    adv_vort = mpcalc.advection(vorticidade, u = u_500, v = v_500, dx = dx, dy = dy, x_dim = -1, y_dim = -2)*100
    
    # Intevalos da advecção de vorticidade
    intervalo_min_1 = -1.6
    intervalo_max_1 = 0
    interval_1 = 0.10
    levels_1 = np.arange(intervalo_min_1, intervalo_max_1, interval_1)
    
    # intevalos da geopotencial
    intervalo_min_2 = np.amin(np.array(geopotencial_500))
    intervalo_max_2 = np.amax(np.array(geopotencial_500))
    interval_2 = 50              # de quanto em quanto voce quer que varie
    levels_2 = np.arange(intervalo_min_2, intervalo_max_2, interval_2)

    # plota a imagem adv vorticidade
    adv_vort_sombreado = ax.contourf(lons, 
                            lats, 
                            adv_vort, 
                            cmap='gist_heat', 
                            levels = levels_1, 
                            extend = 'min'
                            )

    # plota a imagem geopotencial
    adv_vort_contorno = ax.contour(lons,
                            lats, 
                            geopotencial_500, 
                            cmap=cmap_1, 
                            linewidths=3, 
                            levels=levels_2
                            )
    
    ax.clabel(adv_vort_contorno, 
              inline = 1, 
              inline_spacing = 1, 
              fontsize=24,
              fmt = '%3.0f', 
              colors= 'black'
              )
    
    plt.barbs(
      lons_2,
      lats_2,
      u_850_2,
      v_850_2,
      fill_empty=True,
      length=7,
      sizes=dict(emptybarb=0.1, height=0.8),
      barbcolor="black",
      barb_increments=dict(flag=50),)
    
    #adicionando shapefile
    shapefile = list(
        shpreader.Reader(
        '/home/everson/Downloads/GFS-analysis_and_forecast-main/shapefiles/BR_UF_2021/BR_UF_2021.shp'
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
    barra_de_cores = plt.colorbar(adv_vort_sombreado, 
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
    
    # Titulo da previsão
    plt.title('Valid time: {}'.format(vtempo), fontsize=20, loc='right')
   
    
    #--------------------------------------------------------------------------
    # Salva imagem
    plt.savefig(f'/home/everson/Estagio_supervisionado/imagens/adv_vorticidade/adveccao_de_vorticidade_{vtempo}.png', bbox_inches='tight')
    plt.clf()
    
    ######################################################################################################################################
    ######################################################################################################################################
    
    ################################### Corrente de jato ###################################
    
    ###################### Especificações do Plot #############################
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
    
    # Calculo Corrente de jato 250 hPa
    mag = np.sqrt(u_250**2+v_250**2)

    # intevalos da corrente de jato
    intervalo_min_1 = 30
    intervalo_max_1 = 110
    interval_1 = 10            # de quanto em quanto voce quer que varie
    levels_1 = np.arange(intervalo_min_1, intervalo_max_1, interval_1)
    
    # corrente de jato
    corrente_jato_sobreado = plt.contourf(lons,
                       lats, 
                       mag, 
                       cmap=cmocean.cm.amp, 
                       levels = levels_1, 
                       extend='both')
    
    ax.streamplot(lons, 
                  lats, 
                  u_250, 
                  v_250, 
                  density=[4,4], 
                  linewidth=2, 
                  arrowsize=2.5,
                  color='black', 
                  transform=ccrs.PlateCarree())
    
    
    #adicionando shapefile
    shapefile = list(
        shpreader.Reader(
        '/home/everson/Downloads/GFS-analysis_and_forecast-main/shapefiles/BR_UF_2021/BR_UF_2021.shp'
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
    barra_de_cores = plt.colorbar(corrente_jato_sobreado, 
                                  orientation = 'horizontal', 
                                  pad=0.04, 
                                  fraction=0.04
                                  )
    font_size = 20 # Adjust as appropriate.
    barra_de_cores.ax.tick_params(labelsize=font_size)
    
    # Add a title
    plt.title('Corrente de jato (m/s) - 250 hPa',
              fontweight='bold', 
              fontsize=30, 
              loc='left'
              )
    #previsao
    plt.title('Valid time: {}'.format(vtempo), fontsize=20, loc='right')
    
    #--------------------------------------------------------------------------
    # Salva imagem
    plt.savefig(f'/home/everson/Estagio_supervisionado/imagens/corrente_de_jato/corrente_de_jato_250_{format(vtempo)}.png', bbox_inches='tight')
    plt.clf()
    
    ###################### Especificações do Plot #############################
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
    
    # Calculo Corrente de jato 250 hPa
    mag = np.sqrt(u_200**2+v_200**2)

    # intevalos da corrente de jato
    intervalo_min_1 = 30
    intervalo_max_1 = 110
    interval_1 = 10            # de quanto em quanto voce quer que varie
    levels_1 = np.arange(intervalo_min_1, intervalo_max_1, interval_1)
    
    # corrente de jato
    corrente_jato_sobreado = plt.contourf(lons,
                       lats, 
                       mag, 
                       cmap=cmocean.cm.amp, 
                       levels = levels_1, 
                       extend='both')
    
    ax.streamplot(lons, 
                  lats, 
                  u_250, 
                  v_250, 
                  density=[4,4], 
                  linewidth=2, 
                  arrowsize=2.5,
                  color='black', 
                  transform=ccrs.PlateCarree())
    
    
    #adicionando shapefile
    shapefile = list(
        shpreader.Reader(
        '/home/everson/Downloads/GFS-analysis_and_forecast-main/shapefiles/BR_UF_2021/BR_UF_2021.shp'
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
    barra_de_cores = plt.colorbar(corrente_jato_sobreado, 
                                  orientation = 'horizontal', 
                                  pad=0.04, 
                                  fraction=0.04
                                  )
    font_size = 20 # Adjust as appropriate.
    barra_de_cores.ax.tick_params(labelsize=font_size)
    
    # Add a title
    plt.title('Corrente de jato (m/s) - 200 hPa',
              fontweight='bold', 
              fontsize=30, 
              loc='left'
              )
    #previsao
    plt.title('Valid time: {}'.format(vtempo), fontsize=20, loc='right')
    
    #--------------------------------------------------------------------------
    # Salva imagem
    plt.savefig(f'/home/everson/Estagio_supervisionado/imagens/corrente_de_jato/corrente_de_jato_200_{format(vtempo)}.png', bbox_inches='tight')
    plt.clf()
    
    ######################################################################################################################################
    ######################################################################################################################################
    
    #################################### Divergencia ################################################
    ###################### Especificações do Plot #############################
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
    
    # Calculo da divergencia
    divergencia = mpcalc.divergence(u_1000, v_1000, dx = dx, dy = dy, x_dim=- 1, y_dim=- 2) *  1e6
    
    # intevalos da divergencia - umidade
    divq_min = -50
    divq_max = -20
    interval_1 = 2              # de quanto em quanto voce quer que varie
    levels_1 = np.arange(divq_min, divq_max, interval_1)
    
    # plota a imagem divergencia
    div_sombreado = ax.contourf(lons, 
                            lats, 
                            divergencia, 
                            cmap = cmocean.cm.algae, 
                            levels = levels_1, 
                            extend = 'min'
                            )

    ax.streamplot(lons, lats, u_1000, v_1000, density=[4,4], linewidth=2, arrowsize=3, color='black', transform=ccrs.PlateCarree())

    #adicionando shapefile
    shapefile = list(
        shpreader.Reader(
        '/home/everson/Downloads/GFS-analysis_and_forecast-main/shapefiles/BR_UF_2021/BR_UF_2021.shp'
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
    barra_de_cores = plt.colorbar(div_sombreado, 
                                  orientation = 'horizontal', 
                                  pad=0.04, 
                                  fraction=0.04
                                  )
    font_size = 20 # Adjust as appropriate.
    barra_de_cores.ax.tick_params(labelsize=font_size)
    

    
    # Add a title
    plt.title('Convergência em 1000hPa',
              fontweight='bold', 
              fontsize=30, 
              loc='left'
              )
    
    #previsao
    plt.title('Valid time: {}'.format(vtempo), fontsize=20, loc='right')
    
    plt.savefig(f'/home/everson/Estagio_supervisionado/imagens/divergencia/divergencia_{vtempo}.png', bbox_inches='tight')
    plt.clf()
    
    ######################################################################################################################################
    ######################################################################################################################################
    
    #################################### Espessura ################################################
    ###################### Especificações do Plot #############################
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
    
    # Calculo da espessura
    espessura = geopotencial_500 - geopotencial_1000
    
    # intevalos da espessura
    intervalo_min_1 = 4900
    intervalo_max_1 = 5880
    interval_1 = 25             # de quanto em quanto voce quer que varie
    levels_1 = np.arange(intervalo_min_1, intervalo_max_1, interval_1)
    
    # intevalos da pnmm
    intervalo_min_2 = np.amin(np.array(pnmm))
    intervalo_max_2 = np.amax(np.array(pnmm))
    interval_2 = 2              # de quanto em quanto voce quer que varie
    levels_2 = np.arange(intervalo_min_2, intervalo_max_2, interval_2)
    
    # plota a imagem espessura
    espessura_sombreado = ax.contourf(lons, 
                            lats, 
                            espessura, 
                            cmap=cmap_2, 
                            levels = levels_1, 
                            extend = 'neither'
                            )
    
    # plota a imagem pressao
    espessura_contorno = ax.contour(lons,
                          lats, 
                          pnmm, 
                          colors='black', 
                          linewidths=2, 
                          levels=levels_2
                          )
    
    ax.clabel(espessura_contorno, 
              inline = 1, 
              inline_spacing = 1, 
              fontsize=20, 
              fmt = '%3.0f', 
              colors= 'black'
              )
    
    
    #adicionando shapefile
    shapefile = list(
        shpreader.Reader(
        '/home/everson/Downloads/GFS-analysis_and_forecast-main/shapefiles/BR_UF_2021/BR_UF_2021.shp'
        ).geometries()
        )
    
    ax.add_geometries(
        shapefile, ccrs.PlateCarree(), 
        edgecolor = 'black', 
        facecolor='none', 
        linewidth=0.5
        )
    
    # adiciona mascara de terra
    ax.add_feature(cfeature.LAND)
    # adiciona continente e bordas
    ax.coastlines(resolution='10m', color='black', linewidth=3)
    ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=3)
    
    ### script para colocar zoom numa região
    # # inset axes....
    # axins = ax.inset_axes([0.58, 0.58, 0.4, 0.4])
    # axins.contourf(sombreado, cmap=cmap, origin="image")
    # contorno_2 = axins.contour(contorno_1, colors='black', origin="image")
    # axins.clabel(contorno_2, 
    #           inline = 1, 
    #           inline_spacing = 1, 
    #           fontsize=15, 
    #           fmt = '%3.0f', 
    #           colors= 'black'
    #           )
    
    # # sub region of the original image
    # x1, x2, y1, y2 = -55, -40, -35 , -15
    # axins.set_xlim(x1, x2)
    # axins.set_ylim(y1, y2)
    # axins.set_xticklabels([])
    # axins.set_yticklabels([])
    # ax.indicate_inset_zoom(axins, edgecolor="black")
    
    # adiciona legenda 
    barra_de_cores = plt.colorbar(espessura_sombreado,
                                  orientation = 'horizontal', 
                                  pad=0.04, 
                                  fraction=0.04
                                  )
    font_size = 20 # Adjust as appropriate.
    barra_de_cores.ax.tick_params(labelsize=font_size)
    
    	
    # Add a title
    plt.title('Espessura (500-1000)hPa',
              fontweight='bold', 
              fontsize=30, 
              loc='left'
              )
    
    #previsao
    plt.title('Valid time: {}'.format(vtempo), fontsize=20, loc='right')
    
    #--------------------------------------------------------------------------
    # Salva imagem
    plt.savefig(f'/home/everson/Estagio_supervisionado/imagens/espessura/espessura_{vtempo}.png', bbox_inches='tight')
    plt.clf()
    
    ######################################################################################################################################
    ######################################################################################################################################
    
    #################################### Vorticidade ################################################
    ###################### Especificações do Plot #############################
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
    
    # Calculo da vorticidade
    vorticidade = mpcalc.vorticity(u_1000, v_1000, dx=dx, dy=dy, x_dim=- 1, y_dim=- 2)*10**5
    
    # intevalos da pnmm
    intervalo_min_1 = np.amin(np.array(pnmm))
    intervalo_max_1 = np.amax(np.array(pnmm))
    interval_1 = 2              # de quanto em quanto voce quer que varie
    levels_1 = np.arange(intervalo_min_1, intervalo_max_1, interval_1)
    
    # intevalos da vorticidade
    intervalo_min_2 = -20
    intervalo_max_2 = -2
    interval_2 = 2             # de quanto em quanto voce quer que varie
    levels_2 = np.arange(intervalo_min_2, intervalo_max_2, interval_2)
    
    # adiciona mascara de terra
    ax.add_feature(cfeature.LAND)
    
    # plota a imagem da vorticidade
    vort_sombreado = ax.contourf(lons, 
                            lats, 
                            vorticidade, 
                            cmap=cmocean.cm.dense_r, 
                            levels = levels_1, 
                            extend = 'min'
                            )
    
    # plota a imagem pressao
    vort_contorno = ax.contour(lons,
                          lats, 
                          pnmm, 
                          colors='black', 
                          linewidths=2, 
                          levels=levels_2
                          )
    
    ax.clabel(vort_contorno, 
              inline = 1, 
              inline_spacing = 1, 
              fontsize=20, 
              fmt = '%3.0f', 
              colors= 'black'
              )
    
    #adicionando shapefile
    shapefile = list(
        shpreader.Reader(
        '/home/everson/Downloads/GFS-analysis_and_forecast-main/shapefiles/BR_UF_2021/BR_UF_2021.shp'
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
    barra_de_cores = plt.colorbar(vort_sombreado, 
                                  orientation = 'horizontal', 
                                  pad=0.04, 
                                  fraction=0.04
                                  )
    font_size = 20 # Adjust as appropriate.
    barra_de_cores.ax.tick_params(labelsize=font_size)
    
    # Add a title
    plt.title('Vorticidade relativa (1/s) - 1000 hPa',
              fontweight='bold', 
              fontsize=30, 
              loc='left'
              )
    #previsao
    plt.title('Valid time: {}'.format(vtempo), fontsize=20, loc='right')
    
    #--------------------------------------------------------------------------
    # Salva imagem
    plt.savefig(f'/home/everson/Estagio_supervisionado/imagens/vorticidade/vorticidade_relativa_{format(vtempo)}.png', bbox_inches='tight')
    plt.clf()