#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  8 23:17:33 2022

@author: coqueiro
"""

#importando bibliotecas
import cartopy.crs as ccrs #biblioteca de cartografia
import cartopy.feature as cfeature #biblioteca de que coloca shapefile
import matplotlib.pyplot as plt #biblioteca de plots de mapas
import matplotlib.colors #biblioteca de paleta de cores
import metpy.calc as mpcalc #biblioteca de funções físicas meteorológicas
from metpy.units import units #import de unidades físicas
import numpy as np #biblioteca que utiliza arrays 
import xarray as xr #biblioteca que ler dados netcdf
import cartopy.io.shapereader as shpreader # Import shapefiles
from datetime import datetime, timedelta  # basicas datas e tipos de tempo
import cmocean #biblioteca de cores

# input de dados (ERA5, GFS, CFRS)
file_1 = xr.open_dataset(
    '/home/ladsin/Downloads/GFS_analise_11_13.nc4'
    ).metpy.parse_cf()

# conversão de longitude de graus (0°-360°) para radianos (-180 a 180)
file_1 = file_1.assign_coords(dict(
    longitude = (((file_1.longitude.values + 180) % 360) - 180))
    ).sortby('longitude')


# extent - seleciona as latitudes norte (lat_0), sul (lat_1), leste (lon_1), oeste (lon_0)
lon_0 = -120.
lon_1 = -20.
lat_0 = 10.
lat_1 = -55.

# faz uma fatia com as latitude e longitudes para selecionar no dataset
lon_slice = slice(lon_0, lon_1)
lat_slice = slice(lat_0, lat_1)

# seleciona as fatias nas lats e lons do dataset
lats = file_1.latitude.sel(latitude=lat_slice).values
lons = file_1.longitude.sel(longitude=lon_slice).values

# seleciona os niveis verticais em hpa, porém tem diferença para alguns tipos de inputs, PORTANTO COMENTE DE ACORDO COM O SEU DADO!!

## GFS
level_2 = 500 * units('hPa')
level_1 = 250 * units('hPa')

## ERA5
#level_2 = 500 
#level_1 = 250 

# cria uma escala de cores:
colors = ["#0000CD","#FF0000"]
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)
cmap.set_over("#FF0000")
cmap.set_under("#0000CD")


for i in range(len(file_1.variables['time'])):

#ANTES DE SELECIONAR A VARIVAVEL, VEJA QUAIS SÃO OS NOMES REFERENTES AS ELAS COM O COMANDO "file_1" no terminal!!
    
    # VARIAVEL DO GEOPOTENCIAL
    geopotencial = file_1.Geopotential_height_isobaric.metpy.sel(
        time = file_1.time[i], 
        vertical=level_2, 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze()
    
    # COMPONENTE DO VENTO U
    u = file_1['u-component_of_wind_isobaric'].metpy.sel(
        time = file_1.time[i], 
        vertical=level_1, 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze()

    # COMPONENTE DO VENTO V
    v = file_1['v-component_of_wind_isobaric'].metpy.sel(
        time = file_1.time[i], 
        vertical=level_1, 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze()

    # PRESSÃO A NIVEL MÉDIO DO MAR
    pnmm = file_1.Pressure_reduced_to_MSL_msl.metpy.sel(
        time = file_1.time[i], 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze()* 0.01 * units.hPa/units.Pa

    # CALCULA A MAGNITUDE DO VENTO 
    mag = np.sqrt(u**2+v**2)
    
    # SELECIONA A DATA DE ACORDO COM O DADO
    vtime = file_1.time.data[i].astype('datetime64[ms]').astype('O')

    
    # CONFIGURAÇÕES DA GRADE DO PLOT
    # ESCOLHA O TAMANHO DO PLOT EM POLEGADAS (LARGURA X ALTURA)
    plt.figure(figsize=(25,25))
    
    # PROJEÇÃO DA COORDENADA CILINDRICA EQUIDISTANTE 
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
    
    
    # INTERVALOS DA ESCALA DE PRESSÃO
    intervalo_min2 = np.amin(np.array(pnmm))
    intervalo_max2 = np.amax(np.array(pnmm))
    interval_2 = 2              # de quanto em quanto voce quer que varie
    levels_2 = np.arange(intervalo_min2, intervalo_max2, interval_2)
    
    # INTERVALOS DA ESCALA DA CORRENTE DE JATO
    intervalo_min3 = 30
    intervalo_max3 = 110
    interval_3 = 10              # de quanto em quanto voce quer que varie
    levels_3 = np.arange(intervalo_min3, intervalo_max3, interval_3)
    
    # INTERVALOS DA ESCALA DE GEOPOTENCIAL
    intervalo_min1 = np.amin(np.array(geopotencial))
    intervalo_max1 = np.amax(np.array(geopotencial))
    interval_1 = 50              # de quanto em quanto voce quer que varie
    levels_1 = np.arange(intervalo_min1, intervalo_max1, interval_1)

    # SOMBREADO = PLT.CONTOURF()
    # LINHAS = PLT.CONTOUR()
    
    # PLOTANDO A CORRENTE DE JATO
    jato = plt.contourf(lons,
                       lats, 
                       mag, 
                       cmap=cmocean.cm.deep, 
                       levels = levels_3, 
                       extend='neither')

    # PLOTANDO ISOBARAS DE PRESSÃO
    pressao = ax.contour(lons,
                          lats, 
                          pnmm, 
                          colors='black', 
                          linewidths=2, 
                          levels=levels_2
                          )
    # INTRODUZ OS VALORES DAS ISOBARAS
    ax.clabel(pressao, 
              inline = 1, 
              inline_spacing = 1, 
              fontsize=20, 
              fmt = '%3.0f', 
              colors= 'black'
              )
    
    # PLOTANDO O GEOPOTENCIAL
    geopot = ax.contour(lons,
                          lats, 
                          geopotencial,
                          cmap=cmap, 
                          linewidths=3, 
                          linestyles='dashed',
                          levels=levels_1
                          )
    
    # INTRODUZ OS VALORES DOS GEOPOTENCIAIS, SE FOR NECESSÁRIO, PORÉM FICA BEM POLUIDO!
    # ax.clabel(geopot, 
    #           inline = 1, 
    #           inline_spacing = 1, 
    #           fontsize=24,
    #           fmt = '%3.0f', 
    #           colors= 'black'
    #           )

    # ADICIONA O SHAPEFILE DO DA AMÉRICA DO SUL, É NECESSÁRIO BAIXAR NO SEU COMPUTADOR!!
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
    
    # ADICIONA OS CONTINENTES E AS BORDAS
    ax.add_feature(cfeature.LAND)
    ax.coastlines(resolution='10m', color='black', linewidth=3)
    ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=3)
    
    # ADICIONA A LEGENDA
    barra_de_cores = plt.colorbar(jato, 
                                  orientation = 'horizontal', 
                                  pad=0.04, 
                                  fraction=0.04
                                  )
    font_size = 20 # Adjust as appropriate.
    barra_de_cores.ax.tick_params(labelsize=font_size)
    
   
    # ADICIONA O TITULO
    plt.title('Pnmm (black), Geop.500 (red/blue), Jato 250 (shaded)',
              fontweight='bold', 
              fontsize=28, 
              loc='left'
              )
    
    # TITULO PARA PREVISÃO
    plt.title('Valid Time: {}'.format(vtime), fontsize=20, loc='right')
    
    #analise
    #plt.title('Análise: {}'.format(vtime), fontsize=35, loc='right')
    
    #--------------------------------------------------------------------------
    # Salva imagem
    plt.savefig(f'/work/archive/Everson/Coqueiro/Estagio/plots/tripolar/tripolar_250_{vtime}.png', bbox_inches='tight')
