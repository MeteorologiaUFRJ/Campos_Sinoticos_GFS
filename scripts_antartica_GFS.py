#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 28 19:47:34 2022

@author: everson
"""


import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean
import metpy.calc as mpcalc
from metpy.units import units
import cartopy.io.shapereader as shpreader # Import shapefiles
import matplotlib.colors as mcolors


file_1 = xr.open_mfdataset('/home/everson/Trabalho_antartica/dados/analise/*grib2.nc4').metpy.parse_cf()

file_1 = file_1.assign_coords(dict(
    longitude = (((file_1.longitude.values + 180) % 360) - 180))
    ).sortby('longitude')

# file_2 = xr.open_dataset(
#     '/home/everson/Trabalho_antartica/dados/single _var_CE.nc').metpy.parse_cf()

# file_2 = file_2.assign_coords(dict(
#     longitude = (((file_2.longitude.values + 180) % 360) - 180))
#     ).sortby('longitude')

# Extensão (Extent)
lat_slice = slice(-45.,-90.)
lon_slice = slice(-180.,180.)
lon_slice_2 = slice(-180., 180.,15)
lat_slice_2 = slice(-45., -90.,15)

#pega as lat/lon
lats = file_1.latitude.sel(latitude=lat_slice).values
lons = file_1.longitude.sel(longitude=lon_slice).values
lats_2 = file_1.latitude.sel(latitude=lat_slice_2).values
lons_2 = file_1.longitude.sel(longitude=lon_slice_2).values

#################### Seleciona os niveis em Z #########################
level_1 = 1000 * units('hPa')
level_2 = 500 * units('hPa')
level_3 = 850 * units('hPa')
level_4 = 300 * units('hPa')
level_5 = 200 * units('hPa')

#################### Atribui DX e DY ###################################
dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)

######################### Escalas de cores ##############################

# Linhas geopotencial em 500 hPa
colors = ["#0000CD","#FF0000"]
cmap_1 = mcolors.LinearSegmentedColormap.from_list("", colors)
cmap_1.set_over("#FF0000")
cmap_1.set_under("#0000CD")

# Espessura
colors = ["#2d001c", "#5b0351", "#780777", "#480a5e", "#1e1552", 
          "#1f337d", "#214c9f", "#2776c6", "#2fa5f1", "#1bad1d", 
          "#8ad900", "#ffec00", "#ffab00", "#f46300", "#de3b00", 
          "#ab1900", "#6b0200", '#3c0000']

cmap_2 = mcolors.LinearSegmentedColormap.from_list("", colors)
cmap_2.set_over('#3c0000')
cmap_2.set_under('#28000a')

# Umidade especifica
# define os intervalos da legenda
clevs = np.array([0, 2, 4, 6, 8, 10, 12, 14, 16])

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
cmap_3 = mcolors.LinearSegmentedColormap.from_list(
    'specific humidity cmap', colors, clevs.shape[0] - 1)
cmap_3.set_over(np.array([0, 37, 89, 255])/255)
cmap_3.set_under('white')

# nromaliza com base nos intervalos
norm = mcolors.BoundaryNorm(clevs, cmap_3.N) # usado no PColorMesh

# geopotencial_1000 = np.zeros([len(file_1['time']),len(file_1['latitude']),len(file_1['longitude'])])
# geopotencial_500 = np.zeros([len(file_1['time']),len(file_1['latitude']),len(file_1['longitude'])])
# u_500 = np.zeros([len(file_1['time']),len(file_1['latitude']),len(file_1['longitude'])])
# v_500 = np.zeros([len(file_1['time']),len(file_1['latitude']),len(file_1['longitude'])])
# u_850 = np.zeros([len(file_1['time']),len(file_1['latitude']),len(file_1['longitude'])])
# v_850 = np.zeros([len(file_1['time']),len(file_1['latitude']),len(file_1['longitude'])])
# u_1000 = np.zeros([len(file_1['time']),len(file_1['latitude']),len(file_1['longitude'])])
# v_1000 = np.zeros([len(file_1['time']),len(file_1['latitude']),len(file_1['longitude'])])
# u_300 = np.zeros([len(file_1['time']),len(file_1['latitude']),len(file_1['longitude'])])
# v_300 = np.zeros([len(file_1['time']),len(file_1['latitude']),len(file_1['longitude'])])
# pnmm = np.zeros([len(file_1['time']),len(file_1['latitude']),len(file_1['longitude'])])
# temp = np.zeros([len(file_1['time']),len(file_1['latitude']),len(file_1['longitude'])])
# q = np.zeros([len(file_1['time']),len(file_1['latitude']),len(file_1['longitude'])])

################### Looping do calculos ###############################
for i in range(len(file_1.variables['time'])):
    
    args_0 = dict(
        time = file_1.time[i], 
        latitude=lat_slice, 
        longitude=lon_slice
        )
    
    args_1 = dict(
        time = file_1.time[i], 
        latitude=lat_slice, 
        longitude=lon_slice
        )
    
    args_2 = dict(
        time = file_1.time[i], 
        vertical=level_2, 
        latitude=lat_slice, 
        longitude=lon_slice
        )
    
    args_3 = dict(
        time = file_1.time[i], 
        vertical = level_3, 
        latitude = lat_slice, 
        longitude = lon_slice
        )
    
    args_4 = dict(
        time = file_1.time[i], 
        vertical=level_3, 
        latitude=lat_slice_2, 
        longitude=lon_slice_2
        )
    
    args_5 = dict(
        time = file_1.time[i], 
        vertical=level_4, 
        latitude=lat_slice, 
        longitude=lon_slice)
    
    args_6 = dict(
        time = file_1.time[i], 
        vertical=level_1,
        latitude=lat_slice, 
        longitude=lon_slice
        )
    
    args_7 = dict(
        time = file_1.time[i], 
        vertical=level_4, 
        latitude=lat_slice, 
        longitude=lon_slice)
    
    # Geopotencial em 1000 hPa (Adv. temp e vort)
    
    geopotencial_1000 = file_1.Geopotential_height_isobaric.metpy.sel(**args_6).metpy.unit_array.squeeze()
    
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
    
    
    # Componente U em 1000 hPa (Vort rel 1000)
    u_1000 = file_1['u-component_of_wind_isobaric'].metpy.sel(**args_6).metpy.unit_array.squeeze()
    
    # Componente V em 1000 hPa (Vort rel 1000)
    v_1000 = file_1['v-component_of_wind_isobaric'].metpy.sel(**args_6).metpy.unit_array.squeeze()
    
    # Componente U em 250 hPa (Corrente de jato 250)
    u_300 = file_1['u-component_of_wind_isobaric'].metpy.sel(**args_5).metpy.unit_array.squeeze()
    
    # Componente V em 250 hPa (Corrente de jato 250)
    v_300 = file_1['v-component_of_wind_isobaric'].metpy.sel(**args_5).metpy.unit_array.squeeze()
    
    # Pressão ao nivel medio do mar
    pnmm = file_1.Pressure_reduced_to_MSL_msl.metpy.sel(**args_1).metpy.unit_array.squeeze()* 0.01 * units.hPa/units.Pa
    
    # Temperatura em 850 hPa
    temp = file_1.Temperature_isobaric.metpy.sel(**args_3).metpy.unit_array.squeeze().to('degC')
    
    # Umidade especifica em 850hPa
    q = file_1.Specific_humidity_isobaric .metpy.sel(**args_3).metpy.unit_array.squeeze()* 1e3
    
    # # Temperatura a t2m
    # temp_t2m = file_2.t2m.metpy.sel(**args_1).metpy.unit_array.squeeze().to('degC')
    
    # # Temperatura da superficie do mar
    # temp_sst = file_2.sst.metpy.sel(**args_1).metpy.unit_array.squeeze().to('degC')
    
    # # Fluxo médio de calor latente
    # flux_latent = file_2.mslhf.metpy.sel(**args_1).metpy.unit_array.squeeze()
    
    # # Fluxo médio de calor sensivel
    # flux_sensible = file_2.msshf.metpy.sel(**args_1).metpy.unit_array.squeeze()
    
    ############################### DATA ######################################
    vtempo = file_1.time.data[i].astype('datetime64[ms]').astype('O')
    
    ######################## Advecçõa de Temperatura ##########################
    
    ###################### Especificações do Plot #############################
    # escolha o tamanho do plot em polegadas (largura x altura)
    plt.figure(figsize=(25,25))
    ax = plt.axes(projection=ccrs.SouthPolarStereo())
    gl = ax.gridlines(crs=ccrs.PlateCarree(),
                      color='gray',
                      alpha=1.0, 
                      linestyle='--', 
                      linewidth=0.5,
                      xlocs=np.arange(-180, 180, 10), 
                      ylocs=np.arange(-90, -45, 10), 
                      draw_labels=True)
    gl.xlabel_style = {'size': 29, 'color': 'black', 'rotation': 0}
    gl.ylabel_style = {'size': 29, 'color': 'black', 'rotation': 0}
    gl.top_labels = True
    gl.right_labels = False
    gl.bottom_labels = False
    gl.rotate_labels = False
    ax.set_extent([-90, 30, -90, -45], ccrs.PlateCarree())
    
    
    ax.add_feature(cfeature.LAND)
    ax.coastlines(resolution='10m', color='black', linewidth=2)
    ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=2)
    
    
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
                            extend = 'neither',
                            transform=ccrs.PlateCarree()
                            )
    
    # Plota a imagem geopotencial
    geopotencial_500_contorno = ax.contour(lons,
                          lats, 
                          geopotencial_500,
                          colors= 'black',
                          linewidths=1.5, 
                          linestyles='dashed',
                          levels=levels_2,
                          transform=ccrs.PlateCarree()
                          )
    
    ax.clabel(geopotencial_500_contorno, 
              inline = 1, 
              inline_spacing = 1, 
              fontsize=24,
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
    
    plt.title('Antártica - GFS 0.25\n'
              'Adveccao de temperatura (K/s) (shaded) - 850 hPa\n'
              'Valid time: {}'.format(vtempo),
              fontweight='bold', 
              fontsize=35, 
              loc='left'
              )
    
    # Titulo da previsao
    #plt.title('Valid time: {}'.format(vtempo), fontsize=20, loc='right')
    
    #--------------------------------------------------------------------------
    # Salva imagem
    plt.savefig(f'/home/everson/Estagio_supervisionado/imagens/adv_temperatura/adveccao_de_temperatura_{vtempo}.png', bbox_inches='tight')
    plt.close()
    ######################################################################################################################################
    ######################################################################################################################################

    ######################### Advecção de vorticidade relativa negativa ##############################
    
    ###################### Especificações do Plot #############################
    # escolha o tamanho do plot em polegadas (largura x altura)
    plt.figure(figsize=(25,25))
    ax = plt.axes(projection=ccrs.SouthPolarStereo())
    gl = ax.gridlines(crs=ccrs.PlateCarree(),
                      color='gray',
                      alpha=1.0, 
                      linestyle='--', 
                      linewidth=0.5,
                      xlocs=np.arange(-180, 180, 10), 
                      ylocs=np.arange(-90, -45, 10), 
                      draw_labels=True)
    
    gl.xlabel_style = {'size': 29, 'color': 'black', 'rotation': 0}
    gl.ylabel_style = {'size': 29, 'color': 'black', 'rotation': 0}
    gl.top_labels = True
    gl.right_labels = False
    gl.bottom_labels = False
    gl.rotate_labels = False
    ax.set_extent([-90, 30, -90, -45], ccrs.PlateCarree())
    
    # Calculo Adv. vorticidade relativa negativa
    
    vorticidade = mpcalc.vorticity(u_1000, v_1000, dx = dx, dy = dy, x_dim = -1, y_dim = -2)*10**5
    adv_vort = mpcalc.advection(vorticidade, u = u_500, v = v_500, dx = dx, dy = dy, x_dim = -1, y_dim = -2)*100
    
    # Intevalos da advecção de vorticidade
    intervalo_min_3 = -1.6
    intervalo_max_3 = 0
    interval_3 = 0.10
    levels_3 = np.arange(intervalo_min_3, intervalo_max_3, interval_3)
    
    # intevalos da geopotencial
    intervalo_min_4 = np.amin(np.array(geopotencial_500))
    intervalo_max_4 = np.amax(np.array(geopotencial_500))
    interval_4 = 50              # de quanto em quanto voce quer que varie
    levels_4 = np.arange(intervalo_min_4, intervalo_max_4, interval_4)

    # plota a imagem adv vorticidade
    adv_vort_sombreado = ax.contourf(lons, 
                            lats, 
                            adv_vort, 
                            cmap='gist_heat', 
                            levels = levels_3, 
                            extend = 'min',
                            transform=ccrs.PlateCarree()
                            )

    # plota a imagem geopotencial
    adv_vort_contorno = ax.contour(lons,
                            lats, 
                            geopotencial_500, 
                            colors= 'black', 
                            linewidths=1.5, 
                            levels=levels_4,
                            transform=ccrs.PlateCarree()
                            )
    
    ax.clabel(adv_vort_contorno, 
              inline = 1, 
              inline_spacing = 1, 
              fontsize=24,
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
    
    	
    plt.title('Antártica - GFS 0.25\n'
              'Adv de vort relativa (1/s²) (shaded) - 500 hPa\n'
              'Valid time: {}'.format(vtempo),
              fontweight='bold', 
              fontsize=35, 
              loc='left'
              )
    
    # Titulo da previsão
    #plt.title('Valid time: {}'.format(vtempo), fontsize=20, loc='right')
   
    
    #--------------------------------------------------------------------------
    # Salva imagem
    plt.savefig(f'/home/everson/Estagio_supervisionado/imagens/adv_vorticidade/adveccao_de_vorticidade_{vtempo}.png', bbox_inches='tight')
    plt.close()
    
    ######################################################################################################################################
    ######################################################################################################################################
   
    ################################### Corrente de jato ###################################
    
    ###################### Especificações do Plot #############################
    # escolha o tamanho do plot em polegadas (largura x altura)
    plt.figure(figsize=(25,25))
    ax = plt.axes(projection=ccrs.SouthPolarStereo())
    gl = ax.gridlines(crs=ccrs.PlateCarree(),
                      color='gray',
                      alpha=1.0, 
                      linestyle='--', 
                      linewidth=0.5,
                      xlocs=np.arange(-180, 180, 10), 
                      ylocs=np.arange(-90, -45, 10), 
                      draw_labels=True)
    
    gl.xlabel_style = {'size': 29, 'color': 'black', 'rotation': 0}
    gl.ylabel_style = {'size': 29, 'color': 'black', 'rotation': 0}
    gl.top_labels = True
    gl.right_labels = False
    gl.bottom_labels = False
    gl.rotate_labels = False
    ax.set_extent([-90, 30, -90, -45], ccrs.PlateCarree())
    
    
    # Calculo Corrente de jato 250 hPa
    
    mag_1 = np.sqrt(u_300**2+v_300**2)

    # intevalos da corrente de jato
    intervalo_min_5 = 30
    intervalo_max_5 = 110
    interval_5 = 10            # de quanto em quanto voce quer que varie
    levels_5 = np.arange(intervalo_min_5, intervalo_max_5, interval_5)
    
    # corrente de jato
    corrente_jato_sobreado = plt.contourf(lons,
                       lats, 
                       mag_1, 
                       cmap=cmocean.cm.amp, 
                       levels = levels_5, 
                       extend='both',
                       transform=ccrs.PlateCarree())
    
    ax.streamplot(np.array(lons), 
                  np.array(lats), 
                  np.array(u_300), 
                  np.array(v_300), 
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
    plt.title('Antártica - GFS 0.25\n'
              'Corrente de jato (m/s) - 300 hPa (shaded) - 500 hPa\n'
              'Valid time: {}'.format(vtempo),
              fontweight='bold', 
              fontsize=35, 
              loc='left'
              )
    
    
    #previsao
    #plt.title('Valid time: {}'.format(vtempo), fontsize=20, loc='right')
    
    #--------------------------------------------------------------------------
    # Salva imagem
    plt.savefig(f'/home/everson/Estagio_supervisionado/imagens/corrente_de_jato/corrente_de_jato_300_{format(vtempo)}.png', bbox_inches='tight')
    plt.close()
    
    #################################### Espessura ################################################
    ###################### Especificações do Plot #############################
    # escolha o tamanho do plot em polegadas (largura x altura)
    plt.figure(figsize=(25,25))
    ax = plt.axes(projection=ccrs.SouthPolarStereo())
    gl = ax.gridlines(crs=ccrs.PlateCarree(),
                      color='gray',
                      alpha=1.0, 
                      linestyle='--', 
                      linewidth=0.5,
                      xlocs=np.arange(-180, 180, 10), 
                      ylocs=np.arange(-90, -45, 10), 
                      draw_labels=True
                      )
    gl.xlabel_style = {'size': 29, 'color': 'black', 'rotation': 0}
    gl.ylabel_style = {'size': 29, 'color': 'black', 'rotation': 0}
    gl.top_labels = True
    gl.right_labels = False
    gl.bottom_labels = False
    gl.rotate_labels = False
    ax.set_extent([-90, 30, -90, -45], ccrs.PlateCarree())
    
    # Calculo da espessura
    espessura = geopotencial_500 - geopotencial_1000
    
    # intevalos da espessura
    intervalo_min_8 = 4900
    intervalo_max_8 = 5880
    interval_8 = 25             # de quanto em quanto voce quer que varie
    levels_8 = np.arange(intervalo_min_8, intervalo_max_8, interval_8)
    
    # intevalos da pnmm
    intervalo_min_9 = np.amin(np.array(pnmm))
    intervalo_max_9 = np.amax(np.array(pnmm))
    interval_9 = 2              # de quanto em quanto voce quer que varie
    levels_9 = np.arange(intervalo_min_9, intervalo_max_9, interval_9)
    
    # plota a imagem espessura
    espessura_sombreado = ax.contourf(lons, 
                            lats, 
                            espessura, 
                            cmap=cmap_2, 
                            levels = levels_8, 
                            extend = 'neither',
                            transform=ccrs.PlateCarree()
                            )
    
    # plota a imagem pressao
    espessura_contorno = ax.contour(lons,
                          lats, 
                          pnmm, 
                          colors='black', 
                          linewidths=2, 
                          levels=levels_9,
                          transform=ccrs.PlateCarree()
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
    
    plt.title('Antártica - GFS 0.25\n'
              'Espessura (shaded) - 500 hPa\n'
              'Valid time: {}'.format(vtempo),
              fontweight='bold', 
              fontsize=35, 
              loc='left'
              )
    
    
    
    #previsao
    #plt.title('Valid time: {}'.format(vtempo), fontsize=20, loc='right')
    
    #--------------------------------------------------------------------------
    # Salva imagem
    plt.savefig(f'/home/everson/Estagio_supervisionado/imagens/espessura/espessura_{vtempo}.png', bbox_inches='tight')
    plt.close()
    
    ######################################################################################################################################
    ######################################################################################################################################
    
    #################################### Vorticidade ################################################
    ###################### Especificações do Plot #############################
    # escolha o tamanho do plot em polegadas (largura x altura)
    plt.figure(figsize=(25,25))
    ax = plt.axes(projection=ccrs.SouthPolarStereo())
    gl = ax.gridlines(crs=ccrs.PlateCarree(),
                      color='gray',
                      alpha=1.0, 
                      linestyle='--', 
                      linewidth=0.5,
                      xlocs=np.arange(-180, 180, 10), 
                      ylocs=np.arange(-90, -45, 10), 
                      draw_labels=True
                      )
    gl.xlabel_style = {'size': 29, 'color': 'black', 'rotation': 0}
    gl.ylabel_style = {'size': 29, 'color': 'black', 'rotation': 0}
    gl.top_labels = True
    gl.right_labels = False
    gl.bottom_labels = False
    gl.rotate_labels = False
    ax.set_extent([-90, 30, -90, -45], ccrs.PlateCarree())
    
    # Calculo da vorticidade
    vorticidade_1000 = mpcalc.vorticity(u_1000, v_1000, dx = dx, dy = dy, x_dim = -1, y_dim = -2)*10**5
    
    # intevalos da pnmm
    intervalo_min_10 = np.amin(np.array(pnmm))
    intervalo_max_10 = np.amax(np.array(pnmm))
    interval_10 = 2              # de quanto em quanto voce quer que varie
    levels_10 = np.arange(intervalo_min_10, intervalo_max_10, interval_10)
    
    # intevalos da vorticidade
    intervalo_min_11 = -20
    intervalo_max_11 = -2
    interval_11 = 2             # de quanto em quanto voce quer que varie
    levels_11 = np.arange(intervalo_min_11, intervalo_max_11, interval_11)
    
    # adiciona mascara de terra
    ax.add_feature(cfeature.LAND)
    
    # plota a imagem da vorticidade
    vort_sombreado = ax.contourf(lons, 
                            lats, 
                            vorticidade_1000, 
                            cmap=cmocean.cm.dense_r, 
                            levels = levels_11, 
                            extend = 'min',
                            transform=ccrs.PlateCarree()
                            )
    
    # plota a imagem pressao
    pressao_contorno = ax.contour(lons,
                          lats, 
                          pnmm, 
                          colors='black', 
                          linewidths=2, 
                          levels=levels_10,
                          transform=ccrs.PlateCarree()
                          )
    
    ax.clabel(pressao_contorno, 
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
    
    plt.title('Antártica - GFS 0.25\n'
              'Vorticidade relativa (1/s) (shaded) - 1000 hPa\n'
              'Valid time: {}'.format(vtempo),
              fontweight='bold', 
              fontsize=35, 
              loc='left'
              )
    
    
    #previsao
    #plt.title('Valid time: {}'.format(vtempo), fontsize=20, loc='right')
    
    #--------------------------------------------------------------------------
    # Salva imagem
    plt.savefig(f'/home/everson/Estagio_supervisionado/imagens/vorticidade/vorticidade_relativa_{format(vtempo)}.png', bbox_inches='tight')
    plt.close()

    ######################################################################################################################################
    ######################################################################################################################################
    
    
    #################################### Umidade especifica ################################################
    ###################### Especificações do Plot #############################
    # escolha o tamanho do plot em polegadas (largura x altura)
    plt.figure(figsize=(25,25))
    ax = plt.axes(projection=ccrs.SouthPolarStereo())
    gl = ax.gridlines(crs=ccrs.PlateCarree(),
                      color='gray',
                      alpha=1.0, 
                      linestyle='--', 
                      linewidth=0.5,
                      xlocs=np.arange(-180, 180, 10), 
                      ylocs=np.arange(-90, -45, 10), 
                      draw_labels=True
                      )
    gl.xlabel_style = {'size': 29, 'color': 'black', 'rotation': 0}
    gl.ylabel_style = {'size': 29, 'color': 'black', 'rotation': 0}
    gl.top_labels = True
    gl.right_labels = False
    gl.bottom_labels = False
    gl.rotate_labels = False
    ax.set_extent([-90, 30, -90, -45], ccrs.PlateCarree())
    
    # umidade especifica 
    umidade_sombreado = ax.contourf(lons, lats, q, cmap = cmap_3, levels = clevs, extend='both', transform=ccrs.PlateCarree())
    ax.streamplot(np.array(lons), np.array(lats), np.array(u_850), np.array(v_850), density=[4,4], linewidth=2, arrowsize=3, color='black', transform=ccrs.PlateCarree())
    
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
    barra_de_cores = plt.colorbar(umidade_sombreado, 
                                  orientation = 'horizontal', 
                                  pad=0.04, 
                                  fraction=0.04
                                  )
    font_size = 20 # Adjust as appropriate.
    barra_de_cores.ax.tick_params(labelsize=font_size)
    
    	
    # Add a title
    plt.title('Antártica - GFS 0.25\n'
              'Umidade específica 850 hPa\n'
              'Valid time: {}'.format(vtempo),
              fontweight='bold', 
              fontsize=35, 
              loc='left'
              )
    
    #previsao
    ##plt.title('Valid time: {}'.format(vtempo), fontsize=20, loc='right')
    
    #analise
    plt.title('Análise: {}'.format(vtempo), fontsize=35, loc='right')
    
    #--------------------------------------------------------------------------
    # Salva imagem
    plt.savefig(f'/home/everson/Estagio_supervisionado/imagens/umidade_especifica/umidade_especifica_850hPa_{format(vtempo)}.png', bbox_inches='tight')
    
