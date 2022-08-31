# ============================================================================ #
# Imports 
# tambem é necessário cfrgib
# ============================================================================ #

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pygrib
import cartopy.crs as ccrs
import cartopy.feature as cfeature 
import cartopy.io.shapereader as shpreader
import matplotlib.colors as mcolors
from datetime import datetime

# ============================================================================ #
# abrindo arquivos
# ============================================================================ #

# arquivo GFS
fname = r"F:\Lucas\Conteudo\estagio_2\Felipe\Dados\gfs.0p25.2021120712.f000.grib2"
args = dict(engine = 'cfgrib', filter_by_keys={'typeOfLevel': 'isobaricInhPa'})
ds = xr.open_dataset(fname, **args)

# para printar a descrição de todas as variaveis
# print(ds.variables)

# ============================================================================ #
# Variáveis e propriedades
# ============================================================================ #

# Caminho até o shapefile com as fronteiras de estados ou paises ou municipios
fname_shape = r"F:\Lucas\Conteudo\Fisica das nuvens e precipitacao\Dados\BR_UF_2019.shp"

# qualidade da figura
dpi = 100
figsize = (20, 20)

# projeção do mapa
proj = ccrs.PlateCarree()

# limites de lat e lon (Int, Int) NAO USE FLOATS
lat_min = -70
lat_max = 10
lon_min = -105
lon_max = -20

# seleciona uma extensao minima para o recorte de dados
extent = [lon_min, lon_max, lat_min, lat_max]
lon_slice = slice(lon_min, lon_max) # menor pro maior
lat_slice = slice(lat_max, lat_min) # maior pro menor

# nível da atmosfera
level = 850

# convertendo a longitude do grib de [0, 360] para [-180, 180]
transform = dict(longitude = (((ds.longitude.values + 180) % 360) - 180))
ds = ds.assign_coords(transform).sortby('longitude')

# Slicing
ds = ds.sel(isobaricInhPa = level, latitude = lat_slice, longitude = lon_slice)
lons = ds.longitude.values
lats = ds.latitude.values

# selecionando variaveis
u_comp = ds['u'].values
v_comp = ds['v'].values
q = ds['q'].values * 1e3 # esta em kg * kg, converte para g/kg 

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

# ============================================================================ #
# Plot do mapa
# ============================================================================ #

# escolha o tamanho do plot em polegadas (largura x altura) * dpi
fig = plt.figure(figsize=figsize, dpi = dpi)

# usando a projeção da coordenada cilindrica equidistante
ax = plt.axes(projection=proj)
ax.set_extent(extent, proj)

# adiciona continente e bordas
shapefile_geometries = list(shpreader.Reader(fname_shape).geometries())
ax.add_geometries(shapefile_geometries, proj,
                    edgecolor='black', facecolor='none', linewidth=0.5)

ax.coastlines(resolution='50m', color='black', linewidth=1.5)
ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=1)
gl = ax.gridlines(
    crs=proj,
    color='white',
    alpha=1.0,
    linestyle='--',
    linewidth=0.25,
    xlocs=np.arange(-180, 180, 5),
    ylocs=np.arange(-90, 90, 5),
    draw_labels=True
    )

gl.top_labels = False
gl.right_labels = False

# ============================================================================ #
# Plot dos campos meteorologicos
# ============================================================================ #

# corrente de jato
img = ax.contourf(lons, lats, q, cmap = cmap, levels = clevs, extend='both')
# img2 = ax.contour(lons, lats, q, colors='white', linewidths=0.3,  transform=proj)
ax.streamplot(lons, lats, u_comp, v_comp, density=[6,6], linewidth=1.5, color='black', transform=ccrs.PlateCarree())

# adiciona legenda 
cb = fig.colorbar(img, extend ='both', orientation = 'horizontal', pad=0.04, fraction=0.04)
font_size = 20 # Adjust as appropriate.
cb.ax.tick_params(labelsize=font_size)

# ============================================================================ #
# Propriedades da figura
# ============================================================================ #

# Getting the file time and date
datetime64 = ds.valid_time.values
unix_epoch = np.datetime64(0, 's')
one_second = np.timedelta64(1, 's')
seconds_since_epoch = (datetime64 - unix_epoch) / one_second
date = datetime.utcfromtimestamp(seconds_since_epoch)
date_formatted = date.strftime('%Y-%m-%d %HZ')
	
# Add a title
plt.title('Umidade específica (g/kg) & Linhas de corrente - 850 hPa', fontweight='bold', fontsize=20, loc='left')
plt.title(f'{date_formatted}', fontsize=15, loc='right')

#----------------------------------------------------------------------------------------------------------- 

# Salva imagem
plt.savefig(f'umidade_corrente_850 {date_formatted}.png', bbox_inches='tight')