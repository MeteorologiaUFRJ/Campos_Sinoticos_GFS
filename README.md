# 🗺️ Campos Sinóticos de Análise e Previsão do GFS 🌦️

> Campos de análise e previsão do tempo elaborados para o Estágio Supervisionado II / B 
____________________________________________________________________________________________________________________________
## Instalação das bibliotecas necessárias através do terminal 


#### Caso utilize o ambiente Anaconda, instale as bibliotecas necessárias da seguinte forma: (Recomendado)

```
conda install -c conda-forge xarray netCDF4 cartopy matplotlib numpy cmocean metpy
```

#### Do contrário, as bibliotecas podem ser instaladas a partir dos seguintes comandos: 
```
pip install DateTime
```
```
pip install xarray
```
```
pip install netCDF4
```
```
conda install -c conda-forge cartopy
```
> A biblioteca cartopy precisa ser instalada dessa forma, ou não irá instalar/funcionar
```
pip install matplotlib
```
```
pip install numpy
```
```
pip install cmocean
```
```
pip install metpy
```
____________________________________________________________________________________________________________________________

## Passos para baixar os dados de previsão: 

1. Os dados para previsão podem ser acessados através do link: https://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/catalog.html
2. Escolha a rodada de sua preferência
3. Em >>> Access <<< escolha a opção NetcdfSubset para o download dos dados no formato Netcdf (.nc)
![image](https://user-images.githubusercontent.com/91283739/189402418-94e9d495-ffae-4f84-a3fd-ef30f40b3b36.png)
4. Selecione as variaveis em >>> reftime time isobaric latitude longitude <<<
5. Selecione a variavel "Pressure reduced to MSL @ Mean sea level" em >>> reftime time isobaric1 latitude longitude <<<
6. Selecione as latitudes e longitudes de sua preferência em Horizontal subset
7. Selecione o periodo de tempo para análise em Time subset (Stride: 4 - intervalos de 12h em 12h e Stride: 2 - intervalos de 6h em 6h)
8. Selecione a opção "netcdf4_classic" em Output format
____________________________________________________________________________________________________________________________

## Passos para baixar os dados de análise: 

1. Os dados para previsão podem ser acessados através do link:https://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg_ana/catalog.html
2. Escolha a rodada de sua preferência
3. Em >>> Access <<< escolha a opção NetcdfSubset para o download dos dados no formato Netcdf (.nc)
![image](https://user-images.githubusercontent.com/91283739/189402173-d35dfdf3-7fc4-4e59-be96-0634c9da36ad.png)
4. Selecione as variaveis em >>> reftime time isobaric latitude longitude <<<
5. Selecione a variavel "Pressure reduced to MSL @ Mean sea level" em >>> reftime time isobaric1 latitude longitude <<<
6. Selecione as latitudes e longitudes de sua preferência em Horizontal subset
7. Selecione o periodo de tempo para análise em Time subset (Stride: 4 - intervalos de 12h em 12h e Stride: 2 - intervalos de 6h em 6h)
8. Selecione a opção "netcdf4_classic" em Output format

### É possível baixar mais de um dado de análise utilizando "Full Collection Dataset" no início da página.
1. Siga o tutorial acima
2. Ao chegar no passo 7, o periodo de tempo para a análise muda (Stride: 1 - 6h em 6h; Stride: 2 - 12h em 12h) aumentando gradativamente.

____________________________________________________________________________________________________________________________

## Passos para utilizar o script: 
Após o download do dado ajuste o caminho até o diretório onde foi armazenado no seu computador e substitua em: 

```
"file_1 = xr.open_dataset('/diretorio-do-dado/GFS_Global_0p25deg_20220910_0600.grib2.nc4'
    ).metpy.parse_cf()"
```    
    
Além disso, altere o caminho até o diretório de armazenamento do shapefile (arquivo "BR_UF_2021.zip") em: 

```
shapefile = list(
        shpreader.Reader(
        '/caminho-do-shapefile/br_unidades_da_federacao/BR_UF_2019.shp'
        ).geometries()
        )
```
Para alterar a área de plot, altere a Lat/Lon nessa linha do script. 

```
 ax.set_extent([-90, -20, -60, 10], crs=ccrs.PlateCarree())
```
