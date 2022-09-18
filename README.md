# GFS-analysis_and_forecast

Passos para baixar os dados de análise:
1. Baixe os seus dados de previsão aqui:https://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/catalog.html
2. Escolha a ultima rodada de sua preferênciaEscolha o dia das análises 
3. Em >>> Access <<< escolha a opção NetcdfSubset
![image](https://user-images.githubusercontent.com/91283739/189402418-94e9d495-ffae-4f84-a3fd-ef30f40b3b36.png)
4. Selecione as variaveis em >>> reftime time isobaric latitude longitude <<<
5. Selecione a variavel "Pressure reduced to MSL @ Mean sea level" em >>> reftime time isobaric1 latitude longitude <<<
6. Ajeite as latitude e longitudes de sua preferencia
7. Escolha o periodo de tempo para análise
8. Selecione a opção "netcdf4" em >>> Output format <<<


Passos para baixar os dados de análise:
1. Baixe os seus dados de análise aqui:https://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg_ana/catalog.html
2. Escolha o dia das análises
3. Em >>> Access <<< escolha a opção NetcdfSubset
![image](https://user-images.githubusercontent.com/91283739/189402173-d35dfdf3-7fc4-4e59-be96-0634c9da36ad.png)
4. Selecione as variaveis em >>> reftime time isobaric latitude longitude <<
5. Selecione a variavel "Pressure reduced to MSL @ Mean sea level" em >>> reftime time isobaric1 latitude longitude <<<
6. Ajeite as latitude e longitudes de sua preferencia
7. Escolha o periodo de tempo para análise
8. Selecione a opção "netcdf4" em >>> Output format <<<

Após baixar os dados, coloquem os caminhos deles no "file_1"
É necessário alterar o shapefile
Caso o shapefile não seja plotado na imagem, tente colocar a parte que plota o continente no final do script antes do plt.savefig
