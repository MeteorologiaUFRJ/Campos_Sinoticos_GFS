#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 23 11:57:20 2023

@author: everson
"""
from datetime import datetime
import metpy.calc as mpcalc
from metpy.units import units
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import pandas as pd
import cartopy.feature as cfeature
import glob 
import os


# # função que lê as pastas em um diretório
# def get_folders(directory):
#     return [f for f in os.listdir(directory) if os.path.isdir(os.path.join(directory, f))]

# # diretório dos dados .nc de componente u e v
# directories = '/media/bjerknes/HD_todo_pod/Everson/Coqueiro/CENPES/DADOS/Raoni/dados_input_vort'

# # lista de pasta para serem selecionadas
# folders = get_folders(directories)
# folders.sort(reverse=True)

# # Converte as componentes u e v em vorticidade relativa 
# for j in range(len(folders)):
#     # path_name = str(directories) + '/' + str(folders[j])
    
#     path_name = os.path.join(directories, folders[j])
    
#     ## Abre os dados na pasta
#     data = xr.open_mfdataset(f'{path_name}/*grib2.nc')
#     data['u'] = data['U_GRD_L100'].copy()
#     data['v'] = data['V_GRD_L100'].copy()
    
#     Data = str(folders[j])
    
#     # data = data.sel(level=850)
    
#     # Vai preenchendo o np.array...
#     vort_850= mpcalc.vorticity(data.u, data.v)
    
    
#     #sys.exit()    
#     # Cria o datarray dentro da estrutura do data 
#     data['vo'] = data['u'].copy()
#     data['vo'].data = vort_850
#     data['vo'] = xr.DataArray(vort_850, dims = ['time','lat','lon'])   
#     data['vo'].attrs['standard_name'] = 'atmosphere_relative_vorticity'
#     data['vo'].attrs['long_name'] = 'Vorticity (relative)'
#     data['vo'].attrs['units'] = 's**-1'
#     data['lat'].attrs['units'] = 'degrees_north'
#     data['lon'].attrs['units'] = 'degrees_east'
    
    
#     # Deixa somente o vort_850 dentro do data
#     data = data[['vo']]
#     data = data.chunk()

#     # Salva no arquivo netcdf
#     fname = f'/home/everson/Programas/track-master/indat/{(Data)}.nc'
#     data.to_netcdf(f'/home/everson/Programas/track-master/indat/{(Data)}.nc')
    
#     os.system("ncap2 -O -s 'time=double(time)' {} {}".format(fname,fname))
#     print ('finished saving')

# # Defina o caminho para a pasta contendo os dados .nc
# caminho_dados = "/home/everson/Programas/track-master/indat"

# # Percorra todos os arquivos .nc na pasta e execute o comando bin/track.linux com o filtro especificado em specfiltT42.in
# for arquivo in os.listdir(caminho_dados):
#     sorted(arquivo)
#     if arquivo.endswith(".nc"):
#         nome_arquivo_sem_extensao = os.path.splitext(arquivo)[0]
#         comando = f"/home/everson/Programas/track-master/bin/track.linux -i {arquivo} -f filt < /home/everson/Programas/track-master/filtro_gfs.in"
#         os.system(comando)
        
#         caminho_origem = os.path.join("/home/everson/Programas/track-master/outdat", "specfil.filt_band000")
#         caminho_destino = os.path.join("/home/everson/Programas/track-master/outdat", nome_arquivo_sem_extensao + ".dat")
#         os.rename(caminho_origem, caminho_destino)
        
#         caminho_origem = caminho_destino
#         caminho_destino = os.path.join(caminho_dados, nome_arquivo_sem_extensao + ".dat")
#         os.rename(caminho_origem, caminho_destino)
        
# # Percorra todos os arquivos .dat na pasta indat e execute o comando ./master com as opções especificadas
# for arquivo in os.listdir(caminho_dados):
#     if arquivo.endswith(".dat"):
#         nome_arquivo_sem_extensao = os.path.splitext(arquivo)[0]
#         comando = f"./master -c={nome_arquivo_sem_extensao} -d=now -e=track.linux -i={arquivo} -f=vo1991 -n=1,1000,1 -o=/media/bjerknes/HD_todo_pod/Everson/Coqueiro/CENPES/DADOS/Raoni/track_raoni/master_raoni_novo -r=RUN_AT_ -s=RUNDATIN.VOR -j=RUN_AT.in"
#         os.system(comando)
         
# Definindo o caminho para a pasta onde estão os arquivos .gz
caminho_dados = "/home/everson/Documentos/ssd_antigo/maq_virtual/CENPES/tracks_rodadados_tr_passo_3/Extratropical/"

# Percorrer cada diretório e descompactar o arquivo tr_trs_neg.gz
for diretorio in os.listdir(caminho_dados):
    print(diretorio)
    if os.path.isdir(os.path.join(caminho_dados, diretorio)):
        caminho_arquivo_gz = os.path.join(caminho_dados, diretorio, "tr_trs_neg.gz")
        if os.path.exists(caminho_arquivo_gz):
            os.system(f"gunzip {caminho_arquivo_gz}")
            
            dia = diretorio[:2]
            mes = diretorio[2:4]
            ano = diretorio[4:8]
            hora = '00'
            data = ano+mes+dia+hora
            
            caminho_arquivo_descompactado = os.path.join(caminho_dados, diretorio, "tr_trs_neg")
            os.system(f"mv {caminho_arquivo_descompactado} /home/everson/Programas/track-master/utils/bin")
            
            # Muda para o diretrotio /utils/bin/
            os.chdir('/home/everson/Programas/track-master/utils/bin/')    
            
            os.system(f"./count tr_trs_neg 0 0 5 4 0 {data} 3")
        
            os.system(f"mv tr_trs_neg.new {diretorio}")
    
            os.system(f"tr2csv {diretorio}")
            
            os.system(f"mv alltr.csv  {diretorio}.csv")
            
            os.system(f"mv {diretorio}.csv /home/everson/Documentos/ssd_antigo/maq_virtual/CENPES/Track_planilhas_csv_passo_4/{diretorio}.csv")
            
            
# I.2 - Carrega a pasta contendo todas as planilhas dos Tracks

datafiles = glob.glob('/home/everson/Documentos/ssd_antigo/maq_virtual/CENPES/Track_planilhas_csv_passo_4/*.csv')
datafiles.sort()

# range = ao numero de meses, poderia usar o range(len(datafiles))
for i in range(len(datafiles)):
    
    
    # II - Adiciona os rótulos nas colunas, pois quando o track gera as planilhas geram planilhas sem rótulos
    
    #seleciona a planilha
    path = np.array(pd.read_csv(datafiles[i], sep=',', usecols=[0, 1, 2, 3, 4]))
    
    #adiciona os rotulos das colunas
    df = pd.DataFrame(path, columns = ['N_do_Ciclone', 'Data', 'Longitude', 'Latitude', 'Vorticidade'])
    
    data = datafiles[i].split('/')[8].split('.')[0]
    
    
    # III - Converte as longitudes para -180 a 180 graus
    
    # converte e recorta latitude/longitude
    df['Longitude'] = (((df['Longitude']+180) %360) - 180)
    
    df = df[(df.Longitude>=-57) & (df.Longitude<=-40) & (df.Latitude<=-34) & (df.Latitude>=-58)].copy()
    
    
    # IV - Remove a ocorrencia de ciclones a partir de um numero pré-definido (exemplo = 4)
    
    # Conta a ocorrência de cada item
    counts = df['N_do_Ciclone'].value_counts()
    
    # Cria uma lista com os ciclones a serem removidos
    itens_para_remover = counts[counts < 10].index.tolist()
    
    # Remove os ciclones da lista do DataFrame
    df = df[~df['N_do_Ciclone'].isin(itens_para_remover)]
    
    
    # V - Salva a planilha dentro da pasta automaticamente a partir do rotulo 'Data'
    
    #salva a planilha recortada
    df.to_csv(f'/home/everson/Documentos/ssd_antigo/maq_virtual/CENPES/Planilhas_recortadas_padrão_5/{(data)}.csv', index=False)    
    # df.to_csv(f'/work/archive/Everson/Coqueiro/CENPES/DADOS/track/rodada_1/.csv', index=False) 
    
    
    # plota os tracks
    # Caso queira plotar as saídas das planilhas descomente essa parte!
    
    #lat e lon extend da imagem
    latmin, latmax = -50, -20
    lonmin, lonmax = -60, -30
    
    fig = plt.figure(figsize=(10,5))
    args = dict(color='gray',
                alpha=1.0, 
                linestyle='--', 
                linewidth=0.5,
                xlocs=np.arange(lonmin, lonmax, 10), 
                ylocs=np.arange(latmin, latmax, 10), 
                draw_labels=True)
    
    
    # Criação da projeção cartesiana

    ax = plt.axes(projection=ccrs.PlateCarree())
    gl = ax.gridlines(crs=ccrs.PlateCarree(), **args)
    gl.xlabel_style = {'size': 15, 'color': 'black', 'rotation': 0}
    gl.ylabel_style = {'size': 15, 'color': 'Gray', 'rotation': 0}
    gl.top_labels = False
    gl.right_labels = False
    gl.bottom_labels = True
    gl.rotate_labels = True
    ax.set_extent([lonmin, lonmax, latmin, latmax], ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND)
    ax.coastlines(resolution='50m', color='black', linewidth=1)


    #one_ciclone = data[(data.N_do_Ciclone == 200)].copy()
    def plot_cyclones(ciclone):

        for i in range(len(ciclone)-1):
            #print(range(len(ciclone)))
            #print(i)
            ciclone1 = ciclone.iloc[i]
            ciclone2 = ciclone.iloc[i+1]
            if i == 0:
                frist_lat, frist_lon = ciclone1.Latitude, ciclone1.Longitude
                second_lat, second_lon = ciclone2.Latitude, ciclone2.Longitude
                plt.plot([frist_lon, second_lon], [frist_lat, second_lat],
                          color='red',linewidth=0.5,
                          transform=ccrs.PlateCarree()
                          )
            else:
                if i == len(ciclone)-2:
                    frist_lat, frist_lon = ciclone1.Latitude, ciclone1.Longitude
                    second_lat, second_lon = ciclone2.Latitude, ciclone2.Longitude
                    plt.plot([frist_lon, second_lon], [frist_lat, second_lat],
                          color='black',linewidth=0.5,
                          transform=ccrs.PlateCarree()
                          )
                else:
                    frist_lat, frist_lon = ciclone1.Latitude, ciclone1.Longitude
                    second_lat, second_lon = ciclone2.Latitude, ciclone2.Longitude
                    plt.plot([frist_lon, second_lon], [frist_lat, second_lat],
                          color='blue',linewidth=0.5,
                          transform=ccrs.PlateCarree()
                          )



    ciclones = list( dict.fromkeys(df.N_do_Ciclone) )

    for Ciclone in ciclones:
        one_ciclone = df[(df.N_do_Ciclone == Ciclone)].copy()
        plot_cyclones(one_ciclone)
        
        
        
    plt.show()