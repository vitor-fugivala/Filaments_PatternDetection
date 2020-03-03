# Natalia Del Coco

#acha incerteza nos filamentos detectados por SCMS
#========================================================
import sys
import pandas as pd
import numpy as np
import os
import subprocess
from scipy.spatial import distance
from joblib import Parallel, delayed 
from astropy.stats import bootstrap as boots_func

#=======================================================
#FUNCOES
#=======================================================
# def bootstrap(X,n):
#   #X = dados de entrada
#   #n = comprimento do sample gerado
#   if n==None:
#     n=len(X)

#   resample_i = np.floor(np.random.rand(n)*len(X)).astype(int)
#   bootsample=X[resample_i]
#   return bootsample

def bootstrap(X,resample, sample):
  #X = dados de entrada
  #resample -  numeros de boots
  #sample = numeros de ptos dentro de cada resample

  bootsample = boots_func(X,resample,sample)
  return bootsample

# def smoothed_bootstrap(X,n,h):
#   #X = dados de entrada
#   #n = comprimento do sample gerado
#   if n==None:
#     n=len(X)

#   resample_i = np.floor(np.random.rand(n)*len(X)).astype(int)
#   bootsample=X[resample_i]

#   rnd=np.random.normal(loc=bootsample, sclae=h)

#   boot_final=bootsample
#   return bootsample

def closest_dist(node, nodes):
    closest_index = distance.cdist([node], nodes, metric='euclidean').argmin()
    closest_distance = np.sum((node-nodes[closest_index])**2)
    return closest_distance


def SCMS(data,SCMS_code,lin, h, inter, e, erro):
  # SCMS_code='/home/natalia/Codigos/SCMS/Filament_nati/find_filament_test'
  print(SCMS_code)

  if (erro ==0): #SIM, VAI CALCULAR ERRO
    in_btsp=os.getcwd()+'/bootstrap_sample_'+'%s.txt' % lin
    out_btsp=os.getcwd()+'/bootstrap_error_'+'%s.txt' % lin
  else: #ERRO = 1 => VAI CALCULAR PTOS REAIS
    in_btsp=os.getcwd()+'/sample_'+'%s.txt' % lin
    out_btsp=os.getcwd()+'/'+'%s.txt' % lin
  df=pd.DataFrame(data)
  print(df)

  df.to_csv(in_btsp, index=False, sep=' ', header=None)
  subprocess.run([SCMS_code,'-d',in_btsp,'-h',h,'-i',inter,'-e',e,'-n',out_btsp])
  aux_R_btsp=pd.read_csv(out_btsp, header=None, delimiter=',')
  return aux_R_btsp



#=======================================================
#MAIN
#=======================================================
 # ----- LE ENTRADAS ----------------
#RODAS NO DIRETORIO ONDE ESTAO OS DADOS
SCMS_code_fil='/home/natalia/Codigos/SCMS/Filament_nati/find_filament'
SCMS_code_cl='/home/natalia/Codigos/SCMS/Filament_nati/find_cluster'
path=os.getcwd()
# SCMS_out='./testfile_output_err.csv'


in_dt=sys.argv[1] #DADOS ORIGINAIS, COM AS POSICOES DE TODAS AS GALAXIAS
# in_R=sys.argv[2] # filamantos detectados do sample original (ra, dec)
A0_fil=float(sys.argv[2]) #para a largura da gaussiana
A0_cl=float(sys.argv[3]) #para a largura da gaussiana
inter=str(sys.argv[4]) # max iteracoes
e=str(sys.argv[5]) #err minimo
# in_dt='data_med_mapa_2.txt'
# in_R='testfile_output_F_A09.csv'

#le pontos originais
dt=pd.read_csv(in_dt, delimiter=' ', header=None)
ra=dt[0]
dec=dt[1]
ptos=np.asarray([ra,dec]).T

# dR=pd.read_csv(in_R, delimiter=',', header=None)
# ra_r=dR[0]
# dec_r=dR[1]
# R_tot=np.asarray([ra_r,dec_r]).T

#------ DETECTA ESTRUTURAS ------------
sigra=np.std(ra)
sigdec=np.std(dec)
sig=np.min([sigra,sigdec])
d=2
#filametos....................
print('')
print('')
print('====== DETECTA FILAMENTOS ==========')

h_fil=A0_fil*((1/(d+2))**(1/(d+4)))*((len(ra))**(-1/(d+4)))*sig
h_fil=str(h_fil)
detect_fil=SCMS(ptos,SCMS_code_fil,'realfil', h_fil, inter, e, 1)
print('Saída filamentos detectados: realfil.txt')

#clusters....................
print('')
print('')
print('====== DETECTA AGLOMERADOS ==========')

h_cl=A0_cl*((1/(d+2))**(1/(d+4)))*((len(ra))**(-1/(d+4)))*sig
h_cl=str(h_cl)
detect_cl=SCMS(ptos,SCMS_code_cl,'realcl', h_cl, inter, e, 1)
print('Saída aglomerados detectados: realcl.txt')

#le arquivos gerados----------------------------------------
dfil=pd.read_csv('realfil.txt', delimiter=',', header=None)
ra_fil=dfil[0]
dec_fil=dfil[1]
fil_det=np.asarray([ra_fil,dec_fil]).T

dcl=pd.read_csv('realcl.txt', delimiter=',', header=None)
ra_cl=dcl[0]
dec_cl=dcl[1]
cl_det=np.asarray([ra_cl,dec_cl]).T

#------ ERROR - GERA BOOTSTRAP ------------
#se eu quero NBS total de 100, preciso repetir o processo 
print('')
print('')
print('====== SCMS PARA ERROS ==========')

vezes=2
Nbs=2
fil_btsp=[]
cl_btsp=[]
for i in range(vezes):

  print('')
  print('')
  print('====== %s° RUN ==========' % i)

  bs=[]
  for j in range(Nbs):
    aux=smoothed_bootstrap(ptos,len(ptos))
    bs.append(aux.tolist())

  cts=10*i
  fil_btsp_aux = Parallel(n_jobs=11, verbose=1)(delayed(SCMS)(bsi,SCMS_code_fil,(bs.index(bsi)+cts), h_fil, inter, e, 0) for bsi in bs)
  cl_btsp_aux = Parallel(n_jobs=11, verbose=1)(delayed(SCMS)(bsi,SCMS_code_cl,(bs.index(bsi)+cts), h_cl, inter, e, 0) for bsi in bs)

  fil_btsp = fil_btsp + fil_btsp_aux
  cl_btsp = cl_btsp + cl_btsp_aux
  print('comprimento fil_btsp = %s' % len(fil_btsp))
  print('comprimento fil_btsp = %s' % len(fil_btsp))
  print('')

print('SCMS para erros: terminado')
# ----CRIA FILES E RODA SCMS -----
# #filametos....................
# fil_btsp = Parallel(n_jobs=, verbose=1)(delayed(SCMS)(bsi,SCMS_code_fil,(bs.index(bsi)), h_fil, inter, e, 0) for bsi in bs)
# #clusters....................
# cl_btsp = Parallel(n_jobs=2, verbose=1)(delayed(SCMS)(bsi,SCMS_code_cl,(bs.index(bsi)), h_cl, inter, e, 0) for bsi in bs)

#---- CALCULA ERROS ---------------
#filametos....................

print('')
print('')
print('====== CALCULO DOS ERROS ==========')

error_fil=[]

for x in fil_det:
  dist=0
  for i in range(Nbs):
    aux_dist = closest_dist(x,fil_btsp[i].values)
    dist = dist + aux_dist
  aux_err=np.sqrt(dist/Nbs)
  error_fil.append(aux_err)

#clusters....................
error_cl=[]

for x in cl_det:
  dist=0
  for i in range(Nbs):
    aux_dist = closest_dist(x,cl_btsp[i].values)
    dist = dist + aux_dist
  aux_err=np.sqrt(dist/Nbs)
  error_cl.append(aux_err)

#----- INSTABILIDADE -----------
#filametos....................
dif_fil=np.mean(error_fil)+1.69*np.std(error_fil)

estabi_fil=[]
sai_fil=[]
for err in error_fil:
  estabi_fil.append(err-dif_fil)
  if (err-dif_fil >= 0 ):
    sai_fil.append(1) #ou seja, o ponto deve ser desconsiderado
  else:
    sai_fil.append(0)

#clusters....................
dif_cl=np.mean(error_cl)+1.69*np.std(error_cl)

estabi_cl=[]
sai_cl=[]
for err in error_cl:
  estabi_cl.append(err-dif_cl)
  if (err-dif_cl >= 0 ):
    sai_cl.append(1) #ou seja, o ponto deve ser desconsiderado
  else:
    sai_cl.append(0)

#----- ESCREVE OUTPUTS ------
#filametos....................
print('')
print('')
print('====== ESCREVE SAIDAS ==========')

df_out_fil=pd.DataFrame([error_fil,estabi_fil,sai_fil]).T
cols=["Erro", "Estabilidade", "Rejeita"]
df_outname='./errorfil.csv'
df_out_fil.to_csv(df_outname, index=False, sep=',', header=cols)

print('Saída erro filamentos: %s' % df_outname)

#clusters....................
df_out_cl=pd.DataFrame([error_cl,estabi_cl,sai_cl]).T
cols=["Erro", "Estabilidade", "Rejeita"]
df_outname='./errorcl.csv'
df_out_cl.to_csv(df_outname, index=False, sep=',', header=cols)

print('Saída erro aglomerados: %s' % df_outname)


