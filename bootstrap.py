# Natalia Del Coco

# Acha incerteza nos filamentos detectados por SCMS
#========================================================
import sys
import pandas as pd
import numpy as np
import os
import subprocess
from astropy.stats import bootstrap as boots_func
from tqdm import tqdm
import multiprocessing
from functools import partial
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import argparse


#==================================================================
#                             FUNÇÕES
#==================================================================

#-----------------------------------------------------------
# Calcula a distância ao ponto mais próximo de "node".
# Entradas: node (ponto testado), nodes (conjunto de pontos)
def closest_distance(node, nodes):
    nodes = np.asarray(nodes)
  # Calcula distâncias euclidianas e encontra o menor valor
    euclidean_dist = np.sqrt(np.sum((nodes - node)**2, axis=1))
    closest_dist   = np.amin(euclidean_dist)
    return closest_dist

#-----------------------------------------------------------
# Executa o algoritmo Subspace Constrained Mean Shift (SCMS)
# para amostras de bootstrap.
def SCMS(data, SCMS_code, lin, h, inter, e):
    
    # Salva em arquivo .txt a amostra bootstrap de entrada.
    in_btsp  = os.getcwd() + '/bootstrap_sample_'+'%s.txt' % (lin+1)
    np.savetxt(in_btsp, data, fmt=['%f', '%f'])

    # Salva em arquivo .txt os filamentos detectados.
    out_btsp = os.getcwd() + '/bootstrap_fil_'+'%s.txt' % (lin+1)
    subprocess.run([SCMS_code, '-d', in_btsp, '-h', str(h),
                    '-i', str(inter),'-e', str(e), '-n', out_btsp, '-q'])
    aux_R_btsp = np.loadtxt(out_btsp)

    return aux_R_btsp

#--------------------------------------------------------------
# Função auxiliar para gerar amostras bootstrap e rodar o SCMS.
def process_bootstrap(idx, data, SCMS_code_fil, h_fil, iter, e):
    bsi = boots_func(data=data, bootnum=1)[0].tolist()
    result = SCMS(bsi, SCMS_code_fil, idx, h_fil, iter, e)
    return result


#==================================================================
#                               MAIN
#==================================================================
# ---------------- LÊ ENTRADAS ----------------
# RODA NO DIRETÓRIO ONDE ESTÃO OS DADOS
SCMS_code_fil = os.getcwd() + '/find_filament_2d'

# CONFIGURAÇÃO DO PARSER DE ARGUMENTOS
parser = argparse.ArgumentParser(
    description="Executa análise de bootstrap para encontrar incertezas em filamentos cosmológicos.",
    formatter_class=argparse.RawTextHelpFormatter
)

# Argumentos obrigatórios
parser.add_argument("input_file", type=str, help="Caminho para o arquivo de dados original com as posições das galáxias")
parser.add_argument("A0_fil", type=float, help="Parâmetro A_0 para a largura da gaussiana (h²). Ex: 0.8")
parser.add_argument("iterations", type=int, help="Número máximo de iterações para o SCMS. Ex: 1000")
parser.add_argument("min_error", type=float, help="Erro mínimo para convergência do SCMS. Ex: 1E-4")

# Argumentos opcionais
# Número de núcleos para paralelização
n_cores       = os.cpu_count()
default_jobs  = max(1, n_cores - 2)
parser.add_argument(
    "--njobs",
    "-nj",
    type=int,
    default=default_jobs,
    help=f"Número de núcleos a serem usados para paralelização. Padrão: {default_jobs} (total de cores - 2)"
)
# Número de amostras bootstrap
parser.add_argument(
    "--nboots",
    "-nb",
    type=int,
    default=100,
    help="Quantidade de amostras bootstrap desejadas. Padrão: 100"
)

# Processa os argumentos da linha de comando
args = parser.parse_args()

in_dt     = args.input_file
A0_fil    = args.A0_fil
iter      = args.iterations
e         = args.min_error
n_jobs    = args.njobs
Nbs_total = args.nboots

# ........ Lê pontos originais ........
data   = pd.read_csv(in_dt, delimiter='\t')
ra     = data.iloc[:, 0]    # Coordenadas RA.
dec    = data.iloc[:, 1]    # Coordenadas DEC.

# ........ Detecta estruturas .........
sigra  = np.std(ra)         # Desvio padrão em RA.
sigdec = np.std(dec)        # Desvio padrão em DEC.
sig    = np.min([sigra, sigdec]) # Sigma mínimo entre as coordenadas (ver Eq. A1 de Chen, 2015 [https://arxiv.org/abs/1501.05303])
d      = 2                  # Dimensão.

# ........ Filamentos .................
# Cálculo da Eq. (A1) de Chen, 2015.
h_fil = A0_fil*((1/(d + 2))**(1/(d + 4)))*((len(ra))**(-1/(d + 4)))*sig
h_fil = str(h_fil)

# Lê arquivos gerados----------------------------------------
dfil    = np.loadtxt('realfil.txt')
ra_fil  = dfil[:, 0]
dec_fil = dfil[:, 1]
fil_det = np.asarray([ra_fil, dec_fil]).T

print('\nFilamentos de entrada: realfil.txt\n')

#==================================================================
#                            BOOTSTRAP
#==================================================================
print('============= CÁLCULO DE ERROS =============\n')

print(f"USANDO {n_jobs} NÚCLEOS PARA O PARALELISMO.\n")

print(f"PROCESSANDO {Nbs_total} AMOSTRAS BOOTSTRAP.\n")

data = np.array(data)
indices = list(range(Nbs_total))

# Procedimento para rodar o SCMS em paralelo e imprimir
# barra de progresso conforme as amostras de bootstrap
# são processadas e os arquivos de filamentos salvos.
worker_func = partial(process_bootstrap, 
                      data=data, 
                      SCMS_code_fil=SCMS_code_fil, 
                      h_fil=h_fil, 
                      iter=iter, 
                      e=e)

fil_btsp = []

with multiprocessing.Pool(processes=n_jobs) as pool:
    print("Iniciando processamento bootstrap...")
    for resultado in tqdm(pool.imap_unordered(worker_func, indices),
                          total=len(indices),
                          desc="SCMS processed bootstrap samples",
                          ncols=80):
        fil_btsp.append(resultado)

print("\nPROCESSAMENTO BOOTSTRAP CONCLUÍDO!\n")

#==================================================================
#                         CÁLCULO DE ERROS
#==================================================================
# Cálculo da confiabilidade de cada ponto ŷ_0.
# Ver Eq. (4.16) dissertação Del Coco.

# Distâncias mínimas.
D = np.zeros(len(fil_btsp))
# Função de mérito para erros.
C = np.zeros(len(dfil))

# Calcula para todos os pontos.
for i in range(len(dfil)):
    for j in range(len(fil_btsp)):
        D[j] = closest_distance(dfil[i], np.array(fil_btsp[j]))
        #print(i, j, D[j])
    C[i] = np.sqrt(1./Nbs_total*np.sum(D))

# Cria dataframe para retornar filamentos e erros.
df_out_fil  = pd.DataFrame({
    'RA': dfil[:,0],
    'DEC': dfil[:,1],
    'Erro': C
    })

# Um filamento é considerado estável se o parâmetro de
# confiabilidade C é menor que um limiar C_threshold.
C_threshold = df_out_fil['Erro'].mean() + 1.69*df_out_fil['Erro'].std()
print('Limite de confiabilidade (C_threshold) = %f\n' % C_threshold)
# Se Rejeita == 1, o ponto deve ser desconsiderado, Se == 0, é válido.
df_out_fil['Rejeita'] = (df_out_fil['Erro'] - C_threshold >= 0).astype(int)

#----- ESCREVE OUTPUTS ------
df_outname='./finalfil.csv'
df_out_fil.to_csv(df_outname, index=False, sep=',')
print('Escrevendo filamentos finais para: %s' % df_outname)

#----- Plota filamentos confiáveis -----
df_stable = df_out_fil[df_out_fil['Rejeita'] == 0]
fig, ax = plt.subplots()
ax.scatter(data[:,0], data[:,1], color='black', s=1)
h = ax.scatter(df_stable['RA'], df_stable['DEC'], c=df_stable['Erro'], cmap=cm.magma_r, s=60*df_stable['Erro'])
fig.colorbar(h, ax=ax, label='Uncertainty')
plt.xlabel('RA (deg)')
plt.ylabel('DEC (deg)')
ax.set_aspect('equal')
plt.gca().invert_xaxis()
fig.savefig('finalfil.png', bbox_inches='tight')
print('Imagem dos filamentos com erros: ./finalfil.png\n')