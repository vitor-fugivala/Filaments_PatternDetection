# Natalia Del Coco

# Detecta filamentos pelo método SCMS
#==================================================================
import sys
import pandas as pd
import numpy as np
import os
import subprocess
from scipy.spatial import distance
from astropy.stats import bootstrap as boots_func
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import tempfile


#==================================================================
#                             FUNÇÕES
#==================================================================

# Executa o algoritmo Subspace Constrained Mean Shift (SCMS).
def SCMS(data, SCMS_code, lin, h, iter, e):

    # Cria arquivo temporário contendo as coordenadas.
    with tempfile.NamedTemporaryFile(delete=False, suffix=".txt") as tmp_in:
        in_btsp = tmp_in.name
    data2 = np.column_stack([data[0], data[1]])
    np.savetxt(in_btsp, data2, fmt=['%f', '%f'])

    out_btsp = os.getcwd() + f'/{lin}.txt'
    
    try:
        subprocess.run([SCMS_code, '-d', in_btsp, '-h', h, '-i', iter, '-e', e, '-n', out_btsp])
        aux_R_btsp = np.loadtxt(out_btsp)
    # Deleta o arquivo temporário
    finally:
        if os.path.exists(in_btsp):
            os.remove(in_btsp)

    return aux_R_btsp

#==================================================================
#                               MAIN
#==================================================================
# ---------------- LÊ ENTRADAS ----------------
# RODA NO DIRETÓRIO ONDE ESTÃO OS DADOS
SCMS_code_fil = os.getcwd() + '/find_filament_2d'
in_dt  = sys.argv[0]        # Dados originais, com as posições de todas as galáxias.
A0_fil = float(sys.argv[1]) # Para a largura da gaussiana: h².
iter   = str(sys.argv[2])   # Max iterações.
e      = str(sys.argv[3])   # Err mínimo.

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
print('')
print('')
print('========== DETECTA FILAMENTOS ==========')

# Cálculo da Eq. (A1) de Chen, 2015.
h_fil = A0_fil*((1/(d + 2))**(1/(d + 4)))*((len(ra))**(-1/(d + 4)))*sig
h_fil = str(h_fil)

# Roda o SCMS com o valor determinado de h.
detect_fil = SCMS(np.array([ra, dec]), SCMS_code_fil, 'realfil', h_fil, iter, e)
print('Saída filamentos detectados: realfil.txt')

# Plota filamentos detectados
dfil    = np.loadtxt('realfil.txt')
ra_fil  = dfil[:, 0]
dec_fil = dfil[:, 1]
fig, ax = plt.subplots()
ax.scatter(ra, dec, color='darkgrey', s=1)
ax.scatter(ra_fil, dec_fil, color='red', s=1)
ax.set_xlabel('RA (deg)')
ax.set_ylabel('DEC (deg)')
ax.set_aspect('equal')
plt.gca().invert_xaxis()
fig.savefig('filaments.png', bbox_inches='tight')
print('Imagem dos filamentos: filaments.png')
