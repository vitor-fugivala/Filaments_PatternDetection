from numpy import *
import matplotlib.pyplot as plt 
from scipy.spatial import Delaunay as dln
from shapely.geometry.polygon import Polygon
from scipy.interpolate import griddata 
import scipy.ndimage as ndi
import operator
import sys
from matplotlib.colors import LogNorm
from scipy.spatial import cKDTree
import pandas as pd
from matplotlib import patches
from projection_radec_to_xy import proj
#***********************************************************

def rebin(a, shape):
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)

def plota_n(grd,N,dim1,dim2,titulo,save,log,x,xloc,xlab,y,yloc,ylab):
  #log=0 => sem escala log
  #log=1 => com escala log
  plt.style.use('seaborn-muted')
  plt.figure(figsize=(10.24,10.24))

  if (log == 1):
    gmin=nanmin(grd[grd!=nanmin(grd)])
    gmax=nanmax(grd)
    if (gmin < 0):
      gmin=0.000000000001

  for n in range(N):
    plt.subplot(dim1,dim2,n+1)
    if (log == 1):
      plt.imshow(grd[n],cmap='gray_r', origin='lower',norm=LogNorm(vmin=gmin, vmax=gmax),extent=(min(x),max(x),min(y),max(y)),aspect='equal')
    else:
      plt.imshow(grd[n],cmap='gray_r', origin='lower',extent=(min(x),max(x),min(y),max(y)),aspect='equal')
    plt.title('n = {:01d}'.format(n), size='large',y=0.9, x=0.12,weight='bold')
    plt.xlim(min(x), max(x))
    plt.ylim(min(y), max(y))
    if (n == 0 ) or (n == 3) or (n == 6):
      plt.yticks(y_loc,y_labf, size='large')
    else:
      plt.yticks([])
    if (n == 6 ) or (n == 7) or (n == 8):
      plt.xticks(x_loc,x_labf,rotation=90, size='large')
    else:
      plt.xticks([])
    if (n==3):
      plt.ylabel('DEC',labelpad=4,fontsize='large', fontstyle='italic')
    if (n==7):
      plt.xlabel('RA',labelpad=4.5,fontsize='large', fontstyle='italic')

  plt.tick_params(direction='in',labelsize=12,top=True, right=True,labeltop=False, labelright=False)
  plt.tight_layout()
  plt.savefig(save)
  plt.close()

def plota_circmask(grd, xloc,xlab,yloc,ylab,titulo,save):
  #grid1 = blob
  plt.style.use('seaborn-muted')
  circ = patches.Circle((180,180), 180, facecolor='none')
  f=plt.figure(figsize=(10.24, 10.24))
  sb1=f.add_subplot(111)
  sb1.add_patch(circ)
  im1=sb1.imshow(grd, aspect='equal', origin='lower',cmap='gray_r',clip_path=circ, clip_on=True)
  plt.xlabel('RA', fontsize='large', fontstyle='italic')
  plt.ylabel('DEC', labelpad=0.01, fontsize='large', fontstyle='italic')
  plt.xticks(xloc,xlab,fontsize='large')
  plt.yticks(yloc,ylab,fontsize='large',horizontalalignment='right')
  plt.tick_params(direction='in',labelsize=12,top=True, right=True,labeltop=False, labelright=False)
  # plt.suptitle(titulo, size='large')
  plt.savefig(save)
  plt.close()


#***********************************************************

#ENTRADAS

# file_in = sys.argv[1]
# fov = sys.argv[2]
file_in='data_med_mapa_prob.csv' #diametro da imagem em segundos de arco

fil = genfromtxt(file_in, delimiter=',', unpack=True, skip_header=True)
# fov = 21600 #segundos de arco
fov=6 #graus

ra_in = fil[0] # RA
dec_in = fil[1] # DEC
z = fil [4]
massas = fil[5] # SOLAR MASS UNIT
prob = fil[6]

#------------------------------------------------------------
# transforma as entradas e projeta em x,y (pixel coordinates)
ra=[]
for item in ra_in:
  ra.append(float(item))
ra=asarray(ra)


dec=[]
for item in dec_in:
  dec.append(float(item))
dec=asarray(dec)

x, y = proj(ra, dec)

#------------------------------------------------------------
estrutura=[0,1]
le_pronto=0 #1=sim, 0 = nao

#estrutura = 0 => blob
#estrutura = 1 => filament

for est in estrutura:
  if (est == 1):
    x = xx_fi
    y = yy_fi

  ptos = [x,y]
  ptos = asarray(ptos)
  ptos = ptos.T

  tri=dln(ptos)

  trian=tri.simplices
  trian_ptos=tri.points
  num_tri=len(trian)
  num_ptos=len(trian_ptos)

#-----------------------------------------------------------
#axis thick in ra, dec
  x_lab=linspace(max(ra), min(ra), 7)
  x_loc=linspace(min(x), max(x), 7)
  y_lab=linspace(min(dec), max(dec), 7)
  y_loc=linspace(min(y),max(y), 7)

  ult=len(x_lab) - 1

  x_lab=delete(x_lab,ult)
  x_loc=delete(x_loc,ult)
  y_lab=delete(y_lab,ult)
  y_loc=delete(y_loc,ult)

  x_lab=delete(x_lab,0)
  x_loc=delete(x_loc,0)
  y_lab=delete(y_lab,0)
  y_loc=delete(y_loc,0)

  x_labf=[]
  y_labf=[]
  for i in range(len(x_lab)):
    x_labf.append(str(round(x_lab[i],1))+'°')
    y_labf.append(str(round(y_lab[i],1))+'°')


  #plota a tesselagem com os pontos
  if (est == 0):
    titulo='Delaunay tesselation - all'
    save='tessel_all.png'
  else:
    titulo='Delaunay tesselation - no blobs'
    save='tessel_noblobs.png'

  plt.style.use('seaborn-muted')
  plt.figure(figsize=(10.24, 10.24))
  plt.subplot(111)
  plt.triplot(x,y,tri.simplices, color='grey')
  plt.plot(x,y,marker='.',linestyle='',markersize=0.5, color='m' )
  # plt.title(titulo,size='large')
  plt.xlim(min(x), max(x))
  plt.ylim(min(y), max(y))

  plt.xlabel('RA',fontsize='large', fontstyle='italic')
  plt.ylabel('DEC',fontsize='large', fontstyle='italic')
  plt.xticks(x_loc,x_labf,fontsize='large')
  plt.yticks(y_loc,y_labf,fontsize='large')
  plt.tick_params(direction='in',labelsize=12,top=True, right=True,labeltop=False, labelright=False)

  # plt.subplots_adjust(bottom=0.16, right=0.8, top=0.9, left=0.2)

  plt.savefig(save)
  plt.close()
#-----------------------------------------------------------
  if (le_pronto == 1):
    if (est == 0):
      file='area_densi_all.csv'
    else:
      file='area_densi_noblob.csv'
    ja_calc=pd.read_csv(file)
    area_vor=ja_calc['area_voronoi']
    area_vor=area_vor.tolist()
    area_vor=asarray(area_vor)
    pho = ja_calc['densidade']
    pho=pho.tolist()
    pho=asarray(pho)
  else:

  #AREA DOS TRIANGULOS
    coord_groups = [tri.points[x] for x in tri.simplices]
    polygons = [Polygon(x) for x in coord_groups]

    area_tri=[]
    for i in range(num_tri):
      aux=polygons[i].area
      area_tri.append(aux)



  #-----------------------------------------------------------
  #AREA VORONOI
    area_vor=[]

    for p in range(num_ptos):
      contem=[]
      coords=tri.points[p]
      for t in range(num_tri):
        if p in trian[t]:
           contem.append(t)

      num_vor=len(contem)
      aux_vor=0
      for v in contem:
        aux_vor= aux_vor + area_tri[v]
      area_vor.append(aux_vor)
      print(p)

  #-----------------------------------------------------------
  #CALCULA A DENSIDADE DE CADA PONTO 

    pho=[]
    for p in range(num_ptos):
      aux_m=(float(massas[p])*3)/float(area_vor[p])
      aux_m=aux_m*prob[p]
      pho.append(aux_m)
      print(p)

    pho=asarray(pho)

    #salva numa file
    df=pd.DataFrame({'area_voronoi' : area_vor, 'densidade' : pho })
    if (est==0):
      file='area_densi_all.csv'
    else:
      file='area_densi_noblob.csv'
    df.to_csv(file)



#************************************************************
#INTERPOLACAO

#cria grid que sera interpolado
#grid quadrado com lado=diametro da regiao
  n_grid=1440
  grid_x, grid_y = mgrid[min(x):max(x):1440j, min(y):max(y):1440j]


  # #------------------------------------------------------------------
  # #Sbson interpolation

  #cria uma arvore com os pontos novos

  tree = cKDTree(c_[grid_x.ravel(), grid_y.ravel()])

  num_P=n_grid*n_grid

  c_P=zeros(num_P) #somas das areas
  n_P=zeros(num_P)# numero de galaxias perto
  # pho_P=zeros(num_P) # pho_P=n/c


  for i in range(len(trian_ptos)): 
    r, index=tree.query(trian_ptos[i],k=1) #distancia ao ponto P mais proximo
    in_circ = tree.query_ball_point(trian_ptos[i], 2*r) #busca todos Ps dentro de umas esfre
  # de 2r, centrada na galáxia em questão
    print(i)

    for item in in_circ:
      n_P[item]=n_P[item]+1
      c_P[item]= c_P[item] + pho[i]


  pho_P=c_P/n_P

  pho_P=pho_P.reshape(n_grid,n_grid)

  pho_Pc=nan_to_num(pho_P)
  #----------------------------------------------------------

  n_campo=360 #numero de pixeis no campo final
  form=[n_campo,n_campo]

  # campo=rebin(grid_z0, form)

  campo=rebin(pho_Pc, form)
  campo=campo.T

  campo_min=nanmin(campo[campo!=nanmin(campo)])
  campo_max=nanmax(campo)

  #normaliza
  # campo_norm=(campo-campo_min)/(campo_max - campo_min)

  #norm_log
  campo_normlog=log10(campo/campo_min)

  place(campo_normlog, campo_normlog<0, 0)


  #plot
  if (est == 0):
    titulo='Density field - all'
    save='densi_all.png'
  else:
    titulo='Density field - no blobs'
    save='densi_noblobs.png'

  plt.figure(figsize=(10.24, 10.24))
  plt.imshow(campo_normlog,cmap='gray_r', origin='lower',extent=(min(x),max(x),min(y),max(y)), interpolation='None',aspect='equal')
  plt.colorbar()
  # plt.title(titulo)
  plt.xlim(min(x), max(x))
  plt.ylim(min(y), max(y))
  plt.xlabel('RA',fontsize='large', fontstyle='italic')
  plt.ylabel('DEC',fontsize='large', fontstyle='italic')
  plt.xticks(x_loc,x_labf,fontsize='large')
  plt.yticks(y_loc,y_labf,fontsize='large')
  plt.tick_params(direction='in',labelsize=12,top=True, right=True,labeltop=False, labelright=False)
  plt.savefig(save)
  plt.close()

#------------------------------------------
#calcula o pixel scale
 #pizel scale = field of view em arcsec/ num de pixeis da largura da imagem

# num_pix=(shape(img_2d))[0]
  R0=fov/n_campo 
  R0=R0*60


  #R0=0.2

  #--------------------------------------------------
  #FILTRO GAUSSIANO NA IMAGEM ORIGINAL, VARIANDO O SIGMA DA GAUSS
  #sigma=Rn

  N=9
  fn=[0]*N
  fn_xx=[0]*N  #dado convoluido com a segunda derivada da gaussiana
  fn_xy=[0]*N
  fn_yy=[0]*N

  f_norm=[0]*N
  f_norm_xx=[0]*N #normalizado com Rn²
  f_norm_xy=[0]*N
  f_norm_yy=[0]*N



  n=0
  for n in range(N):
    Rn=((sqrt(2))**n)*R0
    #--------------------------------------------------------
    #FILTRO GAUSSIANO
    aux_f=ndi.gaussian_filter(campo_normlog, sigma=Rn)
    fn[n]=aux_f

    aux_fn=ndi.gaussian_filter(campo_normlog, sigma=Rn, order=[2,0])
    fn_xx[n]=aux_fn

    aux_fn=ndi.gaussian_filter(campo_normlog, sigma=Rn, order=[1,1])
    fn_xy[n]=aux_fn

    aux_fn=ndi.gaussian_filter(campo_normlog, sigma=Rn, order=[0,2])
    fn_yy[n]=aux_fn

    #--------------------------------------------------------
    #NORMALIZACAO
    f_aux=fn[n]*(Rn**2)
    f_norm[n]=f_aux

    f_aux=fn_xx[n]*(Rn**2)
    f_norm_xx[n]=f_aux

    f_aux=fn_xy[n]*(Rn**2)
    f_norm_xy[n]=f_aux

    f_aux=fn_yy[n]*(Rn**2)
    f_norm_yy[n]=f_aux

  if (est == 0):
    titulo='$\mathit{f_n}$: Gaussian Smoothed density field - all'
    save='gauss_all.png'
  else:
    titulo='$\mathit{f_n}$: Gaussian Smoothed density field - no blobs'
    save='gauss_noblob.png'

  plota_n(fn,9,3,3,titulo,save,0,x,x_loc,x_labf,y,y_loc,y_labf)



#--------------------------------------------------------
#AUTOVALORES
  l1_map=zeros((N,n_campo,n_campo))
  l2_map=zeros((N,n_campo,n_campo))

  n=0
  i=0
  j=0
  for n in range(N):
    for i in range(n_campo):
      for j in range(n_campo):

          Hxx=f_norm_xx[n][i][j]
          Hxy=f_norm_xy[n][i][j]
          Hyy=f_norm_yy[n][i][j]

          matrix_H=[[Hxx,Hxy],[Hxy,Hyy]]
          autoval=linalg.eigvals(matrix_H)
          autoval=autoval.tolist()

          autoval.sort(reverse=True)
          l1_map[n][i][j]=autoval[0]
          l2_map[n][i][j]=autoval[1]

  if (est == 0):
    titulo='$\mathit{\lambda_n^1}$: ${1^{st}}$ Eigenvalue Maps - all'
    save='lamb1_all.png'
  else:
    titulo='$\mathit{\lambda_n^1}$: ${1^{st}}$ Eigenvalue Maps - no blobs'
    save='lamb1_noblob.png'

  plota_n(l1_map,9,3,3,titulo,save,0,x,x_loc,x_labf,y,y_loc,y_labf)


  if (est == 0):
    titulo='$\mathit{\lambda_n^2}$: ${2^{st}}$ Eigenvalue Maps - all'
    save='lamb2_all.png'
  else:
    titulo='$\mathit{\lambda_n^2}$: ${2^{st}}$ Eigenvalue Maps - no blobs'
    save='lamb1_noblob.png'

  plota_n(l2_map,9,3,3,titulo,save,0,x,x_loc,x_labf,y,y_loc,y_labf)


#--------------------------------------------------------------
#MAPA MORFOLOGICO

#primeiro para esferas: os tres autovals tem que ser menores que 0
#se assim for, o pixel fica com valor 1

  Eps=zeros((N,n_campo,n_campo))

  n=0
  i=0
  j=0
  for n in range(N):
    for i in range(n_campo):
      for j in range(n_campo):
          # vxx_r=autoval_xx_r[n][i]
          # vxy_r=autoval_xy_r[n][i]
          # vyy_r=autoval_yy_r[n][j]
        l1=l1_map[n][i][j]
        l2=l2_map[n][i][j]
        rat=(l1-l2)/l1

              #BLOB
        if (est == 0):
          if (l1 < 0) and (l2 < 0):
            if (rat <=0.1):

                  # if (vyy_i <= 0j) and (vxy_i <= 0) and (vxx_i <= 0):
              Eps[n][i][j] = 1

              #FILAMENTO
        if (est == 1):
          if (l2 < 0):
            if (rat >=1):
              Eps[n][i][j] = 1

  if (est == 0):
    titulo='$\mathit{\epsilon_{blob}}$:  Blob Morphology Mask'
    save='eps_blob.png'
  else:
    titulo='$\mathit{\epsilon_{fi}}$:  Filament Morphology Mask'
    save='eps_fil.png'

  plota_n(Eps,9,3,3,titulo,save,0,x,x_loc,x_labf,y,y_loc,y_labf)


#--------------------------------------------------------------
# FIDELIDADE DA FORMA CARACTERISTICA = S

  S=zeros((N,n_campo,n_campo))

  n=0
  i=0
  j=0
  for n in range(N):
    for i in range(n_campo):
        for j in range(n_campo):

          l1=l1_map[n][i][j]
          l2=l2_map[n][i][j]

          #BLOB
          if (est ==0):
            if (l1 == 0.) and (l2 ==0.):
              S[n][i][j]=0.
            elif (l1 == 0):
              l1=1e-6
              S[n][i][j]=abs(l2/l1)
            else:
              S[n][i][j]=abs(l2/l1)

          #FILAMENTOS
          if (est == 1):
            if (l1 == 0.) and (l2 ==0):
              S[n][i][j]=0.
            elif (l1 == 0):
              l1=1e-6
              S[n][i][j]=(1-abs(l2/l1))
            else:
              S[n][i][j]=(1-abs(l2/l1))


#TA DANDO INFINITO EM ALGUMAS CONTAS

# S_2=nan_to_num(S)


  if (est == 0):
    titulo='$\it{S_{blob}}$: Blob Shape Significance'
    save='S_blob.png'
  else:
    titulo='$\it{S_{fi}}$: Filament Shape Significance'
    save='S_fil.png'

  plota_n(S,9,3,3,titulo,save,1,x,x_loc,x_labf,y,y_loc,y_labf)



#--------------------------------------------------------------
# FILTRO DE RESPOSTA MORFOLOGICA= M

  beta=0.5

  M=1-exp(-S/2*(beta**2))
  M_2=nan_to_num(M)

  if (est == 0):
    titulo='$\it{M_{blob}}$: Blob Morphology Map'
    save='M_blob.png'
  else:
    titulo='$\it{M_{fi}}$: Filament Morphology Map'
    save='M_fil.png'

  plota_n(M,9,3,3,titulo,save,1,x,x_loc,x_labf,y,y_loc,y_labf)

#--------------------------------------------------------------
# MAPA DE INTENSIDADE MORFOLOGICA = I


  I=zeros((N,n_campo,n_campo))

  n=0
  i=0
  j=0
  for n in range(N):
    for i in range(n_campo):
        for j in range(n_campo):
          l1=l1_map[n][i][j]
          l2=l2_map[n][i][j]

          #BLOB
          if (est == 0):
            I[n][i][j]=l2

          #FILAMENTO
          if (est == 1):
            I[n][i][j]=l1


    #normalizando
    maxv=nanmax(I[n])
    minv=nanmin(I[n])

    I[n]=(I[n] - minv)/(maxv - minv)


  if (est == 0):
    titulo='$\it{I_{blob}}$: Blob Morphology Intensity Map'
    save='I_blob.png'
  else:
    titulo='$\it{I_{fi}}$: Filament Morphology Intensity Map'
    save='I_fil.png'

  plota_n(I,9,3,3,titulo,save,1,x,x_loc,x_labf,y,y_loc,y_labf)



#--------------------------------------------------------------
# FILTRO MORFOLOGICO = T

  T=I*M_2
  T_2=nan_to_num(T)


  if (est == 0):
    titulo='$\it{T_{blob}}$: Blob Morphology Filter'
    save='T_blob.png'
  else:
    titulo='$\it{T_{fi}}$: Filament Morphology Filter'
    save='T_fil.png'

  plota_n(T,9,3,3,titulo,save,1,x,x_loc,x_labf,y,y_loc,y_labf)


#--------------------------------------------------------------
# MAPA DE FORMA = F
  F=Eps*T_2

  if (est == 0):
    titulo='$\it{F_{blob}}$: Blob Feature Map'
    save='F_blob.png'
  else:
    titulo='$\it{F_{fi}}$: Filament Feature Map'
    save='F_fil.png'

  plota_n(F,9,3,3,titulo,save,1,x,x_loc,x_labf,y,y_loc,y_labf)


#--------------------------------------------------------------
#MAPA DE SCALE SPACE STACKEADO = PHI

  PHI=zeros((n_campo,n_campo))

  n=0
  i=0
  j=0
  for i in range(n_campo):
    for j in range(n_campo):
      vals=[]
      for n in range(N):
        vals.append(F[n][i][j])

      PHI[i][j]=max(vals)

  maxv=nanmax(PHI)
  minv=nanmin(PHI)

  PHI_n=(PHI - minv)/(maxv - minv)


#ARRUMA POSICAO DOS TICKS-----------------------
  x_loc_p=linspace(0,360, 7)
  ult=len(x_loc_p) - 1
  x_loc_p=delete(x_loc_p,ult)
  x_loc_p=delete(x_loc_p,0)
  y_loc_p=x_loc_p

  if (est == 0):
    titulo='$\it{\Psi_{blob}}$: Blob Scale space map stack'
    save='phi_blob.png'
  else:
    titulo='$\it{\Psi_{fi}}$: Filament Scale space map stack'
    save='phi_fil.png'

  plota_circmask(PHI_n, x_loc_p,x_labf,y_loc_p,y_labf,titulo,save)

#-------------------------------------------
# LIMITES
  if (est == 0):
    save='blob_PHI_hist.png'
  else:
    save='fi_PHI_hist.png'


  hn=PHI_n
  plt.hist(hn)
  plt.savefig(save)
  plt.close()

# depois de analizar essas saidas, devemos escolher os limites

  n=2
  Rn=((sqrt(2))**n)*R0

  lims=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]

  i=0
  j=0

  for a in range(len(lims)-1):
    PHI_final=zeros((n_campo,n_campo))
    for i in range(n_campo):
      for j in range(n_campo):
        if (PHI_n[i][j] > lims[a]) and (PHI_n[i][j] <= lims[a+1]):
          PHI_final[i][j] = (PHI_n[i][j])

          if (est == 0):
            save=str('blob-' + str(lims[a]) + '-' + str(lims[a+1]) + '.png')
          else:
            save=str('fi-' + str(lims[a]) + '-' + str(lims[a+1]) + '.png')

    PHI_final_gaus=ndi.gaussian_filter(PHI_final, sigma=Rn)
    plota_circmask(PHI_final_gaus, x_loc_p,x_labf,y_loc_p,y_labf,titulo,save)



#todos os valores abaixo disso vão virar zero
  if (est == 0): 
    lim_inf=0.
    lim_sup=0.2
  else:
    lim_inf=0.5
    lim_sup=0.7    


  n=2
  Rn=((sqrt(2))**n)*R0

  i=0
  j=0

  PHI_final=zeros((n_campo,n_campo))
  for i in range(n_campo):
    for j in range(n_campo):
      if (PHI_n[i][j] >= lim_inf) and (PHI_n[i][j] <= lim_sup):
        PHI_final[i][j] = (PHI_n[i][j])

  PHI_final_gaus=ndi.gaussian_filter(PHI_final, sigma=Rn)


  if (est == 0):
    titulo='$\it{\Psi_{blob}}$: Clean Blob Scale space map stack'
    save='final_phi_blob.png'
  else:
    titulo='$\it{\Psi_{fi}}$: Clean Filament Scale space map stack'
    save='final_phi_fil.png'

  plota_circmask(PHI_final_gaus, x_loc_p,x_labf,y_loc_p,y_labf,titulo,save)




#---------------------------------------------------------------------
#ASSOCIANDO AOS PONTOS REAIS

#cria grid com 0s e 1s, sendo 1 o pixel que contem galaxia observada
#cria grid com o tamanho da observação e com numero de bins do campo
  if (est == 0):
    ehgal=0.03
    grid_gal_esf=zeros((n_campo,n_campo))
    realdet_esf=zeros(len(x))

  else:
    ehgal=0.25
    grid_gal_fi=zeros((n_campo,n_campo))
    realdet_fi=zeros(len(x))

  grid_x_gal,grid_y_gal=mgrid[min(x):max(x):361j, min(y):max(y):361j]
  x_gal_v=grid_x_gal[:,0]
  y_gal_v=grid_y_gal[0,:]

  #cria novo grid, com o valor medio dos bins anteriores 
  #ao inves de ter os vertices,terei o centro dos bins

  x_gal_c=zeros(n_campo)
  y_gal_c=zeros(n_campo)

  for i in range(n_campo-1):
    x_gal_c[i]=( x_gal_v[i+1] + x_gal_v[i])/2
    y_gal_c[i]=( y_gal_v[i+1] + y_gal_v[i])/2

    # PHI_val=[]
    # for i in range(len(x_gal))

  #pos_grid_gal guarda a posicao no grid de cada galáxia, em index
  pos_grid_gal=zeros([len(x),2])
  for i in range(len(x)):
    minx=(abs(x_gal_c-x[i])).argmin()
    miny=(abs(y_gal_c-y[i])).argmin()

    pos_grid_gal[i,0]=minx
    pos_grid_gal[i,1]=miny

  
  Pfg=PHI_final_gaus.T
  for g in range(len(pos_grid_gal)):
    i=int(pos_grid_gal[g,0])
    j=int(pos_grid_gal[g,1])
    print(i)
    print(j)
    pg=PHI_final_gaus[i][j]
    pgT=Pfg[i][j]
    print(pg)

    if (est == 0):
      if (pg >= ehgal):
        grid_gal_esf[i][j]=1.
      if (pgT >= ehgal):
        realdet_esf[g] = 1
    else:
      if (pg >= ehgal):
        grid_gal_fi[i][j]=1
        realdet_fi[g] = 1


  if (est == 0):
    titulo='$\it{O_{blob}}$: Blob real detection'
    save='realdet_blob.png'
    plota_circmask(grid_gal_esf, x_loc_p,x_labf,y_loc_p,y_labf,titulo,save)
  else:
    titulo='$\it{O_{fi}}$: Filament real detection'
    save='realdet_fil.png'
    plota_circmask(grid_gal_fi, x_loc_p,x_labf,y_loc_p,y_labf,titulo,save)




  
  if (est == 0):
    xx_fi = []
    yy_fi = []
    x_esf = []
    y_esf = []
    raa_esf=[]
    decc_esf=[]
    for r in range(len(realdet_esf)):
      if (realdet_esf[r] == 0.):
        xx_fi.append(x[r])
        yy_fi.append(y[r])
      else:
        x_esf.append(x[r])
        y_esf.append(y[r])
        raa_esf.append(ra[r])
        decc_esf.append(dec[r])
  else:
    x_fi=[]
    y_fi=[]
    raa_fi = []
    decc_fi = []
    for r in range(len(realdet_fi)):
      if (realdet_fi[r] == 1.):
        x_fi.append(x[r])
        y_fi.append(y[r])
        raa_fi.append(ra[r])
        decc_fi.append(dec[r])


#plota junto
plt.figure(figsize=(10.24, 10.24))
plt.subplot(111)
plt.plot(x_esf,y_esf,marker='.',linestyle='',markersize=1, color='grey' )
plt.plot(x_fi,y_fi,marker='.',linestyle='',markersize=1, color='k' )
# plt.title('$\it{O_{tot}}$: Total real detection')
plt.xlim(min(x), max(x))
plt.ylim(min(y), max(y))

plt.xlabel('RA', fontstyle='italic',fontsize='large')
plt.ylabel('DEC', fontstyle='italic',fontsize='large')
# plt.subplots_adjust(bottom=0.16, right=0.8, top=0.9, left=0.2)
plt.tick_params(direction='in',labelsize=12,top=True, right=True,labeltop=False, labelright=False)

plt.xticks(x_loc,x_labf,fontsize='large')
plt.yticks(y_loc,y_labf,fontsize='large')

plt.savefig('realdet_tot.png')
plt.close()

#SALVA PONTOS IDENTIFICADOS
raa_fi=asarray(raa_fi)
decc_fi=asarray(decc_fi)
raa_esf=asarray(raa_esf)
decc_esf=asarray(decc_esf)
x_fi=asarray(x_fi)
y_fi=asarray(y_fi)
x_esf=asarray(x_esf)
y_esf=asarray(y_esf)

columns='ra_fil, dec_fil, x_fil, y_fil'
cols=[raa_fi, decc_fi, x_fi, y_fi]
nome_sai3='Realdet_fil.txt'
savetxt(nome_sai3, transpose(cols), header=columns, delimiter=',', fmt='%1.6e')

columns='ra_cls, dec_cls, x_cls, y_cls'
cols=[raa_esf, decc_esf, x_esf, y_esf]
nome_sai3='Realdet_blob.txt'
savetxt(nome_sai3, transpose(cols), header=columns, delimiter=',', fmt='%1.6e')








