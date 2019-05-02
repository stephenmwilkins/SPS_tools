

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import copy
import pickle

plt.style.use('simple')



fig = plt.figure(figsize = (4.5,3.5))

left  = 0.1
bottom = 0.1  
height = 0.8
width = 0.65

ax = fig.add_axes((left, bottom, width, height))

CloudyVersion = 'C17.00'



SPS = 'BPASSv2.2.binary'
IMF = 'ModSalpeter_300'
model_label = r'$\rm BPASSv2.2\quad Modified\,Salpeter\, (m_{\star}=0.1-300\,{\rm M_{\odot}})$'


line = ['HI4861','OIII4959','OIII5007']



fig = plt.figure(figsize = (3.5,3.5))

left  = 0.15
bottom = 0.15  
height = 0.8
width = 0.8

ax = fig.add_axes((left, bottom, width, height))

ax.axhline(0.0,c='k',alpha=0.1,lw=2)



grid = pickle.load(open('../../BuildGrid/Z_U_nH/grids/'+SPS+'/'+IMF+'/lines_SSP1.p','rb'), encoding='latin1')

inH_def = (np.abs(grid['log10n_H'] - 2.5)).argmin()
iU_def = (np.abs(grid['log10U_S_0'] - -2.0)).argmin()


print(grid['log10U_S_0'], iU_def)
print(grid['log10n_H'], inH_def)


inHs = [(np.abs(grid['log10n_H'] - log10n_H)).argmin() for log10n_H in [1., 2.,2.5, 3., 4.]]      
iUs = [(np.abs(grid['log10U_S_0'] - log10U_S_0)).argmin() for log10U_S_0 in [-1.,-2.,-3.]]      







c = 'k'




L = []
for iZ, Z in enumerate(grid['Z']):
    line_luminosity = 0.0
    for l in line:
        l = l.encode('UTF-8')
        line_luminosity += grid[l]['luminosity'][iZ, iU_def, inH_def]
    L.append(line_luminosity)        
L_def = np.log10(np.array(L))


for inH, lw in zip(inHs, [0.5, 0.75, 1, 1.5, 2.0]):

    L = []
    for iZ, Z in enumerate(grid['Z']):
        line_luminosity = 0.0
        for l in line:
            l = l.encode('UTF-8')
            line_luminosity += grid[l]['luminosity'][iZ, iU_def, inH]
        L.append(line_luminosity)        
    
    L = np.log10(np.array(L))

    ax.plot(grid['log10Z'], (L-L_def), c=c, ls='-', lw=lw, alpha = 0.5)

for iU, ls in zip(iUs,[':','-','--']):

    L = []
    for iZ, Z in enumerate(grid['Z']):
        line_luminosity = 0.0
        for l in line:
            l = l.encode('UTF-8')
            line_luminosity += grid[l]['luminosity'][iZ, iU, inH_def]
        L.append(line_luminosity)        
    
    L = np.log10(np.array(L))

    ax.plot(grid['log10Z'], (L-L_def), c=c, ls=ls, lw=1)

    




dummies = []
labels = []    

labels += [r'$\rm\log_{10}(U_{S,0})='+str(x)+'$' for x in [-1.,-2.,-3.]]
dummies += [ax.plot([], [], ls=ls, c='k', alpha=0.5, lw = 1)[0] for ls in [':','-','--']]        

labels += [r'$\rm\log_{10}(n_{H})='+str(x)+'$' for x in [1., 2.,2.5, 3., 4.]]
dummies += [ax.plot([], [], ls='-', c='k', alpha=0.2, lw = lw)[0] for lw in [0.5, 0.75, 1, 1.5, 2.0]]        


ax.legend(dummies, labels,loc = 'upper left', fontsize = 6, handlelength=2.5)


# fig.text(0.45, 0.04, r'$\rm\log_{10}(Z)$', ha='center', fontsize=6)
# fig.text(0.04, 0.55, r'$\rm\log_{10}(L/L_{default})$', va='center', rotation='vertical', fontsize=6)
# 

        

ax.set_xlabel(r'$\rm\log_{10}(Z)$')
ax.set_ylabel(r'$\rm\log_{10}[(L/L_{ref})_{1\ Myr}]$')


figname = 'figures/line_luminosity_effects_SINGLE.pdf'
print(figname)
fig.savefig(figname)

# ---- define plot



