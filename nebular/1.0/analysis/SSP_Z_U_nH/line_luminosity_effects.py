

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import copy
import pickle

plt.style.use('simple')



CloudyVersion = 'C17.00'



SPS = 'BPASSv2.2.binary'
IMF = 'ModSalpeter_300'
model_label = r'$\rm BPASSv2.2\quad Modified\,Salpeter\, (m_{\star}=0.1-300\,{\rm M_{\odot}})$'


line = ['HI6563']
# line = 'CIII1907,CIII1909'

lines = [['HI1216'],['CIII1909','CIII1907'],['OII3726','OII3729'],['NeIII3869'],['NeIII3967'],['HI4340'],['HI4861'],['OIII4959'],['OIII5007'],['HI6563'], ['NII6583'], ['SII6731', 'SII6716']]


labels = {}
labels['HI1216'] = r'Ly\alpha'
labels['CIII1909-CIII1907'] = r'\regular{[CII]1907,1909}'
labels['OII3726-OII3729'] = r'\regular{[OII]3726,3729}'
labels['NeIII3869'] = r'\regular{[NeIII]3869}'
labels['NeIII3967'] = r'\regular{[NeIII]3967}'
labels['HI4340'] = r'\regular{H\gamma}'
labels['HI4861'] = r'\regular{H\beta}'
labels['OIII5007'] = r'\regular{[OIII]5007}'
labels['OIII4959'] = r'\regular{[OIII]4959}'
labels['HI6563'] = r'\regular{H\alpha}'
labels['NII6583'] = r'\regular{[NII]6584}'
labels['SII6731-SII6716'] = r'\regular{[SII]6717,6731}'



fig, axes = plt.subplots(2, 6, figsize=(7, 2.5), sharey = True, sharex = True)

left  = 0.1
right = 0.8
bottom = 0.15
top = 0.95
wspace = 0.02   # the amount of width reserved for blank space between subplots
hspace = 0.02

fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

grid = pickle.load(open('../../BuildGrid/Z_U_nH/grids/'+SPS+'/'+IMF+'/lines_SSP1.p','rb'), encoding='latin1')


# grid_alt = pickle.load(open('../../BuildGrid/Z/grids/BPASSv2.1.binary/ModSalpeter_100/lines.p','r'))
# grid_alt2 = pickle.load(open('../../BuildGrid/Z/grids/P2/ModSalpeter_100/lines.p','r'))

# print grid['log10Z']


inH_def = (np.abs(grid['log10n_H'] - 2.5)).argmin()
iU_def = (np.abs(grid['log10U_S_0'] - -2.0)).argmin()


print(grid['log10U_S_0'], iU_def)
print(grid['log10n_H'], inH_def)


inHs = [(np.abs(grid['log10n_H'] - log10n_H)).argmin() for log10n_H in [1., 2.,2.5, 3., 4.]]      
iUs = [(np.abs(grid['log10U_S_0'] - log10U_S_0)).argmin() for log10U_S_0 in [-1.,-2.,-3.]]      



axes = axes.flatten()



for il, line in enumerate(lines):

    c = cm.magma(float(il)/float(len(lines)))

    label = '-'.join(line)

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
    
        axes[il].plot(grid['log10Z'], (L-L_def), c=c, ls='-', lw=lw, alpha = 0.5)

    for iU, ls in zip(iUs,[':','-','--']):
    
        L = []
        for iZ, Z in enumerate(grid['Z']):
            line_luminosity = 0.0
            for l in line:
                l = l.encode('UTF-8')
                line_luminosity += grid[l]['luminosity'][iZ, iU, inH_def]
            L.append(line_luminosity)        
        
        L = np.log10(np.array(L))
    
        axes[il].plot(grid['log10Z'], (L-L_def), c=c, ls=ls, lw=1)

    


    # axes[il].set_yscale('log')
    
    axes[il].set_ylim([-1.1,1.1])
    axes[il].set_xlim([-5.,-1.5])
    axes[il].text(-3.25, 0.9, r'$\rm'+labels[label]+'$', ha = 'center', fontsize=5, color=c)











dummies = []
labels = []
# dummies += [ax.plot([], [], ls='-', c='k', alpha=0.5, lw = lw)[0] for lw in [1,2]]        

labels += [r'$\rm\log_{10}(U_{S,0})='+str(x)+'$' for x in [-1.,-2.,-3.]]
dummies += [axes[-1].plot([], [], ls=ls, c='k', alpha=0.5, lw = 1)[0] for ls in [':','-','--']]        

labels += [r'$\rm\log_{10}(n_{H})='+str(x)+'$' for x in [1., 2.,2.5, 3., 4.]]
dummies += [axes[-1].plot([], [], ls='-', c='k', alpha=0.2, lw = lw)[0] for lw in [0.5, 0.75, 1, 1.5, 2.0]]        


axes[-1].legend(dummies, labels,loc = 'center left', bbox_to_anchor = (1.0, 1.1), fontsize = 6, handlelength=2.5)



fig.text(0.45, 0.04, r'$\rm\log_{10}(Z)$', ha='center', fontsize=6)
fig.text(0.04, 0.55, r'$\rm\log_{10}(L/L_{default})$', va='center', rotation='vertical', fontsize=6)


        
# ax.legend(loc = 'center left', fontsize=7, bbox_to_anchor = (1.0, 0.5))
# 
# 
# ax.set_xlabel(r'$\rm\log_{10}(Z)$')
# ax.set_ylabel(r'$\rm\log_{10}[(L/L_{H\beta})_{1\ Myr}]$')


figname = 'figures/line_luminosity_effects.pdf'
print(figname)
fig.savefig(figname)

# ---- define plot



